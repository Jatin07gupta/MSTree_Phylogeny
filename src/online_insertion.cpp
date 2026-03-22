#include "online_insertion.h"
#include "distance.h"
#include "distance_oracle.h"
#include "mdl_clustering.h"
#include "nj.h"
#include "mfnj.h"
#include "bionj.h"
#include "one_hot.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

namespace clnj {

using LocalFn = std::function<LocalResult(
    const std::vector<int>&, DistanceOracle&, int, int, double,
    NegCounts*, std::vector<ZeroEdgeEntry>*)>;

static LocalFn get_local_fn(TreeAlgo algo) {
    switch (algo) {
        case TreeAlgo::BIONJ: return bionj_local;
        case TreeAlgo::MFNJ:  return mfnj_local;
        default:              return nj_local;
    }
}

static void remap_hidden_nodes(TreeState& state, int shift_from, int shift_amount) {
    if (shift_amount == 0) return;
    auto remap = [&](int id) -> int {
        return (id >= shift_from) ? id + shift_amount : id;
    };

    Adjacency new_adj;
    for (const auto& [node, nbs] : state.adjacency) {
        int nn = remap(node);
        auto& new_nbs = new_adj[nn];
        for (int nb : nbs) new_nbs.insert(remap(nb));
    }
    state.adjacency = std::move(new_adj);

    EdgeWeights new_ew;
    for (const auto& [ek, w] : state.edge_weights) {
        new_ew[normalize_edge(remap(ek.first), remap(ek.second))] = w;
    }
    state.edge_weights = std::move(new_ew);

    std::unordered_map<int, int> new_hi;
    for (const auto& [hid, val] : state.hidden_info)
        new_hi[remap(hid)] = val;
    state.hidden_info = std::move(new_hi);

    state.next_hidden_id += shift_amount;
}

static std::set<int> find_steiner_tree(
    const std::set<int>& members,
    const Adjacency& adjacency
) {
    if (members.size() <= 1) return members;

    std::unordered_map<int, std::unordered_set<int>> adj;
    for (const auto& [node, nbs] : adjacency) adj[node] = nbs;

    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<int> to_remove;
        for (const auto& [node, nbs] : adj) {
            if (nbs.size() <= 1 && !members.count(node))
                to_remove.push_back(node);
        }
        for (int node : to_remove) {
            for (int nb : adj[node]) adj[nb].erase(node);
            adj.erase(node);
            changed = true;
        }
    }

    std::set<int> steiner;
    for (const auto& [node, _] : adj) steiner.insert(node);
    return steiner;
}

struct BoundaryEdge {
    int steiner_node;
    int external_node;
    double weight;
};

static std::vector<BoundaryEdge> find_boundary_edges(
    const std::set<int>& steiner_nodes,
    const Adjacency& adjacency,
    const EdgeWeights& edge_weights
) {
    std::vector<BoundaryEdge> result;
    for (int node : steiner_nodes) {
        auto it = adjacency.find(node);
        if (it == adjacency.end()) continue;
        for (int nb : it->second) {
            if (!steiner_nodes.count(nb)) {
                auto ek = normalize_edge(node, nb);
                double w = MIN_BRANCH;
                auto ew_it = edge_weights.find(ek);
                if (ew_it != edge_weights.end()) w = ew_it->second;
                result.push_back({node, nb, w});
            }
        }
    }
    return result;
}

static void remove_steiner_tree(
    const std::set<int>& steiner_nodes,
    const std::vector<BoundaryEdge>& boundary_edges,
    Adjacency& adjacency,
    EdgeWeights& edge_weights,
    std::unordered_map<int, int>& hidden_info,
    int n_observed
) {
    for (const auto& be : boundary_edges) {
        auto ek = normalize_edge(be.steiner_node, be.external_node);
        edge_weights.erase(ek);
        adjacency[be.steiner_node].erase(be.external_node);
        adjacency[be.external_node].erase(be.steiner_node);
    }

    std::vector<std::pair<int,int>> internal_edges;
    for (int node : steiner_nodes) {
        auto it = adjacency.find(node);
        if (it == adjacency.end()) continue;
        for (int nb : it->second) {
            if (steiner_nodes.count(nb) && node < nb)
                internal_edges.push_back({node, nb});
        }
    }
    for (const auto& [a, b] : internal_edges) {
        adjacency[a].erase(b);
        adjacency[b].erase(a);
        edge_weights.erase(normalize_edge(a, b));
    }

    for (int node : steiner_nodes) {
        if (node >= n_observed) {
            auto it = adjacency.find(node);
            if (it != adjacency.end() && it->second.empty()) {
                adjacency.erase(node);
                hidden_info.erase(node);
            }
        }
    }
}

static void install_nj_and_reconnect(
    LocalResult& local_result,
    const std::vector<BoundaryEdge>& boundary_edges,
    Adjacency& adjacency,
    EdgeWeights& edge_weights,
    std::unordered_map<int, int>& hidden_info,
    int& next_hidden_id,
    bool verbose
) {
    for (int h : local_result.new_hidden_ids) {
        hidden_info[h] = 0;
        adjacency[h] = {};
    }

    if (boundary_edges.empty()) {
        for (const auto& [n1, n2, w] : local_result.edges) {
            double d = std::max(MIN_BRANCH, w);
            adjacency[n1].insert(n2);
            adjacency[n2].insert(n1);
            edge_weights[normalize_edge(n1, n2)] = d;
        }
        return;
    }

    if (local_result.edges.empty()) {
        if (verbose)
            std::cerr << "  WARNING: NJ produced no edges for cluster rebuild\n";
        return;
    }

    auto& final_edge = local_result.edges.back();
    int fe_u = std::get<0>(final_edge);
    int fe_v = std::get<1>(final_edge);
    double fe_w = std::get<2>(final_edge);

    for (size_t i = 0; i + 1 < local_result.edges.size(); ++i) {
        const auto& [n1, n2, w] = local_result.edges[i];
        double d = std::max(MIN_BRANCH, w);
        adjacency[n1].insert(n2);
        adjacency[n2].insert(n1);
        edge_weights[normalize_edge(n1, n2)] = d;
    }

    int h_connector = next_hidden_id++;
    hidden_info[h_connector] = 0;
    adjacency[h_connector] = {};

    double w_u = std::max(fe_w / 2.0, MIN_BRANCH);
    double w_v = std::max(fe_w - w_u, MIN_BRANCH);

    adjacency[fe_u].insert(h_connector);
    adjacency[h_connector].insert(fe_u);
    edge_weights[normalize_edge(fe_u, h_connector)] = w_u;

    adjacency[h_connector].insert(fe_v);
    adjacency[fe_v].insert(h_connector);
    edge_weights[normalize_edge(h_connector, fe_v)] = w_v;

    const auto& primary = boundary_edges[0];
    adjacency[h_connector].insert(primary.external_node);
    adjacency[primary.external_node].insert(h_connector);
    edge_weights[normalize_edge(h_connector, primary.external_node)] = primary.weight;

    for (size_t i = 1; i < boundary_edges.size(); ++i) {
        const auto& be = boundary_edges[i];
        int attach_node = fe_u;
        adjacency[attach_node].insert(be.external_node);
        adjacency[be.external_node].insert(attach_node);
        edge_weights[normalize_edge(attach_node, be.external_node)] = be.weight;
        if (verbose)
            std::cout << "    Extra boundary: " << attach_node
                      << " -> " << be.external_node << "\n";
    }
}

static int find_cluster_hierarchical(
    const ClusterTree& ct,
    const VectorXd& new_dist_row,
    bool verbose
) {
    std::unordered_map<int, const Cluster*> id_map;
    for (const auto& c : ct.clusters)
        id_map[c.id] = &c;

    std::vector<const Cluster*> roots;
    for (const auto& c : ct.clusters)
        if (c.parent_cluster_id == -1)
            roots.push_back(&c);

    if (roots.empty()) return -1;

    int best_leaf = -1;
    double best_dist = std::numeric_limits<double>::max();

    for (const Cluster* root : roots) {
        const Cluster* current = root;

        if (verbose)
            (void)current;

        while (current) {
            if (current->child_cluster_ids.empty()) {
                if (current->center_node >= 0 &&
                    current->center_node < (int)new_dist_row.size()) {
                    double d = new_dist_row(current->center_node);
                    if (d < best_dist) {
                        best_dist = d;
                        best_leaf = current->id;
                    }
                }
                break;
            }

            const Cluster* nearest_child = nullptr;
            double nearest_dist = std::numeric_limits<double>::max();

            for (int child_id : current->child_cluster_ids) {
                auto it = id_map.find(child_id);
                if (it == id_map.end()) continue;
                const Cluster* child = it->second;
                if (child->center_node < 0 ||
                    child->center_node >= (int)new_dist_row.size()) continue;
                if (child->observed_members.empty()) continue;
                double d = new_dist_row(child->center_node);
                if (d < nearest_dist) {
                    nearest_dist = d;
                    nearest_child = child;
                }
            }

            if (!nearest_child) {
                if (current->center_node >= 0 &&
                    current->center_node < (int)new_dist_row.size()) {
                    double d = new_dist_row(current->center_node);
                    if (d < best_dist) {
                        best_dist = d;
                        best_leaf = current->id;
                    }
                }
                break;
            }

            current = nearest_child;
        }
    }

    if (best_leaf < 0) {
        for (const auto& c : ct.clusters) {
            if (!c.child_cluster_ids.empty()) continue;
            if (c.center_node < 0 ||
                c.center_node >= (int)new_dist_row.size()) continue;
            double d = new_dist_row(c.center_node);
            if (d < best_dist) {
                best_dist = d;
                best_leaf = c.id;
            }
        }
        if (verbose && best_leaf >= 0)
            std::cout << "    (fallback to flat scan -> cluster "
                      << best_leaf << ")\n";
    }

    return best_leaf;
}

static bool check_tree_invariant(const Adjacency& adjacency, bool verbose) {
    int n_adj = 0, e_count = 0;
    for (const auto& [v, nbs] : adjacency) { ++n_adj; e_count += (int)nbs.size(); }
    e_count /= 2;
    bool ok = (n_adj == 0) || (e_count == n_adj - 1);
    if (!ok && verbose)
        std::cout << "  WARNING: Tree invariant violated: |E|="
                  << e_count << " != |V|-1=" << n_adj - 1 << "\n";
    else if (verbose)
        std::cout << "  Tree invariant: PASSED (|V|=" << n_adj
                  << ", |E|=" << e_count << ")\n";
    return ok;
}

InsertionResult insert_taxon(
    TreeState& state,
    const std::string& new_name,
    const std::string& new_sequence,
    bool verbose
) {
    auto results = insert_batch(state, {new_name}, {new_sequence}, verbose);
    return results.empty()
        ? InsertionResult{-1, -1, 0, false, "insert_batch returned empty"}
        : results[0];
}

std::vector<InsertionResult> insert_batch(
    TreeState& state,
    const std::vector<std::string>& new_names,
    const std::vector<std::string>& new_sequences,
    bool verbose
) {
    int batch_size = (int)new_names.size();
    if (batch_size == 0 || new_names.size() != new_sequences.size()) {
        if (batch_size != 0)
            std::cerr << "ERROR: Mismatched names/sequences count\n";
        return {};
    }
    if (state.cluster_tree.clusters.empty()) {
        std::cerr << "ERROR: No clusters in tree state\n";
        return {InsertionResult{-1, -1, 0, false, "No clusters"}};
    }

    if (verbose)
        std::cout << "\n=== NJ-Rebuild Batch Insertion of " << batch_size << " taxa (CPU) ===\n";

    auto t_total_start = std::chrono::high_resolution_clock::now();

    std::set<int> center_ids;
    for (const auto& c : state.cluster_tree.clusters)
        if (c.center_node >= 0)
            center_ids.insert(c.center_node);
    std::vector<int> center_list(center_ids.begin(), center_ids.end());

    int old_n = state.n_observed;
    int m = state.aln_len;

    bool using_cache = state.has_cached_encoding();
    OneHotData oh;
    if (!using_cache) {
        oh = one_hot_encode(state.sequences);
        if (verbose)
            std::cout << "  (Cache miss: re-encoded " << old_n << " sequences)\n";
    }

    std::vector<int8_t>  all_new_nuc(batch_size * m);
    std::vector<uint8_t> all_new_valid(batch_size * m);
    for (int i = 0; i < batch_size; ++i)
        encode_sequence_flat(new_sequences[i], m,
                             all_new_nuc.data() + i * m,
                             all_new_valid.data() + i * m);

    std::vector<std::unordered_map<int, double>> center_dists(batch_size);

    for (int i = 0; i < batch_size; ++i) {
        if (using_cache) {
            center_dists[i] = compute_distance_to_subset_cached(
                state.cached_nuc_idx.data(), state.cached_valid.data(),
                old_n, m,
                all_new_nuc.data() + i * m, all_new_valid.data() + i * m,
                center_list, state.model, state.gamma_alpha,
                state.gA, state.gC, state.gG, state.gT, state.min_shared_sites
            );
        } else {
            center_dists[i] = compute_distance_to_subset(
                oh, new_sequences[i], center_list, state.model, state.gamma_alpha,
                state.gA, state.gC, state.gG, state.gT, state.min_shared_sites
            );
        }
    }

    std::vector<VectorXd> sparse_rows(batch_size);
    for (int i = 0; i < batch_size; ++i) {
        sparse_rows[i] = VectorXd::Constant(old_n, 1e30);
        for (const auto& [idx, d] : center_dists[i])
            if (idx < old_n) sparse_rows[i](idx) = d;
    }

    if (verbose)
        std::cout << "  Phase 1: Computed distances to " << center_list.size()
                  << " cluster centers" << (using_cache ? " [CACHED]" : " [re-encoded]") << "\n";

    std::vector<int> assignments(batch_size, -1);
    for (int i = 0; i < batch_size; ++i) {
        if (verbose)
            std::cout << "  Assigning '" << new_names[i] << "':\n";
        assignments[i] = find_cluster_hierarchical(
            state.cluster_tree, sparse_rows[i], verbose
        );
    }

    std::map<int, std::vector<int>> cluster_to_new_indices;
    for (int i = 0; i < batch_size; ++i)
        if (assignments[i] >= 0)
            cluster_to_new_indices[assignments[i]].push_back(i);

    if (verbose) {
        std::cout << "  Cluster assignments:\n";
        for (const auto& [cid, idxs] : cluster_to_new_indices)
            std::cout << "    Cluster " << cid << ": " << idxs.size() << " new taxa\n";
    }

    remap_hidden_nodes(state, old_n, batch_size);

    int new_n = old_n + batch_size;
    MatrixXd new_D = MatrixXd::Zero(new_n, new_n);
    new_D.topLeftCorner(old_n, old_n) = state.D;

    int total_dist_pairs = 0;

    for (const auto& [cid, new_indices] : cluster_to_new_indices) {
        Cluster* cluster_ptr = nullptr;
        for (auto& c : state.cluster_tree.clusters)
            if (c.id == cid) { cluster_ptr = &c; break; }
        if (!cluster_ptr) continue;

        std::vector<int> old_members(cluster_ptr->observed_members.begin(),
                                     cluster_ptr->observed_members.end());

        for (int idx : new_indices) {
            std::unordered_map<int, double> dists;
            if (using_cache) {
                dists = compute_distance_to_subset_cached(
                    state.cached_nuc_idx.data(), state.cached_valid.data(),
                    old_n, m,
                    all_new_nuc.data() + idx * m, all_new_valid.data() + idx * m,
                    old_members, state.model, state.gamma_alpha,
                    state.gA, state.gC, state.gG, state.gT, state.min_shared_sites
                );
            } else {
                dists = compute_distance_to_subset(
                    oh, new_sequences[idx], old_members, state.model, state.gamma_alpha,
                    state.gA, state.gC, state.gG, state.gT, state.min_shared_sites
                );
            }
            int new_id = old_n + idx;
            for (const auto& [target, d] : dists) {
                new_D(target, new_id) = d;
                new_D(new_id, target) = d;
            }
            total_dist_pairs += (int)dists.size();
        }

        if (new_indices.size() > 1) {
            int local_n = (int)new_indices.size();
            for (int a = 0; a < local_n; ++a) {
                for (int b = a + 1; b < local_n; ++b) {
                    int idx_a = new_indices[a], idx_b = new_indices[b];
                    const int8_t*  nuc_a = all_new_nuc.data() + idx_a * m;
                    const uint8_t* val_a = all_new_valid.data() + idx_a * m;
                    const int8_t*  nuc_b = all_new_nuc.data() + idx_b * m;
                    const uint8_t* val_b = all_new_valid.data() + idx_b * m;

                    int n_valid = 0, n_diff = 0, n_ts_AG = 0, n_ts_CT = 0;
                    for (int s = 0; s < m; ++s) {
                        if (!val_a[s] || !val_b[s]) continue;
                        ++n_valid;
                        if (nuc_a[s] != nuc_b[s]) {
                            ++n_diff;
                            int sum = nuc_a[s] + nuc_b[s];
                            if (sum == 2) ++n_ts_AG;
                            else if (sum == 4) ++n_ts_CT;
                        }
                    }

                    double d = 0.0;
                    if (n_valid < std::max(state.min_shared_sites, 1)) {
                        d = MAX_DIST;
                    } else {
                        double nv = (double)n_valid;
                        int n_tv = n_diff - n_ts_AG - n_ts_CT;
                        if (state.model == DistModel::JC69)
                            d = jc69_correct((double)n_diff / nv, state.gamma_alpha);
                        else if (state.model == DistModel::K2P)
                            d = k2p_correct((double)(n_ts_AG + n_ts_CT) / nv,
                                            (double)n_tv / nv, state.gamma_alpha);
                        else if (state.model == DistModel::TN93)
                            d = tn93_correct((double)n_ts_AG / nv, (double)n_ts_CT / nv,
                                             (double)n_tv / nv,
                                             state.gA, state.gC, state.gG, state.gT,
                                             state.gamma_alpha);
                    }
                    d = std::max(d, 0.0);

                    int new_id_a = old_n + idx_a;
                    int new_id_b = old_n + idx_b;
                    new_D(new_id_a, new_id_b) = d;
                    new_D(new_id_b, new_id_a) = d;
                    ++total_dist_pairs;
                }
            }
        }
    }

    state.D = std::move(new_D);

    for (int i = 0; i < batch_size; ++i) {
        state.names.push_back(new_names[i]);
        state.sequences.push_back(new_sequences[i]);
        int new_id = old_n + i;
        state.adjacency[new_id] = {};
    }
    state.n_observed = new_n;

    if (using_cache) {
        state.cached_nuc_idx.resize((size_t)new_n * m);
        state.cached_valid.resize((size_t)new_n * m);
        std::copy(all_new_nuc.begin(), all_new_nuc.end(),
                  state.cached_nuc_idx.begin() + (size_t)old_n * m);
        std::copy(all_new_valid.begin(), all_new_valid.end(),
                  state.cached_valid.begin() + (size_t)old_n * m);
    }

    std::vector<InsertionResult> results(batch_size);
    LocalFn local_fn = get_local_fn(state.tree_algo);

    for (const auto& [cid, new_indices] : cluster_to_new_indices) {
        Cluster* cluster_ptr = nullptr;
        for (auto& c : state.cluster_tree.clusters)
            if (c.id == cid) { cluster_ptr = &c; break; }
        if (!cluster_ptr) {
            for (int idx : new_indices) {
                results[idx] = {old_n + idx, cid, 0, false, "Cluster not found"};
            }
            continue;
        }

        std::set<int> old_members = cluster_ptr->observed_members;
        std::set<int> all_members = old_members;
        for (int idx : new_indices) all_members.insert(old_n + idx);

        if (verbose) {
            std::cout << "\n  Rebuilding cluster " << cid << ": "
                      << old_members.size() << " old + "
                      << new_indices.size() << " new = "
                      << all_members.size() << " total members\n";
        }

        auto steiner = find_steiner_tree(old_members, state.adjacency);

        if (verbose) {
            int n_hidden_in_steiner = 0;
            for (int s : steiner)
                if (s >= new_n) ++n_hidden_in_steiner;
            std::cout << "    Steiner tree: " << steiner.size() << " nodes ("
                      << n_hidden_in_steiner << " hidden)\n";
        }

        auto boundary = find_boundary_edges(steiner, state.adjacency, state.edge_weights);
        if (verbose)
            std::cout << "    Boundary edges: " << boundary.size() << "\n";

        remove_steiner_tree(steiner, boundary, state.adjacency,
                           state.edge_weights, state.hidden_info, new_n);

        std::vector<int> member_list(all_members.begin(), all_members.end());
        std::sort(member_list.begin(), member_list.end());

        DistanceOracle oracle(state.D, new_n, state.adjacency, state.edge_weights);
        auto local_result = local_fn(
            member_list, oracle, new_n, state.next_hidden_id,
            MIN_BRANCH, nullptr, nullptr
        );

        for (int h : local_result.new_hidden_ids)
            state.next_hidden_id = std::max(state.next_hidden_id, h + 1);

        if (verbose) {
            std::cout << "    NJ produced " << local_result.edges.size()
                      << " edges, " << local_result.new_hidden_ids.size()
                      << " hidden nodes\n";
        }

        if (old_members.size() == 1 && boundary.empty()) {
            for (const auto& [n1, n2, w] : local_result.edges) {
                double d = std::max(MIN_BRANCH, w);
                state.adjacency[n1].insert(n2);
                state.adjacency[n2].insert(n1);
                state.edge_weights[normalize_edge(n1, n2)] = d;
            }
            for (int h : local_result.new_hidden_ids) {
                state.hidden_info[h] = 0;
                if (state.adjacency.find(h) == state.adjacency.end())
                    state.adjacency[h] = {};
            }
        } else {
            install_nj_and_reconnect(
                local_result, boundary,
                state.adjacency, state.edge_weights,
                state.hidden_info, state.next_hidden_id,
                verbose
            );
        }

        cluster_ptr->observed_members = all_members;
        for (int idx : new_indices) {
            int nid = old_n + idx;
            state.cluster_tree.node_to_cluster[nid] = cid;
        }

        for (int idx : new_indices) {
            int nid = old_n + idx;
            results[idx] = {
                nid, cid, (int)all_members.size(), true,
                "Inserted '" + new_names[idx] + "' via NJ rebuild in cluster "
                + std::to_string(cid)
            };
        }
    }

    for (int i = 0; i < batch_size; ++i) {
        if (assignments[i] < 0) {
            results[i] = {old_n + i, -1, 0, false,
                          "Could not assign to any cluster"};
        }
    }

    auto t_total_end = std::chrono::high_resolution_clock::now();
    double total_ms = std::chrono::duration<double, std::milli>(t_total_end - t_total_start).count();

    if (verbose) {
        check_tree_invariant(state.adjacency, verbose);
        int ok = 0;
        for (const auto& r : results) if (r.success) ++ok;
        std::cout << "\n=== Batch insertion complete: " << ok << "/"
                  << batch_size << " successful ===\n";
        std::cout << "  Total insertion: " << total_ms << " ms\n";
    }

    return results;
}

}  // namespace clnj
