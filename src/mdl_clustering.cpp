#include "mdl_clustering.h"
#include "distance_oracle.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <vector>

namespace clnj {

TreeDistCache precompute_tree_distances(
    const Adjacency& adjacency,
    const EdgeWeights& edge_weights,
    int n_observed
) {
    TreeDistCache cache;
    cache.n_obs = n_observed;
    cache.data.assign((size_t)n_observed * n_observed, 0.0);

    for (int src = 0; src < n_observed; ++src) {
        std::vector<double> dist(adjacency.size() + n_observed, -1.0);
        dist[src] = 0.0;

        std::queue<int> bfs;
        bfs.push(src);

        while (!bfs.empty()) {
            int node = bfs.front(); bfs.pop();

            auto it = adjacency.find(node);
            if (it == adjacency.end()) continue;

            for (int neighbor : it->second) {
                if (neighbor < (int)dist.size() && dist[neighbor] >= 0.0) continue;
                if (neighbor >= (int)dist.size()) dist.resize(neighbor + 1, -1.0);
                if (dist[neighbor] >= 0.0) continue;

                auto ek = normalize_edge(node, neighbor);
                double w = 0.0;
                auto ew_it = edge_weights.find(ek);
                if (ew_it != edge_weights.end())
                    w = std::max(0.0, ew_it->second);

                dist[neighbor] = dist[node] + w;
                bfs.push(neighbor);
            }
        }

        for (int dst = 0; dst < n_observed; ++dst) {
            double d = (dst < (int)dist.size() && dist[dst] >= 0.0) ? dist[dst] : 1e30;
            cache.data[(size_t)src * n_observed + dst] = d;
        }
    }

    return cache;
}

static std::set<int> collect_observed_in_subtree(
    int start, int blocked_node,
    const Adjacency& adj,
    int n_obs
) {
    std::set<int> result;
    std::queue<int> q;
    std::unordered_set<int> visited;
    q.push(start);
    visited.insert(start);
    visited.insert(blocked_node);

    while (!q.empty()) {
        int node = q.front(); q.pop();
        if (node < n_obs) result.insert(node);
        auto it = adj.find(node);
        if (it == adj.end()) continue;
        for (int nb : it->second) {
            if (!visited.count(nb)) {
                visited.insert(nb);
                q.push(nb);
            }
        }
    }
    return result;
}

static int find_central_node(
    const std::set<int>& members,
    const TreeDistCache& dc
) {
    if (members.size() <= 1) return members.empty() ? -1 : *members.begin();

    int best_node = -1;
    double best_sum_sq = std::numeric_limits<double>::max();

    for (int candidate : members) {
        double sum_sq = 0.0;
        for (int other : members) {
            if (candidate == other) continue;
            double d = dc(candidate, other);
            sum_sq += d * d;
        }
        if (sum_sq < best_sum_sq) {
            best_sum_sq = sum_sq;
            best_node = candidate;
        }
    }
    return best_node;
}

static double compute_DL(
    const std::set<int>& members,
    int center,
    double mj,
    const TreeDistCache& dc
) {
    int n = (int)members.size();
    if (n <= 1) return 0.0;

    double variance = 0.0;
    for (int m : members) {
        double d = dc(m, center);
        variance += d * d;
    }
    variance /= n;

    double log_var = std::log(std::max(variance, 1e-300));
    return n * log_var + mj * std::log((double)n) + n * (std::log(2.0 * M_PI) + 1.0);
}

static double compute_DL_split(
    const std::vector<std::set<int>>& parts,
    const std::vector<int>& centers,
    double mj,
    int n_total,
    const TreeDistCache& dc
) {
    double total = 0.0;
    for (size_t k = 0; k < parts.size(); ++k) {
        int nk = (int)parts[k].size();
        if (nk <= 1) continue;
        double variance = 0.0;
        for (int m : parts[k]) {
            double d = dc(m, centers[k]);
            variance += d * d;
        }
        variance /= nk;
        total += nk * std::log(std::max(variance, 1e-300));
    }
    total += (double)parts.size() * mj * std::log((double)n_total);
    total += n_total * (std::log(2.0 * M_PI) + 1.0);
    return total;
}

static std::vector<std::set<int>> split_at_node(
    int split_node,
    const std::set<int>& members,
    const Adjacency& adjacency,
    int n_obs
) {
    std::vector<std::set<int>> parts;
    auto it = adjacency.find(split_node);
    if (it == adjacency.end()) return parts;

    for (int nb : it->second) {
        auto sub = collect_observed_in_subtree(nb, split_node, adjacency, n_obs);
        std::set<int> filtered;
        for (int node : sub)
            if (members.count(node)) filtered.insert(node);
        if (!filtered.empty()) parts.push_back(std::move(filtered));
    }
    return parts;
}

static int find_best_split_node(
    const std::set<int>& members,
    const Adjacency& adjacency,
    const EdgeWeights& edge_weights,
    int n_obs,
    const TreeDistCache& dc
) {
    std::set<int> all_nodes_in_subtree;
    for (int m : members) all_nodes_in_subtree.insert(m);

    std::queue<int> bfs;
    std::unordered_set<int> visited;
    for (int m : members) { bfs.push(m); visited.insert(m); }
    while (!bfs.empty()) {
        int node = bfs.front(); bfs.pop();
        all_nodes_in_subtree.insert(node);
        auto it = adjacency.find(node);
        if (it == adjacency.end()) continue;
        for (int nb : it->second) {
            if (visited.count(nb)) continue;
            bool has_member_beyond = false;
            auto sub = collect_observed_in_subtree(nb, node, adjacency, n_obs);
            for (int s : sub)
                if (members.count(s)) { has_member_beyond = true; break; }
            if (has_member_beyond) {
                visited.insert(nb);
                bfs.push(nb);
            }
        }
    }

    int best_node = -1;
    double best_sum_sq = std::numeric_limits<double>::max();

    for (int candidate : all_nodes_in_subtree) {
        auto it = adjacency.find(candidate);
        if (it == adjacency.end()) continue;
        if (it->second.size() < 2) continue;

        double sum_sq = 0.0;
        if (candidate < n_obs) {
            for (int m : members) {
                double d = dc(candidate, m);
                sum_sq += d * d;
            }
        } else {
            for (int m : members) {
                double d = get_tree_distance(candidate, m, adjacency, edge_weights);
                sum_sq += d * d;
            }
        }
        if (sum_sq < best_sum_sq) {
            best_sum_sq = sum_sq;
            best_node = candidate;
        }
    }

    if (best_node == -1 && !members.empty())
        best_node = *members.begin();

    return best_node;
}

struct ClusterBuildCtx {
    const Adjacency& adjacency;
    const EdgeWeights& edge_weights;
    int n_obs;
    double mj;
    std::vector<Cluster>& clusters;
    int next_cluster_id = 0;
    const TreeDistCache& dc;
};

static int recursive_cluster(
    const std::set<int>& members,
    int parent_id,
    ClusterBuildCtx& ctx
) {
    if (members.size() <= 1) {
        Cluster c;
        c.id = ctx.next_cluster_id++;
        c.observed_members = members;
        c.center_node = members.empty() ? -1 : *members.begin();
        c.description_length = 0.0;
        c.parent_cluster_id = parent_id;
        ctx.clusters.push_back(std::move(c));
        return c.id;
    }

    int center = find_central_node(members, ctx.dc);
    double dl_unsplit = compute_DL(members, center, ctx.mj, ctx.dc);

    int split_node = find_best_split_node(members, ctx.adjacency, ctx.edge_weights,
                                          ctx.n_obs, ctx.dc);
    auto parts = split_at_node(split_node, members, ctx.adjacency, ctx.n_obs);

    if (parts.size() < 2) {
        Cluster c;
        c.id = ctx.next_cluster_id++;
        c.observed_members = members;
        c.center_node = center;
        c.description_length = dl_unsplit;
        c.parent_cluster_id = parent_id;
        ctx.clusters.push_back(std::move(c));
        return c.id;
    }

    std::vector<int> part_centers;
    for (auto& p : parts) {
        int pc = find_central_node(p, ctx.dc);
        part_centers.push_back(pc);
    }

    double dl_split = compute_DL_split(parts, part_centers, ctx.mj,
                                        (int)members.size(), ctx.dc);

    if (dl_split >= dl_unsplit) {
        Cluster c;
        c.id = ctx.next_cluster_id++;
        c.observed_members = members;
        c.center_node = center;
        c.description_length = dl_unsplit;
        c.parent_cluster_id = parent_id;
        ctx.clusters.push_back(std::move(c));
        return c.id;
    }

    Cluster parent_cluster;
    parent_cluster.id = ctx.next_cluster_id++;
    parent_cluster.observed_members = members;
    parent_cluster.center_node = center;
    parent_cluster.description_length = dl_split;
    parent_cluster.parent_cluster_id = parent_id;
    int parent_cid = parent_cluster.id;
    ctx.clusters.push_back(std::move(parent_cluster));

    for (auto& p : parts) {
        int child_id = recursive_cluster(p, parent_cid, ctx);
        for (auto& cl : ctx.clusters) {
            if (cl.id == parent_cid) {
                cl.child_cluster_ids.push_back(child_id);
                break;
            }
        }
    }

    return parent_cid;
}

static void merge_small_clusters(
    ClusterTree& ct,
    int merge_threshold,
    const TreeDistCache& dc
) {
    std::vector<int> leaf_cluster_ids;
    for (auto& c : ct.clusters)
        if (c.child_cluster_ids.empty())
            leaf_cluster_ids.push_back(c.id);

    for (int cid : leaf_cluster_ids) {
        Cluster* small = nullptr;
        for (auto& c : ct.clusters)
            if (c.id == cid) { small = &c; break; }
        if (!small || (int)small->observed_members.size() >= merge_threshold) continue;
        if (small->center_node < 0) continue;

        std::set<int> tried_ids;
        bool merged = false;

        for (int attempt = 0; attempt < 10 && !merged; ++attempt) {
            double best_dist = std::numeric_limits<double>::max();
            int best_merge_id = -1;

            for (auto& candidate : ct.clusters) {
                if (candidate.id == cid) continue;
                if (!candidate.child_cluster_ids.empty()) continue;
                if (candidate.center_node < 0) continue;
                if (candidate.observed_members.empty()) continue;
                if (tried_ids.count(candidate.id)) continue;

                double d = dc(small->center_node, candidate.center_node);
                if (d < best_dist) {
                    best_dist = d;
                    best_merge_id = candidate.id;
                }
            }
            if (best_merge_id < 0) break;
            tried_ids.insert(best_merge_id);

            Cluster* target = nullptr;
            for (auto& c : ct.clusters)
                if (c.id == best_merge_id) { target = &c; break; }
            if (!target) break;

            std::set<int> merged_members = small->observed_members;
            merged_members.insert(target->observed_members.begin(),
                                  target->observed_members.end());
            int merged_center = find_central_node(merged_members, dc);
            double dl_merged = compute_DL(merged_members, merged_center, ct.mj, dc);
            double dl_separate = small->description_length + target->description_length;

            if (dl_merged < dl_separate) {
                target->observed_members = merged_members;
                target->center_node = merged_center;
                target->description_length = dl_merged;
                small->observed_members.clear();
                small->center_node = -1;
                merged = true;
            }
        }
    }

    ct.clusters.erase(
        std::remove_if(ct.clusters.begin(), ct.clusters.end(),
                        [](const Cluster& c) { return c.observed_members.empty(); }),
        ct.clusters.end()
    );
}

ClusterTree mdl_cluster_tree(
    const Adjacency& adjacency,
    const EdgeWeights& edge_weights,
    const std::unordered_map<int, int>& hidden_info,
    int n_observed,
    double mj,
    int merge_threshold,
    bool verbose
) {
    (void)hidden_info;
    ClusterTree ct;
    ct.mj = mj;

    std::set<int> all_observed;
    for (int i = 0; i < n_observed; ++i)
        all_observed.insert(i);

    if (verbose)
        std::cout << "  MDL Clustering: " << n_observed << " observed nodes, mj=" << mj << "\n";

    auto t_cache_start = std::chrono::high_resolution_clock::now();
    TreeDistCache dc = precompute_tree_distances(adjacency, edge_weights, n_observed);
    auto t_cache_end = std::chrono::high_resolution_clock::now();
    double cache_ms = std::chrono::duration<double, std::milli>(t_cache_end - t_cache_start).count();
    if (verbose)
        std::cout << "  Pre-computed " << n_observed << "x" << n_observed
                  << " tree distance cache in " << cache_ms << " ms\n";

    ClusterBuildCtx ctx{adjacency, edge_weights, n_observed, mj, ct.clusters, 0, dc};

    recursive_cluster(all_observed, -1, ctx);

    if (verbose) {
        int n_leaf = 0;
        for (auto& c : ct.clusters)
            if (c.child_cluster_ids.empty()) ++n_leaf;
        std::cout << "  After recursive split: " << ct.clusters.size()
                  << " clusters (" << n_leaf << " leaf clusters)\n";
    }

    merge_small_clusters(ct, merge_threshold, dc);

    ct.node_to_cluster.clear();
    for (auto& c : ct.clusters) {
        if (c.child_cluster_ids.empty()) {
            for (int m : c.observed_members)
                ct.node_to_cluster[m] = c.id;
        }
    }

    if (verbose) {
        int n_leaf = 0;
        int min_size = n_observed, max_size = 0;
        for (auto& c : ct.clusters) {
            if (c.child_cluster_ids.empty()) {
                ++n_leaf;
                int sz = (int)c.observed_members.size();
                min_size = std::min(min_size, sz);
                max_size = std::max(max_size, sz);
            }
        }
        std::cout << "  After merging small clusters: " << ct.clusters.size()
                  << " total, " << n_leaf << " leaf clusters\n";
        if (n_leaf > 0)
            std::cout << "  Cluster sizes: min=" << min_size << ", max=" << max_size << "\n";
    }

    return ct;
}

}  // namespace clnj
