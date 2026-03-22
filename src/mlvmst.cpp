#include "mlvmst.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

namespace clnj {

static std::vector<int> compute_delta_max_int(
    int n,
    const std::vector<std::set<int>>& F_C,
    const std::map<int, std::set<int>>& G_U_adj
) {
    // F_C sorted by size descending
    auto F_C_ge = F_C;
    std::sort(F_C_ge.begin(), F_C_ge.end(),
              [](const std::set<int>& a, const std::set<int>& b) {
                  return a.size() > b.size();
              });

    std::vector<int> delta(n, 0);
    for (int i = 0; i < n; ++i) {
        auto it = G_U_adj.find(i);
        if (it == G_U_adj.end()) continue;
        std::set<int> N_i = it->second;

        for (auto& C : F_C_ge) {
            if (N_i.empty()) break;
            if (C.count(i)) continue;
            bool has_common = false;
            for (int v : N_i)
                if (C.count(v)) { has_common = true; break; }
            if (!has_common) continue;

            ++delta[i];
            for (auto nit = N_i.begin(); nit != N_i.end(); ) {
                if (C.count(*nit)) nit = N_i.erase(nit);
                else ++nit;
            }
        }
    }
    return delta;
}

MlvmstResult build_mlvmst(
    const MatrixXd& D,
    const Algo2Result& algo2,
    const std::string& ordering,
    const std::string& mst_algorithm,
    bool verbose
) {
    int n = (int)D.rows();
    auto t0 = std::chrono::steady_clock::now();

    // Build G_U adjacency
    std::map<int, std::set<int>> G_U_adj;
    for (auto& [u, v] : algo2.G_U_edges) {
        G_U_adj[u].insert(v);
        G_U_adj[v].insert(u);
    }

    auto delta = compute_delta_max_int(n, algo2.F_C, G_U_adj);

    // Vertex ordering
    std::vector<int> vertex_order(n);
    std::iota(vertex_order.begin(), vertex_order.end(), 0);
    if (ordering == "increasing") {
        std::sort(vertex_order.begin(), vertex_order.end(),
                  [&](int a, int b) {
                      return delta[a] < delta[b] || (delta[a] == delta[b] && a < b);
                  });
    } else {
        std::sort(vertex_order.begin(), vertex_order.end(),
                  [&](int a, int b) {
                      return delta[a] > delta[b] || (delta[a] == delta[b] && a < b);
                  });
    }

    std::vector<int> rank(n);
    for (int i = 0; i < n; ++i) rank[vertex_order[i]] = i;

    // Collect edges, sort by (weight, min_rank, max_rank)
    struct MSTEdge {
        int u, v;
        double w;
        int min_rank, max_rank;
    };

    std::vector<MSTEdge> edges;
    edges.reserve(n * (n - 1) / 2);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            double w = D(i, j);
            if (!std::isfinite(w)) w = 1e30;
            int r1 = std::min(rank[i], rank[j]);
            int r2 = std::max(rank[i], rank[j]);
            edges.push_back({i, j, w, r1, r2});
        }

    std::sort(edges.begin(), edges.end(),
              [](const MSTEdge& a, const MSTEdge& b) {
                  if (a.w != b.w) return a.w < b.w;
                  if (a.min_rank != b.min_rank) return a.min_rank < b.min_rank;
                  return a.max_rank < b.max_rank;
              });

    // Shared DSU primitives for MST backends
    std::vector<int> parent(n);
    std::vector<int> rnk(n, 0);
    std::iota(parent.begin(), parent.end(), 0);

    auto find = [&](int x) -> int {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto unite = [&](int a, int b) -> bool {
        int ra = find(a), rb = find(b);
        if (ra == rb) return false;
        if (rnk[ra] < rnk[rb]) std::swap(ra, rb);
        parent[rb] = ra;
        if (rnk[ra] == rnk[rb]) ++rnk[ra];
        return true;
    };

    MlvmstResult result;
    result.adjacency.resize(n, n);
    result.adjacency.setZero();
    result.delta_max = delta;

    auto add_tree_edge = [&](int u, int v, int& edges_added) {
        result.adjacency(u, v) = 1;
        result.adjacency(v, u) = 1;
        ++edges_added;
    };

    int edges_added = 0;
    std::string mst_algo = mst_algorithm;
    std::transform(mst_algo.begin(), mst_algo.end(), mst_algo.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });

    if (mst_algo == "boruvka") {
        // Boruvka: in each round, every component chooses its lightest outgoing edge.
        int components = n;
        while (components > 1 && edges_added < n - 1) {
            std::vector<int> cheapest(n, -1);

            for (int idx = 0; idx < (int)edges.size(); ++idx) {
                const auto& e = edges[idx];
                if (e.w > 1e29) continue;
                int ru = find(e.u), rv = find(e.v);
                if (ru == rv) continue;

                auto better = [&](int a, int b) {
                    if (b < 0) return true;
                    const auto& ea = edges[a];
                    const auto& eb = edges[b];
                    if (ea.w != eb.w) return ea.w < eb.w;
                    if (ea.min_rank != eb.min_rank) return ea.min_rank < eb.min_rank;
                    if (ea.max_rank != eb.max_rank) return ea.max_rank < eb.max_rank;
                    if (ea.u != eb.u) return ea.u < eb.u;
                    return ea.v < eb.v;
                };

                if (better(idx, cheapest[ru])) cheapest[ru] = idx;
                if (better(idx, cheapest[rv])) cheapest[rv] = idx;
            }

            bool merged_any = false;
            for (int c = 0; c < n; ++c) {
                int idx = cheapest[c];
                if (idx < 0) continue;
                const auto& e = edges[idx];
                if (unite(e.u, e.v)) {
                    add_tree_edge(e.u, e.v, edges_added);
                    --components;
                    merged_any = true;
                    if (edges_added >= n - 1) break;
                }
            }
            if (!merged_any) break;
        }
    } else {
        // Default: Kruskal
        if (mst_algo != "kruskal" && verbose) {
            std::cout << "  WARNING: Unknown mst_algorithm='" << mst_algorithm
                      << "'. Falling back to kruskal.\n";
        }
        for (auto& e : edges) {
            if (e.w > 1e29) break;
            if (unite(e.u, e.v)) {
                add_tree_edge(e.u, e.v, edges_added);
                if (edges_added >= n - 1) break;
            }
        }
    }

    // Compute degree and leaf count
    int max_deg = 0, min_deg = n;
    result.leaf_count = 0;
    for (int i = 0; i < n; ++i) {
        int deg = 0;
        for (int j = 0; j < n; ++j)
            if (result.adjacency(i, j)) ++deg;
        if (deg == 1) ++result.leaf_count;
        if (deg > 0) {
            max_deg = std::max(max_deg, deg);
            min_deg = std::min(min_deg, deg);
        }
    }

    if (verbose) {
        auto t1 = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();
        std::cout << "  MLVMST built in " << elapsed << "s\n";
        std::cout << "  MST backend: " << mst_algo << "\n";
        std::cout << "  Edges: " << edges_added << ", Leaves: " << result.leaf_count << "\n";
        std::cout << "  Degree range: " << min_deg << "-" << max_deg << "\n";

        int dmin = *std::min_element(delta.begin(), delta.end());
        int dmax = *std::max_element(delta.begin(), delta.end());
        std::cout << "  Delta_max range: " << dmin << "-" << dmax << "\n";
    }

    return result;
}

}  // namespace clnj
