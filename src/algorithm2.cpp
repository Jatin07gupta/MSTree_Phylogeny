#include "algorithm2.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace clnj {

// Disjoint Set Union with path compression and union by rank
class DSU {
public:
    std::vector<int> parent;
    std::vector<int> rank_;
    std::vector<std::set<int>> members;

    explicit DSU(int n) : parent(n), rank_(n, 0), members(n) {
        std::iota(parent.begin(), parent.end(), 0);
        for (int i = 0; i < n; ++i)
            members[i].insert(i);
    }

    int find(int x) {
        int r = x;
        while (parent[r] != r) r = parent[r];
        while (parent[x] != r) {
            int next = parent[x];
            parent[x] = r;
            x = next;
        }
        return r;
    }

    int unite(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return ra;
        if (rank_[ra] < rank_[rb]) std::swap(ra, rb);
        parent[rb] = ra;
        members[ra].insert(members[rb].begin(), members[rb].end());
        members[rb].clear();
        if (rank_[ra] == rank_[rb]) ++rank_[ra];
        return ra;
    }
};

Algo2Result construct_FC_and_GU(const MatrixXd& D, bool verbose) {
    int n = (int)D.rows();
    auto t0 = std::chrono::steady_clock::now();

    // Collect all upper-triangle pairs and sort by distance
    struct Edge {
        int u, v;
        double w;
    };
    int n_pairs = n * (n - 1) / 2;
    std::vector<Edge> edges;
    edges.reserve(n_pairs);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            double w = D(i, j);
            if (!std::isfinite(w)) w = 1e30;
            edges.push_back({i, j, w});
        }
    std::sort(edges.begin(), edges.end(),
              [](const Edge& a, const Edge& b) { return a.w < b.w; });

    DSU dsu(n);

    std::set<std::set<int>> FC_set;
    // Add full vertex set and all singletons
    std::set<int> full;
    for (int i = 0; i < n; ++i) {
        full.insert(i);
        FC_set.insert({i});
    }
    FC_set.insert(full);

    std::set<std::pair<int,int>> GU_edges;

    // Tier buffers (Algorithm 2 style):
    // - E_w stores edges in current weight tier
    // - V_w stores vertices touched by those edges
    std::vector<std::pair<int,int>> E_w;
    std::set<int> V_w;

    auto flush_tier = [&]() {
        // 1) Union all buffered edges
        for (auto& [u, v] : E_w) {
            int ru = dsu.find(u), rv = dsu.find(v);
            if (ru != rv) dsu.unite(u, v);
        }

        // 2) Add resulting components for touched vertices
        for (int v : V_w) {
            int r = dsu.find(v);
            auto& m = dsu.members[r];
            if (m.size() > 1) FC_set.insert(m);
        }

        // 3) Clear buffers
        E_w.clear();
        V_w.clear();
    };

    if (edges.empty()) {
        Algo2Result result;
        result.F_C.assign(FC_set.begin(), FC_set.end());
        std::sort(result.F_C.begin(), result.F_C.end(),
                  [](const std::set<int>& a, const std::set<int>& b) {
                      return a.size() < b.size();
                  });
        result.G_U_edges = GU_edges;
        return result;
    }

    double prev_w = edges.front().w;

    for (auto& e : edges) {
        if (!std::isfinite(e.w) || e.w > 1e29) break;

        // Tier boundary: flush previous tier before processing this edge
        if (e.w > prev_w + 1e-15) {
            flush_tier();
        }

        int ru = dsu.find(e.u), rv = dsu.find(e.v);
        if (ru != rv) {
            GU_edges.insert({std::min(e.u, e.v), std::max(e.u, e.v)});
            E_w.push_back({e.u, e.v});
            V_w.insert(e.u);
            V_w.insert(e.v);
        }
        prev_w = e.w;
    }

    // Final tier flush
    flush_tier();

    // Convert to sorted vector
    Algo2Result result;
    result.F_C.assign(FC_set.begin(), FC_set.end());
    std::sort(result.F_C.begin(), result.F_C.end(),
              [](const std::set<int>& a, const std::set<int>& b) {
                  return a.size() < b.size();
              });
    result.G_U_edges = GU_edges;

    if (verbose) {
        auto t1 = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();
        std::cout << "  Algorithm 2 done in " << elapsed << "s\n";
        std::cout << "  F_C: " << result.F_C.size() << " components, G_U: "
                  << result.G_U_edges.size() << " edges\n";
    }

    return result;
}

}  // namespace clnj
