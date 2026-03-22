#include "distance_oracle.h"

#include <queue>

namespace clnj {

double get_tree_distance(int a, int b,
                         const Adjacency& adjacency,
                         const EdgeWeights& edge_weights) {
    if (a == b) return 0.0;

    std::unordered_set<int> visited;
    visited.insert(a);
    std::queue<std::pair<int, double>> bfs;
    bfs.push({a, 0.0});

    while (!bfs.empty()) {
        auto [node, dist_so_far] = bfs.front();
        bfs.pop();

        auto it = adjacency.find(node);
        if (it == adjacency.end()) continue;

        for (int neighbor : it->second) {
            if (visited.count(neighbor)) continue;
            visited.insert(neighbor);

            auto ek = normalize_edge(node, neighbor);
            double w = 0.0;
            auto ew_it = edge_weights.find(ek);
            if (ew_it != edge_weights.end())
                w = std::max(0.0, ew_it->second);
            else
                w = 1e30;  // inf

            double new_dist = dist_so_far + w;
            if (neighbor == b) return new_dist;
            bfs.push({neighbor, new_dist});
        }
    }
    return 1e30;  // unreachable
}

DistanceOracle::DistanceOracle(const MatrixXd& D_obs, int m,
                               Adjacency& adjacency, EdgeWeights& edge_weights)
    : D_obs_(D_obs), m_(m), adjacency_(adjacency), edge_weights_(edge_weights) {}

double DistanceOracle::dist(int i, int j) const {
    if (i == j) return 0.0;
    if (i < m_ && j < m_) return D_obs_(i, j);
    return std::max(0.0, get_tree_distance(i, j, adjacency_, edge_weights_));
}

double DistanceOracle::dist_bfs(int i, int j) const {
    if (i == j) return 0.0;
    return std::max(0.0, get_tree_distance(i, j, adjacency_, edge_weights_));
}

double DistanceOracle::dist_for_neighbourhood(int i, int j, bool has_hidden) const {
    if (i == j) return 0.0;
    if (has_hidden) return dist_bfs(i, j);
    return dist(i, j);
}

}  // namespace clnj
