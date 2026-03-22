#pragma once

#include "types.h"
#include <vector>

namespace clnj {

struct TreeDistCache {
    std::vector<double> data;
    int n_obs = 0;

    double operator()(int i, int j) const {
        return data[(size_t)i * n_obs + j];
    }

    bool empty() const { return data.empty(); }
};

TreeDistCache precompute_tree_distances(
    const Adjacency& adjacency,
    const EdgeWeights& edge_weights,
    int n_observed
);

ClusterTree mdl_cluster_tree(
    const Adjacency& adjacency,
    const EdgeWeights& edge_weights,
    const std::unordered_map<int, int>& hidden_info,
    int n_observed,
    double mj = 7.0,
    int merge_threshold = 100,
    bool verbose = true
);

}  // namespace clnj
