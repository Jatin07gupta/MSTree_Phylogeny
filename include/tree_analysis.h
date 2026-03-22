#pragma once

#include "types.h"
#include <vector>
#include <string>

namespace clnj {

TreeStats analyze_tree(
    const Adjacency& adjacency,
    const EdgeWeights& edge_weights,
    const std::unordered_map<int, int>& hidden_info,
    int n_obs,
    const std::vector<std::string>& names,
    bool print_report = true
);

}  // namespace clnj
