#include "tree_analysis.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <vector>

namespace clnj {

TreeStats analyze_tree(
    const Adjacency& adjacency,
    const EdgeWeights& edge_weights,
    const std::unordered_map<int, int>& hidden_info,
    int n_obs,
    const std::vector<std::string>& names,
    bool print_report
) {
    TreeStats stats;
    stats.n_obs = n_obs;
    stats.n_hidden = (int)hidden_info.size();
    stats.total_nodes = (int)adjacency.size();
    stats.total_edges = (int)edge_weights.size();
    stats.valid_tree = (stats.total_edges == stats.total_nodes - 1);

    if (edge_weights.empty()) {
        if (print_report) std::cout << "  No edges in tree.\n";
        return stats;
    }

    std::vector<double> weights;
    weights.reserve(edge_weights.size());
    for (auto& [ek, w] : edge_weights) weights.push_back(w);
    std::sort(weights.begin(), weights.end());

    stats.min_weight = weights.front();
    stats.max_weight = weights.back();
    double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    stats.mean_weight = sum / weights.size();
    if (weights.size() % 2 == 0)
        stats.median_weight = (weights[weights.size()/2 - 1] + weights[weights.size()/2]) / 2.0;
    else
        stats.median_weight = weights[weights.size()/2];

    double sq_sum = 0;
    for (double w : weights) sq_sum += (w - stats.mean_weight) * (w - stats.mean_weight);
    stats.std_weight = std::sqrt(sq_sum / weights.size());

    for (double w : weights) {
        if (w < -EPS) ++stats.negative_edges;
        else if (std::abs(w) < EPS) ++stats.zero_edges;
        else ++stats.positive_edges;
    }

    for (auto& [ek, w] : edge_weights) {
        auto [a, b] = ek;
        if (a < n_obs && b < n_obs) ++stats.obs_obs;
        else if (a < n_obs || b < n_obs) ++stats.obs_hidden;
        else ++stats.hidden_hidden;
    }

    stats.max_degree = 0;
    stats.leaves = 0;
    for (auto& [nd, neighbors] : adjacency) {
        int deg = (int)neighbors.size();
        if (deg == 1) ++stats.leaves;
        stats.max_degree = std::max(stats.max_degree, deg);
    }

    if (print_report) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "FINAL TREE ANALYSIS\n";
        std::cout << std::string(70, '=') << "\n";
        std::cout << "  Observed nodes:     " << stats.n_obs << "\n";
        std::cout << "  Hidden nodes:       " << stats.n_hidden << "\n";
        std::cout << "  Total nodes:        " << stats.total_nodes << "\n";
        std::cout << "  Total edges:        " << stats.total_edges << "\n";
        std::cout << "  Tree valid (|E|=|V|-1): "
                  << (stats.valid_tree ? "YES" : "NO") << "\n";

        std::cout << "\n  Edge weights:\n";
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "    Min:     " << stats.min_weight << "\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "    Max:     " << stats.max_weight << "\n";
        std::cout << "    Mean:    " << stats.mean_weight << "\n";
        std::cout << "    Median:  " << stats.median_weight << "\n";
        std::cout << "    Std:     " << stats.std_weight << "\n";
        std::cout << "    Positive: " << stats.positive_edges
                  << ",  Zero: " << stats.zero_edges
                  << ",  Negative: " << stats.negative_edges << "\n";

        std::cout << "\n  Edge types:\n";
        std::cout << "    Obs-Obs: " << stats.obs_obs
                  << ",  Obs-Hidden: " << stats.obs_hidden
                  << ",  Hidden-Hidden: " << stats.hidden_hidden << "\n";

        std::cout << "\n  Degrees:\n";
        std::cout << "    Leaves: " << stats.leaves
                  << ",  Max degree: " << stats.max_degree << "\n";
    }

    return stats;
}

}  // namespace clnj
