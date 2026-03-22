#include "clnj.h"
#include "distance_oracle.h"
#include "nj.h"
#include "mfnj.h"
#include "bionj.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

namespace clnj {

using LocalFn = std::function<LocalResult(
    const std::vector<int>&, DistanceOracle&, int, int, double,
    NegCounts*, std::vector<ZeroEdgeEntry>*)>;

static void assert_tree_invariant(const Adjacency& adj, const EdgeWeights& ew) {
    int n = (int)adj.size();
    int num_edges = 0;
    for (auto& [v, neighbors] : adj)
        num_edges += (int)neighbors.size();
    num_edges /= 2;
    if (n > 0 && num_edges != n - 1)
        std::cerr << "  WARNING: Tree invariant violated: |E|=" << num_edges
                  << " != |V|-1=" << n - 1 << "\n";
}

ClnjResult clnj_clean(
    const MatrixXd& D,
    const Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>& mst_adj,
    TreeAlgo algo,
    double min_branch_length,
    bool verbose,
    bool report_negatives,
    bool trace_zero_edges
) {
    int m = (int)D.rows();

    ClnjResult result;
    for (int i = 0; i < m; ++i)
        result.adjacency[i] = {};

    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            if (mst_adj(i, j) > 0) {
                result.adjacency[i].insert(j);
                result.adjacency[j].insert(i);
                result.edge_weights[{i, j}] = D(i, j);
            }
        }
    }

    DistanceOracle oracle(D, m, result.adjacency, result.edge_weights);

    NegCounts* neg_ptr = report_negatives ? &result.negative_counts : nullptr;
    std::vector<ZeroEdgeEntry>* zel_ptr = trace_zero_edges ? &result.zero_edge_log : nullptr;

    LocalFn local_fn;
    if (algo == TreeAlgo::BIONJ)
        local_fn = bionj_local;
    else if (algo == TreeAlgo::MFNJ)
        local_fn = mfnj_local;
    else
        local_fn = nj_local;

    std::vector<int> degrees(m, 0);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
            if (mst_adj(i, j) > 0) ++degrees[i];

    std::vector<int> internal_nodes;
    for (int i = 0; i < m; ++i)
        if (degrees[i] > 1) internal_nodes.push_back(i);

    std::sort(internal_nodes.begin(), internal_nodes.end(),
              [&](int a, int b) {
                  return degrees[a] > degrees[b] || (degrees[a] == degrees[b] && a < b);
              });

    if (verbose) {
        std::cout << "CLNJ Clean: Processing " << m << " observed nodes\n";
        std::cout << "Found " << internal_nodes.size() << " internal nodes to process\n";
    }

    for (int iter = 0; iter < (int)internal_nodes.size(); ++iter) {
        int center = internal_nodes[iter];
        auto& neighbors = result.adjacency[center];

        std::set<int> S_v;
        S_v.insert(center);
        for (int nb : neighbors) S_v.insert(nb);
        std::vector<int> S_v_list(S_v.begin(), S_v.end());
        std::sort(S_v_list.begin(), S_v_list.end());

        int next_hidden_id = 0;
        for (auto& [k, _] : result.adjacency)
            next_hidden_id = std::max(next_hidden_id, k + 1);

        if (verbose) {
            std::cout << "  Iteration " << iter + 1 << "/" << internal_nodes.size()
                      << ": node " << center
                      << " (degree " << neighbors.size()
                      << ", size " << S_v.size() << ")\n";
        }

        auto local_result = local_fn(
            S_v_list, oracle, m, next_hidden_id, min_branch_length,
            neg_ptr, zel_ptr
        );

        if (verbose && !local_result.new_hidden_ids.empty()) {
            std::cout << "    Created " << local_result.new_hidden_ids.size()
                      << " hidden node(s)\n";
        }

        for (int h_id : local_result.new_hidden_ids) {
            result.hidden_info[h_id] = 0;
            result.adjacency[h_id] = {};
        }

        for (int n1 : S_v) {
            for (int n2 : S_v) {
                if (n1 != n2 && result.adjacency.count(n1) &&
                    result.adjacency[n1].count(n2)) {
                    result.adjacency[n1].erase(n2);
                    result.adjacency[n2].erase(n1);
                    auto ek = normalize_edge(n1, n2);
                    result.edge_weights.erase(ek);
                }
            }
        }

        for (auto& [n1, n2, dist] : local_result.edges) {
            double d = (dist >= 0) ? std::max(min_branch_length, dist) : min_branch_length;
            result.adjacency[n1].insert(n2);
            result.adjacency[n2].insert(n1);
            auto ek = normalize_edge(n1, n2);
            result.edge_weights[ek] = d;
        }

        assert_tree_invariant(result.adjacency, result.edge_weights);
    }

    if (verbose) {
        int total_nodes = (int)result.adjacency.size();
        int total_edges = (int)result.edge_weights.size();
        int num_hidden = (int)result.hidden_info.size();
        std::cout << "CLNJ complete: " << total_nodes
                  << " nodes (" << m << " observed + " << num_hidden << " hidden)\n";
        std::cout << "  Total edges: " << total_edges << "\n";
    }

    return result;
}

}  // namespace clnj
