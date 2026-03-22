#pragma once

#include "types.h"

namespace clnj {

double get_tree_distance(int a, int b,
                         const Adjacency& adjacency,
                         const EdgeWeights& edge_weights);

class DistanceOracle {
public:
    DistanceOracle(const MatrixXd& D_obs, int m,
                   Adjacency& adjacency, EdgeWeights& edge_weights);

    double dist(int i, int j) const;
    double dist_bfs(int i, int j) const;
    double dist_for_neighbourhood(int i, int j, bool has_hidden) const;

private:
    const MatrixXd& D_obs_;
    int m_;
    Adjacency& adjacency_;
    EdgeWeights& edge_weights_;
};

}  // namespace clnj
