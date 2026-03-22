#pragma once

#include "types.h"
#include "distance_oracle.h"
#include <vector>

namespace clnj {

LocalResult bionj_local(
    const std::vector<int>& nodes,
    DistanceOracle& oracle,
    int m,
    int next_hidden_id,
    double min_branch_length = 0.0,
    NegCounts* neg_cnt = nullptr,
    std::vector<ZeroEdgeEntry>* zero_edge_log = nullptr
);

}  // namespace clnj
