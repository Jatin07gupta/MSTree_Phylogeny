#pragma once

#include "types.h"

namespace clnj {

ClnjResult clnj_clean(
    const MatrixXd& D,
    const Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>& mst_adj,
    TreeAlgo algo = TreeAlgo::MFNJ,
    double min_branch_length = MIN_BRANCH,
    bool verbose = true,
    bool report_negatives = false,
    bool trace_zero_edges = false
);

}  // namespace clnj
