#pragma once

#include "types.h"
#include <string>
#include <vector>

namespace clnj {

MatrixXd kmer_distance_matrix(
    const std::vector<std::string>& sequences,
    int k = 16,
    int sketch_size = 1000,
    bool verbose = true
);

}  // namespace clnj
