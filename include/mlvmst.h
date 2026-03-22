#pragma once

#include "types.h"
#include <string>

namespace clnj {

MlvmstResult build_mlvmst(
    const MatrixXd& D,
    const Algo2Result& algo2,
    const std::string& ordering = "increasing",
    const std::string& mst_algorithm = "kruskal",
    bool verbose = true
);

}  // namespace clnj
