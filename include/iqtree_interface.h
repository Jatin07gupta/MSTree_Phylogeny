#pragma once

#include "types.h"
#include <string>

namespace clnj {

IqtreeResult run_iqtree_model_selection(
    const std::string& fasta_path,
    bool verbose = true
);

}  // namespace clnj
