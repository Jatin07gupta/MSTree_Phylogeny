#pragma once

#include "types.h"
#include <string>
#include <vector>

namespace clnj {

InsertionResult insert_taxon(
    TreeState& state,
    const std::string& new_name,
    const std::string& new_sequence,
    bool verbose = true
);

std::vector<InsertionResult> insert_batch(
    TreeState& state,
    const std::vector<std::string>& new_names,
    const std::vector<std::string>& new_sequences,
    bool verbose = true
);

}  // namespace clnj
