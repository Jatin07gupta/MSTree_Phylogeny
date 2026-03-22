#pragma once

#include "types.h"
#include <string>

namespace clnj {

bool save_tree_state(const TreeState& state, const std::string& path);
bool load_tree_state(TreeState& state, const std::string& path);

}  // namespace clnj
