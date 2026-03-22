#pragma once

#include "types.h"
#include <string>
#include <vector>

namespace clnj {

struct FastaRecord {
    std::string name;
    std::string sequence;
};

std::vector<FastaRecord> parse_fasta(const std::string& path);

CleanResult load_clean_fasta(
    const std::string& fasta_path,
    int min_non_gap = 100,
    double max_ambiguity = 0.5,
    double length_percentile_cutoff = 5.0,
    bool skip_length_check = false
);

}  // namespace clnj
