#pragma once

#include <string>

namespace clnj {

bool is_aligned(const std::string& fasta_path);

std::string align_with_mafft(
    const std::string& input_fasta,
    const std::string& output_fasta = "",
    int n_seqs = 0,
    bool rna_struct = false
);

}  // namespace clnj
