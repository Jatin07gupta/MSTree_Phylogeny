#include "alignment_utils.h"
#include "fasta_parser.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace clnj {

bool is_aligned(const std::string& fasta_path) {
    auto records = parse_fasta(fasta_path);
    if (records.empty()) return true;
    size_t first_len = records[0].sequence.size();
    for (size_t i = 1; i < records.size(); ++i)
        if (records[i].sequence.size() != first_len) return false;
    return true;
}

std::string align_with_mafft(
    const std::string& input_fasta,
    const std::string& output_fasta_arg,
    int n_seqs,
    bool rna_struct
) {
    std::string output_fasta = output_fasta_arg;
    if (output_fasta.empty()) {
        auto dot = input_fasta.rfind('.');
        if (dot != std::string::npos)
            output_fasta = input_fasta.substr(0, dot) + "_aligned" + input_fasta.substr(dot);
        else
            output_fasta = input_fasta + "_aligned.fasta";
    }

    std::string cmd = "mafft";
    if (rna_struct)
        cmd += " --kimura 1 --xinsi";
    else if (n_seqs < 500)
        cmd += " --auto";
    else if (n_seqs < 10000)
        cmd += " --retree 2";
    else
        cmd += " --retree 1";

    cmd += " --preservecase " + input_fasta + " > " + output_fasta + " 2>/dev/null";

    std::cout << "  Running MAFFT: " << cmd << "\n";
    auto t0 = std::chrono::steady_clock::now();
    int ret = std::system(cmd.c_str());
    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    if (ret != 0)
        std::cerr << "  WARNING: MAFFT returned exit code " << ret << "\n";

    std::cout << "  Alignment written to " << output_fasta
              << "  (" << elapsed << "s)\n";
    return output_fasta;
}

}  // namespace clnj
