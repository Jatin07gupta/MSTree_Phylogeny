#pragma once

#include <string>

namespace clnj {

struct PipelineArgs {
    std::string fasta_path;
    int min_non_gap = 100;
    int subsample_n = 0;
    std::string tree_algo = "mfnj";
    std::string model = "JC69";
    double gamma_alpha = -1.0;
    bool no_align = false;
    bool rna_struct = false;
    int min_shared_sites = 100;
    double max_ambiguity = 0.5;
    double length_percentile_cutoff = 5.0;
    std::string distance_method = "alignment";
    std::string mst_algorithm = "kruskal";
    int kmer_size = 16;
    int sketch_size = 1000;
    bool report_negatives = false;
    bool trace_zero_edges = false;

    std::string save_state_path;
    double mdl_mj = 7.0;
    int mdl_merge_threshold = 100;
};

struct InsertionArgs {
    std::string state_path;
    std::string fasta_path;
    std::string output_state_path;
    bool verbose = true;
};

int run_pipeline(const PipelineArgs& args);
int run_insertion_pipeline(const InsertionArgs& args);

/** Load state from path and print key stats (for validation/debugging). */
int dump_state(const std::string& path);

}  // namespace clnj
