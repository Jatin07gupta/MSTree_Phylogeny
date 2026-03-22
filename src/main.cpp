#include "pipeline.h"
#include <CLI/CLI.hpp>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    CLI::App app{"Multi-model phylogenetic pipeline: Distance -> MLVMST -> CLNJ (C++)"};
    app.require_subcommand(0, 1);

    auto* build_cmd = app.add_subcommand("build", "Build phylogenetic tree from FASTA");
    clnj::PipelineArgs build_args;

    build_cmd->add_option("fasta", build_args.fasta_path,
                          "Path to FASTA file (aligned or unaligned)")
             ->required();
    build_cmd->add_option("min_non_gap", build_args.min_non_gap,
                          "Minimum valid nucleotides per sequence (default: 100)");
    build_cmd->add_option("subsample_n", build_args.subsample_n,
                          "Subsample to N sequences; 0 = use all (default: 0)");
    build_cmd->add_option("tree_algo", build_args.tree_algo,
                          "Tree reconstruction algorithm: nj, mfnj, bionj (default: mfnj)")
             ->check(CLI::IsMember({"nj", "mfnj", "bionj"}));

    build_cmd->add_option("--model,-m", build_args.model,
                          "Distance model: JC69, K2P, TN93, LOGDET, AUTO, IQTREE (default: JC69)")
             ->check(CLI::IsMember({"JC69","jc69","K2P","k2p","TN93","tn93",
                                    "LOGDET","logdet","AUTO","auto","IQTREE","iqtree"}));
    build_cmd->add_option("--gamma-alpha,-g", build_args.gamma_alpha,
                          "Gamma shape parameter for rate heterogeneity (default: none)");
    build_cmd->add_flag("--no-align", build_args.no_align,
                        "Skip automatic MAFFT alignment even if input is unaligned");
    build_cmd->add_flag("--rna-struct", build_args.rna_struct,
                        "Use MAFFT X-INS-i for RNA structure-aware alignment");
    build_cmd->add_option("--min-shared-sites", build_args.min_shared_sites,
                          "Minimum shared valid sites for a valid pairwise distance (default: 100)");
    build_cmd->add_option("--max-ambiguity", build_args.max_ambiguity,
                          "Max fraction of ambiguity codes (default: 0.5)");
    build_cmd->add_option("--length-percentile-cutoff", build_args.length_percentile_cutoff,
                          "Remove sequences below this percentile (default: 5.0)");
    build_cmd->add_option("--distance-method", build_args.distance_method,
                          "Distance computation method: alignment, kmer (default: alignment)")
             ->check(CLI::IsMember({"alignment", "kmer"}));
    build_cmd->add_option("--mst-algorithm", build_args.mst_algorithm,
                          "MST backend for MLVMST stage: kruskal, boruvka (default: kruskal)")
             ->check(CLI::IsMember({"kruskal", "boruvka"}));
    build_cmd->add_option("--kmer-size", build_args.kmer_size,
                          "K-mer size for alignment-free distance (default: 16)");
    build_cmd->add_option("--sketch-size", build_args.sketch_size,
                          "MinHash sketch size (default: 1000)");
    build_cmd->add_flag("--report-negatives", build_args.report_negatives,
                        "Track and report negative-value diagnostic counters");
    build_cmd->add_flag("--trace-zero-edges", build_args.trace_zero_edges,
                        "Track and report zero-length edges created during CLNJ");
    build_cmd->add_option("--save-state", build_args.save_state_path,
                          "Save tree state (including MDL clusters) for online insertion");
    build_cmd->add_option("--mdl-mj", build_args.mdl_mj,
                          "MDL free parameter mj (default: 7.0)");
    build_cmd->add_option("--mdl-merge-threshold", build_args.mdl_merge_threshold,
                          "Merge clusters smaller than this (default: 100)");

    auto* insert_cmd = app.add_subcommand("insert", "Insert new taxa into existing tree");
    clnj::InsertionArgs insert_args;

    insert_cmd->add_option("state", insert_args.state_path,
                           "Path to saved tree state file")
              ->required();
    insert_cmd->add_option("fasta", insert_args.fasta_path,
                           "Path to FASTA file with new sequences to insert")
              ->required();
    insert_cmd->add_option("--output-state,-o", insert_args.output_state_path,
                           "Output path for updated state (default: overwrite input)");
    insert_cmd->add_flag("--verbose,!--quiet", insert_args.verbose,
                         "Verbose output (default: true)");

    auto* dump_cmd = app.add_subcommand("dump", "Dump tree state stats (for validation)");
    std::string dump_path;
    dump_cmd->add_option("state", dump_path, "Path to tree state file")->required();

    CLI11_PARSE(app, argc, argv);

    if (dump_cmd->parsed()) {
        return clnj::dump_state(dump_path);
    }
    if (build_cmd->parsed()) {
        return clnj::run_pipeline(build_args);
    }
    if (insert_cmd->parsed()) {
        return clnj::run_insertion_pipeline(insert_args);
    }

    if (app.get_subcommands().empty()) {
        std::cerr << "Usage: clnj_pipeline <build|insert|dump> [options]\n"
                  << "Run 'clnj_pipeline build --help', 'clnj_pipeline insert --help', or 'clnj_pipeline dump --help' for details.\n";
        return 1;
    }

    return 0;
}
