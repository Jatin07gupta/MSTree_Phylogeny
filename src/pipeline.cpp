#include "pipeline.h"

#include "alignment_utils.h"
#include "algorithm2.h"
#include "clnj.h"
#include "distance.h"
#include "fasta_parser.h"
#include "iqtree_interface.h"
#include "kmer_distance.h"
#include "mdl_clustering.h"
#include "mlvmst.h"
#include "online_insertion.h"
#include "one_hot.h"
#include "tree_analysis.h"
#include "tree_state.h"
#include "types.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace clnj {

static std::string sep(int n = 70) { return std::string(n, '\xe2' == '-' ? '-' : '-'); }

int run_pipeline(const PipelineArgs& args) {
    auto fasta_path = args.fasta_path;
    int min_non_gap = args.min_non_gap;
    int subsample_n = args.subsample_n;
    auto tree_algo = parse_tree_algo(args.tree_algo);
    std::string dist_model_str = args.model;
    for (auto& c : dist_model_str) c = (char)std::toupper(c);
    double gamma_alpha = args.gamma_alpha;
    auto distance_method = args.distance_method;
    bool is_iqtree = (dist_model_str == "IQTREE");

    auto algo_label = tree_algo_to_string(tree_algo);
    std::string gamma_label = (gamma_alpha > 0)
        ? " + Gamma(alpha=" + std::to_string(gamma_alpha) + ")" : "";
    std::string model_label = dist_model_str + gamma_label;
    if (distance_method == "kmer")
        model_label = "K-mer(k=" + std::to_string(args.kmer_size)
                    + ", sketch=" + std::to_string(args.sketch_size) + ")";

    std::mt19937 rng(SEED);

    std::cout << std::string(80, '=') << "\n";
    std::cout << "INTEGRATED PIPELINE: Multi-Model Distance -> MLVMST -> CLNJ\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << "  Input:             " << fasta_path << "\n";
    std::cout << "  Distance method:   " << distance_method << "\n";
    std::cout << "  min_non_gap:       " << min_non_gap << "\n";
    std::cout << "  min_shared_sites:  " << args.min_shared_sites << "\n";
    std::cout << "  max_ambiguity:     " << args.max_ambiguity << "\n";
    std::cout << "  subsample:         " << (subsample_n > 0 ? std::to_string(subsample_n) : "ALL") << "\n";
    std::cout << "  min_branch_length: " << MIN_BRANCH << "\n";
    std::cout << "  Distance model:    " << model_label << "\n";
    std::cout << "  Tree algorithm:    " << algo_label << "\n";
    std::cout << "  MST backend:       " << args.mst_algorithm << "\n";

    auto T0 = std::chrono::steady_clock::now();

    // ── STAGE 0: Alignment Detection ─────────────────────────────
    std::string aligned_path = fasta_path;

    if (distance_method == "kmer") {
        std::cout << "\n" << std::string(70, '-') << "\n";
        std::cout << "STAGE 0: Alignment Detection  [SKIPPED -- k-mer mode]\n";
        std::cout << std::string(70, '-') << "\n";
    } else {
        std::cout << "\n" << std::string(70, '-') << "\n";
        std::cout << "STAGE 0: Alignment Detection\n";
        std::cout << std::string(70, '-') << "\n";
        bool already_aligned = is_aligned(fasta_path);
        std::cout << "  Input aligned: " << (already_aligned ? "True" : "False") << "\n";

        if (!already_aligned && !args.no_align) {
            auto records = parse_fasta(fasta_path);
            int n_raw = (int)records.size();
            std::cout << "  Detected " << n_raw << " unaligned sequences -- calling MAFFT\n";
            aligned_path = align_with_mafft(fasta_path, "", n_raw, args.rna_struct);
        } else if (!already_aligned && args.no_align) {
            std::cout << "  WARNING: Input is unaligned but --no-align was set. "
                      << "Distance computation may fail.\n";
        }
    }

    // ── STAGE 1: Load and Clean ──────────────────────────────────
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "STAGE 1: Load, Clean, Name\n";
    std::cout << std::string(70, '-') << "\n";
    auto clean = load_clean_fasta(
        aligned_path, min_non_gap,
        args.max_ambiguity, args.length_percentile_cutoff,
        (distance_method == "kmer")
    );
    int n_total = (int)clean.names.size();

    if (subsample_n > 0 && n_total > subsample_n) {
        std::vector<int> indices(n_total);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), rng);
        indices.resize(subsample_n);
        std::sort(indices.begin(), indices.end());

        std::vector<std::string> new_names, new_seqs;
        for (int i : indices) {
            new_names.push_back(clean.names[i]);
            new_seqs.push_back(clean.sequences[i]);
        }
        clean.names = std::move(new_names);
        clean.sequences = std::move(new_seqs);
        std::cout << "  Subsampled: " << n_total << " -> " << subsample_n
                  << " (seed=" << SEED << ")\n";
    }

    int n = (int)clean.names.size();
    if (n == 0) {
        std::cout << "  ERROR: No sequences remain after filtering. Exiting.\n";
        return 1;
    }
    int m_len = (int)clean.sequences[0].size();
    if (distance_method == "kmer") {
        int min_len = m_len, max_len = m_len;
        long long sum_len = 0;
        for (auto& s : clean.sequences) {
            int l = (int)s.size();
            min_len = std::min(min_len, l);
            max_len = std::max(max_len, l);
            sum_len += l;
        }
        std::cout << "  Final: n=" << n << ", seq lengths: min=" << min_len
                  << ", max=" << max_len << ", mean=" << sum_len / n << "\n";
    } else {
        std::cout << "  Final: n=" << n << ", alignment length m=" << m_len << "\n";
    }

    MatrixXd D;
    double t_enc = 0.0, t_dist = 0.0;
    double saved_gA = 0.25, saved_gC = 0.25, saved_gG = 0.25, saved_gT = 0.25;

    // ── K-MER DISTANCE PATH ──────────────────────────────────────
    if (distance_method == "kmer") {
        std::cout << "\n" << std::string(70, '-') << "\n";
        std::cout << "STAGE 3: K-mer Distance Matrix (" << n << "x" << n
                  << ")  [" << model_label << "]\n";
        std::cout << std::string(70, '-') << "\n";
        auto t0 = std::chrono::steady_clock::now();
        D = kmer_distance_matrix(clean.sequences, args.kmer_size, args.sketch_size);
        auto t1 = std::chrono::steady_clock::now();
        t_dist = std::chrono::duration<double>(t1 - t0).count();
    } else {
        // ── STAGE 2: One-Hot Encoding ────────────────────────────
        std::cout << "\n" << std::string(70, '-') << "\n";
        std::cout << "STAGE 2: One-Hot Encoding\n";
        std::cout << std::string(70, '-') << "\n";
        auto t0 = std::chrono::steady_clock::now();
        auto oh = one_hot_encode(clean.sequences);
        auto t1 = std::chrono::steady_clock::now();
        t_enc = std::chrono::duration<double>(t1 - t0).count();
        std::cout << "  Shape:  (" << oh.n << ", " << oh.m << ", 4)  dtype: float32\n";
        std::cout << "  Memory: " << std::fixed << std::setprecision(1)
                  << oh.nbytes() / 1e6 << " MB\n";
        std::cout << "  Time:   " << std::setprecision(2) << t_enc << "s\n";

        // Valid sites per sequence
        double min_valid = 1e18, max_valid = 0, sum_valid = 0;
        std::vector<double> valids;
        for (int i = 0; i < oh.n; ++i) {
            double v = 0;
            for (int j = 0; j < oh.m; ++j)
                if (oh.is_valid(i, j)) ++v;
            min_valid = std::min(min_valid, v);
            max_valid = std::max(max_valid, v);
            sum_valid += v;
            valids.push_back(v);
        }
        std::sort(valids.begin(), valids.end());
        double median_valid = (valids.size() % 2 == 0)
            ? (valids[valids.size()/2-1] + valids[valids.size()/2]) / 2.0
            : valids[valids.size()/2];
        std::cout << "  Valid sites/seq: min=" << (int)min_valid
                  << ", max=" << (int)max_valid
                  << ", mean=" << (int)(sum_valid / oh.n)
                  << ", median=" << (int)median_valid << "\n";

        // ── STAGE 2b: Alignment Diagnostics ──────────────────────
        std::cout << "\n" << std::string(70, '-') << "\n";
        std::cout << "STAGE 2b: Alignment Diagnostics\n";
        std::cout << std::string(70, '-') << "\n";
        auto pair_stats = compute_pair_counts(oh);
        saved_gA = pair_stats.gA; saved_gC = pair_stats.gC;
        saved_gG = pair_stats.gG; saved_gT = pair_stats.gT;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  Base frequencies: A=" << pair_stats.gA << "  C=" << pair_stats.gC
                  << "  G=" << pair_stats.gG << "  T=" << pair_stats.gT << "\n";
        std::cout << std::setprecision(3);
        std::cout << "  Ts/Tv ratio:     " << pair_stats.ts_tv_ratio << "\n";
        std::cout << std::setprecision(4);
        std::cout << "  Mean p-distance: " << pair_stats.mean_p << "\n";
        std::cout << "  Freq uniform:    " << (pair_stats.freq_uniform ? "True" : "False") << "\n";

        // ── IQ-TREE Model Selection ──────────────────────────────
        DistModel dist_model = DistModel::JC69;
        if (is_iqtree) {
            std::cout << "\n" << std::string(70, '-') << "\n";
            std::cout << "STAGE 2c: IQ-TREE Model Selection (external)\n";
            std::cout << std::string(70, '-') << "\n";
            auto iq_result = run_iqtree_model_selection(aligned_path);
            dist_model = parse_model(iq_result.pipeline_model);
            if (iq_result.has_gamma() && gamma_alpha < 0)
                gamma_alpha = iq_result.gamma_alpha;
            gamma_label = (gamma_alpha > 0)
                ? " + Gamma(alpha=" + std::to_string(gamma_alpha) + ")" : "";
            model_label = model_to_string(dist_model) + gamma_label
                        + " [via IQ-TREE: " + iq_result.iqtree_model + "]";
            std::cout << "  Mapped to pipeline: " << model_label << "\n";
        } else {
            dist_model = parse_model(dist_model_str);
        }

        // ── STAGE 3: Distance Matrix ─────────────────────────────
        std::cout << "\n" << std::string(70, '-') << "\n";
        std::cout << "STAGE 3: Distance Matrix (" << n << "x" << n
                  << ")  [" << model_label << "]\n";
        std::cout << std::string(70, '-') << "\n";
        auto td0 = std::chrono::steady_clock::now();
        D = compute_distance_matrix(oh, dist_model, gamma_alpha, args.min_shared_sites, true);
        auto td1 = std::chrono::steady_clock::now();
        t_dist = std::chrono::duration<double>(td1 - td0).count();
    }

    // Distance matrix statistics
    std::vector<double> upper_vals;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            upper_vals.push_back(D(i, j));

    int total_pairs = (int)upper_vals.size();
    int n_zero = 0;
    std::set<double> unique_rounded;
    for (double v : upper_vals) {
        if (std::abs(v) < EPS) ++n_zero;
        unique_rounded.insert(std::round(v * 1e10) / 1e10);
    }
    int n_unique = (int)unique_rounded.size();

    std::sort(upper_vals.begin(), upper_vals.end());
    double upper_min = upper_vals.front();
    double upper_max = upper_vals.back();
    double upper_mean = std::accumulate(upper_vals.begin(), upper_vals.end(), 0.0) / total_pairs;
    double upper_median = upper_vals[total_pairs / 2];

    bool symmetric = true, diag_zero = true, no_nan = true;
    for (int i = 0; i < n && symmetric; ++i)
        for (int j = 0; j < n && symmetric; ++j)
            if (std::abs(D(i,j) - D(j,i)) > 1e-12) symmetric = false;
    for (int i = 0; i < n; ++i) if (D(i,i) != 0.0) diag_zero = false;
    for (double v : upper_vals) if (std::isnan(v)) no_nan = false;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Time:       " << t_dist << "s\n";
    std::cout << "  Symmetric:  " << (symmetric ? "True" : "False") << "\n";
    std::cout << "  Diag zero:  " << (diag_zero ? "True" : "False") << "\n";
    std::cout << "  No NaN:     " << (no_nan ? "True" : "False") << "\n";
    std::cout << std::setprecision(6);
    std::cout << "  Range:      [" << upper_min << ", " << upper_max << "]\n";
    std::cout << "  Mean:       " << upper_mean << "\n";
    std::cout << "  Median:     " << upper_median << "\n";
    std::cout << "  Zero dists: " << n_zero << "/" << total_pairs
              << " (" << std::setprecision(2) << 100.0 * n_zero / total_pairs << "%)\n";
    std::cout << "  Unique:     " << n_unique << "/" << total_pairs
              << " (" << std::setprecision(1) << 100.0 * n_unique / total_pairs << "%)\n";

    // Distance distribution histogram
    double bins[] = {0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 5.0, 10.0};
    int n_bins = 10;
    std::vector<int> hist(n_bins, 0);
    for (double v : upper_vals) {
        for (int b = n_bins - 1; b >= 0; --b) {
            if (v >= bins[b]) { ++hist[b]; break; }
        }
    }
    int max_hist = *std::max_element(hist.begin(), hist.end());
    std::cout << "\n  Distance distribution:\n";
    for (int b = 0; b < n_bins; ++b) {
        double pct = 100.0 * hist[b] / total_pairs;
        int bar_len = std::max(1, (int)(40.0 * hist[b] / std::max(max_hist, 1)));
        std::string bar(bar_len, '#');
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "    [" << bins[b] << ", " << bins[b+1] << ")  ";
        std::cout << std::setw(8) << hist[b] << "  ("
                  << std::setprecision(1) << std::setw(5) << pct << "%)  " << bar << "\n";
    }

    // ── STAGE 4: Algorithm 2 ─────────────────────────────────────
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "STAGE 4: Algorithm 2 (F_C, G_U)\n";
    std::cout << std::string(70, '-') << "\n";
    auto algo2 = construct_FC_and_GU(D);

    // ── STAGE 5: MLVMST ─────────────────────────────────────────
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "STAGE 5: MLVMST (Algorithm 3)\n";
    std::cout << std::string(70, '-') << "\n";
    auto mlvmst = build_mlvmst(D, algo2, "increasing", args.mst_algorithm);

    // ── STAGE 6: CLNJ ───────────────────────────────────────────
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "STAGE 6: CLNJ Tree Reconstruction  [" << algo_label << "]\n";
    std::cout << std::string(70, '-') << "\n";

    std::cout << "\n  [A] min_branch_length = " << MIN_BRANCH << "\n";
    auto tc0 = std::chrono::steady_clock::now();
    auto clnj_result = clnj_clean(D, mlvmst.adjacency, tree_algo, MIN_BRANCH, true,
                                  args.report_negatives, args.trace_zero_edges);
    auto tc1 = std::chrono::steady_clock::now();
    double t_clnj = std::chrono::duration<double>(tc1 - tc0).count();
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  CLNJ completed in " << t_clnj << "s\n";

    // Print diagnostic results if enabled
    if (args.report_negatives && !clnj_result.negative_counts.empty()) {
        std::cout << "\n  Negative-value diagnostic counters:\n";
        std::vector<std::pair<std::string, int>> sorted_neg(
            clnj_result.negative_counts.begin(), clnj_result.negative_counts.end());
        std::sort(sorted_neg.begin(), sorted_neg.end());
        for (auto& [key, cnt] : sorted_neg)
            std::cout << "    " << key << ": " << cnt << "\n";
    } else if (args.report_negatives) {
        std::cout << "\n  Negative-value diagnostic counters: (none detected)\n";
    }

    if (args.trace_zero_edges) {
        std::cout << "\n  Zero-edge log (" << clnj_result.zero_edge_log.size() << " entries):\n";
        if (!clnj_result.zero_edge_log.empty()) {
            std::map<std::string, int> case_counts;
            for (auto& entry : clnj_result.zero_edge_log)
                case_counts[entry.case_type]++;
            for (auto& [c, cnt] : case_counts)
                std::cout << "    " << c << ": " << cnt << "\n";
        }
    }

    std::cout << "\n  [B] min_branch_length = 0 (raw, for comparison)\n";
    auto clnj0 = clnj_clean(D, mlvmst.adjacency, tree_algo, 0.0, false,
                             args.report_negatives, args.trace_zero_edges);

    // ── STAGE 7: Tree Analysis ───────────────────────────────────
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "STAGE 7: Tree Analysis\n";
    std::cout << std::string(70, '-') << "\n";

    std::cout << "\n  WITH min_branch_length = 1e-6:";
    auto stats = analyze_tree(clnj_result.adjacency, clnj_result.edge_weights,
                              clnj_result.hidden_info, n, clean.names);

    std::cout << "\n  WITHOUT min_branch_length (raw):";
    auto stats0 = analyze_tree(clnj0.adjacency, clnj0.edge_weights,
                               clnj0.hidden_info, n, clean.names);

    // Edge weight comparison
    std::cout << "\n" << std::string(70, '-') << "\n";
    std::cout << "Edge Weight Comparison\n";
    std::cout << std::string(70, '-') << "\n";

    auto print_edge_stats = [](const std::string& label, const EdgeWeights& ew) {
        std::vector<double> w;
        for (auto& [ek, v] : ew) w.push_back(v);
        int neg = 0, zero = 0, pos = 0;
        for (double v : w) {
            if (v < -EPS) ++neg;
            else if (std::abs(v) < EPS) ++zero;
            else ++pos;
        }
        std::sort(w.begin(), w.end());
        double mean = std::accumulate(w.begin(), w.end(), 0.0) / w.size();
        std::cout << "\n  " << label << ":\n";
        std::cout << "    Negative: " << neg << "  Zero: " << zero
                  << "  Positive: " << pos << "  Total: " << w.size() << "\n";
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "    Min: " << w.front();
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  Max: " << w.back() << "  Mean: " << mean << "\n";
    };

    print_edge_stats("With floor (1e-6)", clnj_result.edge_weights);
    print_edge_stats("Raw (no floor)", clnj0.edge_weights);

    if (args.report_negatives && !clnj0.negative_counts.empty()) {
        std::cout << "\n  Negative counters (raw, no floor):\n";
        std::vector<std::pair<std::string, int>> sorted_neg0(
            clnj0.negative_counts.begin(), clnj0.negative_counts.end());
        std::sort(sorted_neg0.begin(), sorted_neg0.end());
        for (auto& [key, cnt] : sorted_neg0)
            std::cout << "    " << key << ": " << cnt << "\n";
    }

    if (args.trace_zero_edges) {
        std::cout << "\n  Zero-edge log (raw, no floor): "
                  << clnj0.zero_edge_log.size() << " entries\n";
        if (!clnj0.zero_edge_log.empty()) {
            std::map<std::string, int> cc0;
            for (auto& entry : clnj0.zero_edge_log)
                cc0[entry.case_type]++;
            for (auto& [c, cnt] : cc0)
                std::cout << "    " << c << ": " << cnt << "\n";
        }
    }

    // ── STAGE 8: Pipeline Summary ────────────────────────────────
    auto T1 = std::chrono::steady_clock::now();
    double total_time = std::chrono::duration<double>(T1 - T0).count();

    std::cout << std::fixed;
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "PIPELINE SUMMARY\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "  Dataset:              " << fasta_path.substr(fasta_path.rfind('/') + 1) << "\n";
    std::cout << "  Sequences used:       " << n << " (from " << n_total << " usable)\n";
    std::cout << "  Alignment length:     " << m_len << "\n";
    std::cout << "  Distance model:       " << model_label << "\n";
    std::cout << "  Pairwise deletion:    Yes\n";
    std::cout << "  Tree algorithm:       " << algo_label << "\n";
    std::cout << "  MST backend:          " << args.mst_algorithm << "\n";
    std::cout << std::setprecision(2);
    std::cout << "  Encoding time:        " << t_enc << "s\n";
    std::cout << "  Distance matrix time: " << t_dist << "s\n";
    std::cout << std::setprecision(1);
    std::cout << "  CLNJ time:            " << t_clnj << "s\n";
    std::cout << "  Total wall time:      " << total_time << "s ("
              << total_time / 60.0 << " min)\n";
    std::cout << "\n  --- Distance Matrix ---\n";
    std::cout << "  Zero distances:       " << n_zero << "/" << total_pairs << "\n";
    std::cout << std::setprecision(6);
    std::cout << "  Distance range:       [" << upper_min << ", " << upper_max << "]\n";
    std::cout << "\n  --- Tree (with floor) ---\n";
    std::cout << "  Hidden nodes:         " << stats.n_hidden << "\n";
    std::cout << "  Total edges:          " << stats.total_edges << "\n";
    std::cout << "  Tree valid:           " << (stats.valid_tree ? "True" : "False") << "\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  Min edge weight:      " << stats.min_weight << "\n";
    std::cout << "  Zero-length edges:    " << stats.zero_edges << "\n";
    std::cout << "\n  --- Tree (raw, no floor) ---\n";
    std::cout << "  Zero-length edges:    " << stats0.zero_edges << "\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  Min edge weight:      " << stats0.min_weight << "\n";
    std::cout << std::string(70, '=') << "\n";

    if (!args.save_state_path.empty()) {
        std::cout << "\n" << std::string(70, '-') << "\n";
        std::cout << "STAGE 9: MDL Clustering & Save State\n";
        std::cout << std::string(70, '-') << "\n";

        auto tc_mdl0 = std::chrono::steady_clock::now();
        auto ct = mdl_cluster_tree(
            clnj_result.adjacency, clnj_result.edge_weights,
            clnj_result.hidden_info, n,
            args.mdl_mj, args.mdl_merge_threshold, true
        );
        auto tc_mdl1 = std::chrono::steady_clock::now();
        std::cout << "  MDL clustering time: "
                  << std::chrono::duration<double>(tc_mdl1 - tc_mdl0).count() << "s\n";

        DistModel resolved_model = parse_model(dist_model_str);

        TreeState ts;
        ts.names = clean.names;
        ts.sequences = clean.sequences;
        ts.aln_len = m_len;
        ts.D = D;
        ts.adjacency = clnj_result.adjacency;
        ts.edge_weights = clnj_result.edge_weights;
        ts.hidden_info = clnj_result.hidden_info;
        ts.cluster_tree = std::move(ct);
        ts.model = resolved_model;
        ts.gamma_alpha = gamma_alpha;
        ts.gA = saved_gA; ts.gC = saved_gC;
        ts.gG = saved_gG; ts.gT = saved_gT;
        ts.min_shared_sites = args.min_shared_sites;
        ts.n_observed = n;
        ts.next_hidden_id = 0;
        for (const auto& [k, _] : ts.adjacency)
            ts.next_hidden_id = std::max(ts.next_hidden_id, k + 1);
        ts.distance_method = distance_method;
        ts.kmer_size = args.kmer_size;
        ts.sketch_size = args.sketch_size;
        ts.tree_algo = tree_algo;

        if (distance_method == "alignment") {
            int cn = n, cm = m_len;
            ts.cached_nuc_idx.resize((size_t)cn * cm);
            ts.cached_valid.resize((size_t)cn * cm);
            for (int i = 0; i < cn; ++i) {
                const auto& seq = ts.sequences[i];
                int slen = std::min((int)seq.size(), cm);
                for (int j = 0; j < slen; ++j) {
                    int8_t nuc = -1;
                    char c = seq[j];
                    if (c == 'A' || c == 'a') nuc = 0;
                    else if (c == 'C' || c == 'c') nuc = 1;
                    else if (c == 'G' || c == 'g') nuc = 2;
                    else if (c == 'T' || c == 't' || c == 'U' || c == 'u') nuc = 3;
                    ts.cached_nuc_idx[(size_t)i * cm + j] = nuc;
                    ts.cached_valid[(size_t)i * cm + j] = (nuc >= 0) ? 1 : 0;
                }
                for (int j = slen; j < cm; ++j) {
                    ts.cached_nuc_idx[(size_t)i * cm + j] = -1;
                    ts.cached_valid[(size_t)i * cm + j] = 0;
                }
            }
            std::cout << "  Cached nuc_idx + valid arrays: "
                      << cn << " x " << cm << " ("
                      << (cn * cm * 2) / 1024 << " KB)\n";
        }

        if (save_tree_state(ts, args.save_state_path))
            std::cout << "  State saved to: " << args.save_state_path << "\n";
        else
            std::cerr << "  ERROR: Failed to save state to " << args.save_state_path << "\n";
    }

    return 0;
}

int run_insertion_pipeline(const InsertionArgs& args) {
    std::cout << std::string(80, '=') << "\n";
    std::cout << "ONLINE INSERTION PIPELINE (CPU)\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << "  State file:  " << args.state_path << "\n";
    std::cout << "  New FASTA:   " << args.fasta_path << "\n";

    TreeState state;
    if (!load_tree_state(state, args.state_path)) {
        std::cerr << "ERROR: Failed to load tree state from " << args.state_path << "\n";
        return 1;
    }

    std::cout << "  Loaded state: " << state.n_observed << " observed nodes, "
              << state.hidden_info.size() << " hidden, "
              << state.cluster_tree.clusters.size() << " clusters\n";

    int n_leaf_clusters = 0;
    for (const auto& c : state.cluster_tree.clusters)
        if (c.child_cluster_ids.empty()) ++n_leaf_clusters;
    std::cout << "  Leaf clusters: " << n_leaf_clusters << "\n";

    auto records = parse_fasta(args.fasta_path);
    if (records.empty()) {
        std::cerr << "ERROR: No sequences found in " << args.fasta_path << "\n";
        return 1;
    }

    std::vector<std::string> new_names, new_seqs;
    for (const auto& rec : records) {
        new_names.push_back(rec.name);
        new_seqs.push_back(rec.sequence);
    }
    std::cout << "  New sequences: " << new_names.size() << "\n\n";

    auto t0 = std::chrono::steady_clock::now();
    auto results = insert_batch(state, new_names, new_seqs, args.verbose);
    auto t1 = std::chrono::steady_clock::now();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    int success = 0;
    for (const auto& r : results) if (r.success) ++success;

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "INSERTION SUMMARY\n";
    std::cout << std::string(70, '=') << "\n";
    std::cout << "  Total insertions attempted: " << results.size() << "\n";
    std::cout << "  Successful: " << success << "\n";
    std::cout << "  Final tree: " << state.adjacency.size() << " nodes, "
              << state.edge_weights.size() << " edges\n";
    std::cout << "  Final observed: " << state.n_observed << "\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Time: " << elapsed << "s\n";
    std::cout << std::string(70, '=') << "\n";

    std::string out_path = args.output_state_path.empty()
        ? args.state_path : args.output_state_path;
    if (save_tree_state(state, out_path))
        std::cout << "\n  Updated state saved to: " << out_path << "\n";
    else
        std::cerr << "\n  ERROR: Failed to save updated state to " << out_path << "\n";

    return 0;
}

int dump_state(const std::string& path) {
    TreeState state;
    if (!load_tree_state(state, path)) return 1;
    int n_hidden = (int)state.hidden_info.size();
    int n_edges = (int)state.edge_weights.size();
    int n_leaf = 0;
    for (const auto& [node, nbs] : state.adjacency)
        if (nbs.size() == 1) ++n_leaf;
    int leaf_clusters = 0;
    for (const auto& c : state.cluster_tree.clusters)
        if (c.child_cluster_ids.empty()) ++leaf_clusters;
    std::cout << "State: " << path << "\n";
    std::cout << "  Observed: " << state.n_observed << "  Hidden: " << n_hidden
              << "  Total nodes: " << (int)state.adjacency.size()
              << "  Edges: " << n_edges << "  Leaves: " << n_leaf << "\n";
    std::cout << "  Clusters: " << state.cluster_tree.clusters.size()
              << "  Leaf clusters: " << leaf_clusters << "\n";
    std::cout << "  Cached encoding: " << (state.has_cached_encoding() ? "yes" : "no") << "\n";
    return 0;
}

}  // namespace clnj
