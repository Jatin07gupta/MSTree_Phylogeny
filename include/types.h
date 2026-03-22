#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace clnj {

// ── Pipeline constants (match Python exactly) ────────────────────

constexpr double MAX_DIST       = 10.0;
constexpr double MIN_BRANCH     = 1e-6;
constexpr double EPS            = 1e-12;
constexpr double TIE_TOL        = 1e-9;
constexpr int    SEED           = 42;
constexpr double LOG_FLOOR      = 1e-300;
constexpr double FREQ_FLOOR     = 1e-10;

// ── Eigen type aliases ───────────────────────────────────────────

using MatrixXd = Eigen::MatrixXd;
using VectorXd = Eigen::VectorXd;
using VectorXi = Eigen::VectorXi;

// ── Tree data structures ─────────────────────────────────────────

using Adjacency  = std::unordered_map<int, std::unordered_set<int>>;
using EdgeWeights = std::map<std::pair<int,int>, double>;

inline std::pair<int,int> normalize_edge(int a, int b) {
    return {std::min(a,b), std::max(a,b)};
}

// ── Result structures ────────────────────────────────────────────

struct CleanResult {
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    int aln_len = 0;
};

struct PairStats {
    double gA = 0.25, gC = 0.25, gG = 0.25, gT = 0.25;
    double ts_tv_ratio = 1.0;
    double mean_p = 0.0;
    double mean_P1 = 0.0, mean_P2 = 0.0, mean_Q = 0.0;
    bool   freq_uniform = true;
    int    n_sampled = 0;
};

struct ModelSelection {
    std::string model;
    double gamma_alpha = -1.0;  // negative means "not set"
    bool has_gamma() const { return gamma_alpha > 0.0; }
};

struct IqtreeResult {
    std::string iqtree_model;
    std::string pipeline_model;
    double gamma_alpha = -1.0;
    double bic = std::numeric_limits<double>::infinity();
    bool has_gamma() const { return gamma_alpha > 0.0; }
};

struct Algo2Result {
    std::vector<std::set<int>> F_C;
    std::set<std::pair<int,int>> G_U_edges;
};

struct MlvmstResult {
    Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic> adjacency;
    std::vector<int> delta_max;
    int leaf_count = 0;
};

// ── Diagnostic structures ────────────────────────────────────────

using NegCounts = std::unordered_map<std::string, int>;

struct ZeroEdgeEntry {
    std::string case_type;
    int a = -1;
    int b = -1;
    double w = 0.0;
    double d_raw = 0.0;
    double d_ij = 0.0;
    double delta_i_raw = 0.0;
    double delta_j_raw = 0.0;
    std::string reason;
    std::vector<int> nodes;
};

struct LocalResult {
    std::vector<std::tuple<int, int, double>> edges;
    std::vector<int> new_hidden_ids;
};

struct ClnjResult {
    Adjacency adjacency;
    EdgeWeights edge_weights;
    std::unordered_map<int, int> hidden_info;  // hidden_id -> 0 (placeholder)
    NegCounts negative_counts;
    std::vector<ZeroEdgeEntry> zero_edge_log;
};

struct TreeStats {
    int n_obs = 0;
    int n_hidden = 0;
    int total_nodes = 0;
    int total_edges = 0;
    bool valid_tree = false;
    double min_weight = 0.0;
    double max_weight = 0.0;
    double mean_weight = 0.0;
    double median_weight = 0.0;
    double std_weight = 0.0;
    int zero_edges = 0;
    int negative_edges = 0;
    int positive_edges = 0;
    int obs_obs = 0;
    int obs_hidden = 0;
    int hidden_hidden = 0;
    int leaves = 0;
    int max_degree = 0;
};

// ── Distance model enum ──────────────────────────────────────────

enum class DistModel {
    JC69, K2P, TN93, LOGDET, AUTO
};

inline DistModel parse_model(const std::string& s) {
    std::string u = s;
    for (auto& c : u) c = static_cast<char>(std::toupper(c));
    if (u == "JC69")   return DistModel::JC69;
    if (u == "K2P")    return DistModel::K2P;
    if (u == "TN93")   return DistModel::TN93;
    if (u == "LOGDET") return DistModel::LOGDET;
    if (u == "AUTO")   return DistModel::AUTO;
    return DistModel::JC69;
}

inline std::string model_to_string(DistModel m) {
    switch (m) {
        case DistModel::JC69:   return "JC69";
        case DistModel::K2P:    return "K2P";
        case DistModel::TN93:   return "TN93";
        case DistModel::LOGDET: return "LOGDET";
        case DistModel::AUTO:   return "AUTO";
    }
    return "JC69";
}

enum class TreeAlgo { NJ, MFNJ, BIONJ };

inline TreeAlgo parse_tree_algo(const std::string& s) {
    if (s == "nj")    return TreeAlgo::NJ;
    if (s == "mfnj")  return TreeAlgo::MFNJ;
    if (s == "bionj") return TreeAlgo::BIONJ;
    return TreeAlgo::MFNJ;
}

inline std::string tree_algo_to_string(TreeAlgo a) {
    switch (a) {
        case TreeAlgo::NJ:    return "Standard NJ";
        case TreeAlgo::MFNJ:  return "MFNJ";
        case TreeAlgo::BIONJ: return "BioNJ";
    }
    return "MFNJ";
}

// ── MDL Clustering structures (paper Algorithm 3) ────────────────

struct Cluster {
    int id = -1;
    std::set<int> observed_members;
    int center_node = -1;
    double description_length = 0.0;
    std::vector<int> child_cluster_ids;
    int parent_cluster_id = -1;
};

struct ClusterTree {
    std::vector<Cluster> clusters;
    std::unordered_map<int, int> node_to_cluster;
    double mj = 7.0;

    int find_cluster(int obs_node) const {
        auto it = node_to_cluster.find(obs_node);
        return (it != node_to_cluster.end()) ? it->second : -1;
    }
};

// ── Dynamic insertion result ─────────────────────────────────────

struct InsertionResult {
    int new_obs_id = -1;
    int target_cluster_id = -1;
    int nodes_in_local_rebuild = 0;
    bool success = false;
    std::string message;
};

// ── Persistent tree state for online insertion ───────────────────

struct TreeState {
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    int aln_len = 0;

    MatrixXd D;

    Adjacency adjacency;
    EdgeWeights edge_weights;
    std::unordered_map<int, int> hidden_info;

    ClusterTree cluster_tree;

    DistModel model = DistModel::JC69;
    double gamma_alpha = -1.0;
    double gA = 0.25, gC = 0.25, gG = 0.25, gT = 0.25;
    int min_shared_sites = 100;

    int n_observed = 0;
    int next_hidden_id = 0;

    std::string distance_method = "alignment";
    int kmer_size = 16;
    int sketch_size = 1000;
    TreeAlgo tree_algo = TreeAlgo::MFNJ;

    // Cached flat arrays for fast distance computation (avoid re-encoding).
    std::vector<int8_t>  cached_nuc_idx;
    std::vector<uint8_t> cached_valid;

    bool has_cached_encoding() const {
        return !cached_nuc_idx.empty() &&
               (int)cached_nuc_idx.size() == n_observed * aln_len;
    }
};

}  // namespace clnj
