#pragma once

#include "types.h"
#include "one_hot.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace clnj {

double jc69_correct(double p, double alpha);
double k2p_correct(double P, double Q, double alpha);
double tn93_correct(double P1, double P2, double Qv,
                    double gA, double gC, double gG, double gT, double alpha);

PairStats compute_pair_counts(const OneHotData& oh,
                              int n_sample_pairs = 2000, int seed = 42);

ModelSelection select_model(const OneHotData& oh, bool verbose = true);

MatrixXd compute_distance_matrix(
    const OneHotData& oh,
    DistModel model = DistModel::JC69,
    double gamma_alpha = -1.0,   // negative = not set
    int min_shared_sites = 1,
    bool verbose = false
);

std::unordered_map<int, double> compute_distance_to_subset(
    const OneHotData& oh_existing,
    const std::string& new_sequence,
    const std::vector<int>& target_indices,
    DistModel model = DistModel::JC69,
    double gamma_alpha = -1.0,
    double gA = 0.25, double gC = 0.25,
    double gG = 0.25, double gT = 0.25,
    int min_shared_sites = 1
);

std::unordered_map<int, double> compute_distance_to_subset_cached(
    const int8_t* nuc_idx,
    const uint8_t* valid,
    int n_existing, int m,
    const int8_t* new_nuc,
    const uint8_t* new_valid,
    const std::vector<int>& target_indices,
    DistModel model = DistModel::JC69,
    double gamma_alpha = -1.0,
    double gA = 0.25, double gC = 0.25,
    double gG = 0.25, double gT = 0.25,
    int min_shared_sites = 1
);

void encode_sequence_flat(
    const std::string& seq, int m,
    int8_t* out_nuc, uint8_t* out_valid
);

}  // namespace clnj
