#include "distance.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <vector>

namespace clnj {

// ═══════════════════════════════════════════════════════════════════
//  GAMMA TRANSFORM
// ═══════════════════════════════════════════════════════════════════

static inline double gamma_transform(double x, double alpha) {
    x = std::max(x, LOG_FLOOR);
    if (alpha < 0.0)  // alpha not set => standard log
        return -std::log(x);
    return alpha * (std::pow(x, -1.0 / alpha) - 1.0);
}

// ═══════════════════════════════════════════════════════════════════
//  MODEL CORRECTION FUNCTIONS (single-value versions)
// ═══════════════════════════════════════════════════════════════════

double jc69_correct(double p, double alpha) {
    p = std::clamp(p, 0.0, 0.74999);
    double inner = 1.0 - (4.0 / 3.0) * p;
    double d = 0.75 * gamma_transform(inner, alpha);
    if (!std::isfinite(d)) d = 0.0;
    return std::clamp(d, 0.0, MAX_DIST);
}

double k2p_correct(double P, double Q, double alpha) {
    double w1 = std::max(1.0 - 2.0 * P - Q, FREQ_FLOOR);
    double w2 = std::max(1.0 - 2.0 * Q, FREQ_FLOOR);
    double d = 0.5 * gamma_transform(w1, alpha) + 0.25 * gamma_transform(w2, alpha);
    if (!std::isfinite(d)) d = 0.0;
    return std::clamp(d, 0.0, MAX_DIST);
}

double tn93_correct(double P1, double P2, double Qv,
                           double gA, double gC, double gG, double gT,
                           double alpha) {
    double gR = std::max(gA + gG, FREQ_FLOOR);
    double gY = std::max(gC + gT, FREQ_FLOOR);
    double _gA = std::max(gA, FREQ_FLOOR);
    double _gC = std::max(gC, FREQ_FLOOR);
    double _gG = std::max(gG, FREQ_FLOOR);
    double _gT = std::max(gT, FREQ_FLOOR);

    double K1 = 2.0 * _gA * _gG / gR;
    double K2 = 2.0 * _gT * _gC / gY;
    double K3 = 2.0 * (gR * gY - _gA * _gG * gY / gR - _gT * _gC * gR / gY);

    double w1 = std::max(1.0 - P1 / K1 - Qv / (2.0 * gR), FREQ_FLOOR);
    double w2 = std::max(1.0 - P2 / K2 - Qv / (2.0 * gY), FREQ_FLOOR);
    double w3 = std::max(1.0 - Qv / (2.0 * gR * gY), FREQ_FLOOR);

    double d = K1 * gamma_transform(w1, alpha)
             + K2 * gamma_transform(w2, alpha)
             + K3 * gamma_transform(w3, alpha);
    if (!std::isfinite(d)) d = 0.0;
    return std::clamp(d, 0.0, MAX_DIST);
}

// ═══════════════════════════════════════════════════════════════════
//  ALIGNMENT DIAGNOSTICS
// ═══════════════════════════════════════════════════════════════════

PairStats compute_pair_counts(const OneHotData& oh, int n_sample_pairs, int seed) {
    PairStats stats;
    int n = oh.n, m = oh.m;

    // Global base frequencies
    double total_valid = 0.0;
    double count_A = 0, count_C = 0, count_G = 0, count_T = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (oh.is_valid(i, j)) {
                total_valid += 1.0;
                int nuc = oh.argmax(i, j);
                switch (nuc) {
                    case 0: count_A++; break;
                    case 1: count_C++; break;
                    case 2: count_G++; break;
                    case 3: count_T++; break;
                }
            }
        }
    }
    if (total_valid > 0) {
        stats.gA = count_A / total_valid;
        stats.gC = count_C / total_valid;
        stats.gG = count_G / total_valid;
        stats.gT = count_T / total_valid;
    }

    // Sample pairs
    std::mt19937 rng(seed);
    int max_pairs = n * (n - 1) / 2;
    int n_sample = std::min(n_sample_pairs, max_pairs);

    std::set<std::pair<int,int>> pairs;
    std::uniform_int_distribution<int> dist(0, n - 1);
    int attempts = 0;
    while ((int)pairs.size() < n_sample && attempts < n_sample * 10) {
        int a = dist(rng), b = dist(rng);
        ++attempts;
        if (a != b) pairs.insert({std::min(a,b), std::max(a,b)});
    }

    std::vector<double> p_vals, P1_vals, P2_vals, Q_vals;
    for (auto& [pi, pj] : pairs) {
        int nv = 0;
        for (int s = 0; s < m; ++s)
            if (oh.is_valid(pi, s) && oh.is_valid(pj, s)) ++nv;
        if (nv < 10) continue;

        int nd = 0, n_ts_AG = 0, n_ts_CT = 0;
        for (int s = 0; s < m; ++s) {
            if (!oh.is_valid(pi, s) || !oh.is_valid(pj, s)) continue;
            int ni = oh.argmax(pi, s);
            int nj_v = oh.argmax(pj, s);
            if (ni != nj_v) {
                ++nd;
                int sum = ni + nj_v;
                if (sum == 2) ++n_ts_AG;      // A(0) <-> G(2)
                else if (sum == 4) ++n_ts_CT;  // C(1) <-> T(3)
            }
        }
        double nv_f = (double)nv;
        int n_tv = nd - n_ts_AG - n_ts_CT;
        p_vals.push_back(nd / nv_f);
        P1_vals.push_back(n_ts_AG / nv_f);
        P2_vals.push_back(n_ts_CT / nv_f);
        Q_vals.push_back(n_tv / nv_f);
    }

    if (p_vals.empty()) return stats;

    double sum_ts = 0, sum_tv = 0;
    for (size_t i = 0; i < p_vals.size(); ++i) {
        sum_ts += P1_vals[i] + P2_vals[i];
        sum_tv += Q_vals[i];
    }
    double mean_ts = sum_ts / p_vals.size();
    double mean_tv = sum_tv / p_vals.size();
    stats.ts_tv_ratio = mean_ts / std::max(mean_tv, 1e-10);

    double sum_p = 0, sum_P1 = 0, sum_P2 = 0, sum_Q = 0;
    for (size_t i = 0; i < p_vals.size(); ++i) {
        sum_p += p_vals[i];
        sum_P1 += P1_vals[i];
        sum_P2 += P2_vals[i];
        sum_Q += Q_vals[i];
    }
    stats.mean_p  = sum_p / p_vals.size();
    stats.mean_P1 = sum_P1 / p_vals.size();
    stats.mean_P2 = sum_P2 / p_vals.size();
    stats.mean_Q  = sum_Q / p_vals.size();
    stats.n_sampled = (int)p_vals.size();

    stats.freq_uniform = (std::abs(stats.gA - 0.25) < 0.05 &&
                          std::abs(stats.gC - 0.25) < 0.05 &&
                          std::abs(stats.gG - 0.25) < 0.05 &&
                          std::abs(stats.gT - 0.25) < 0.05);
    return stats;
}

// ═══════════════════════════════════════════════════════════════════
//  AUTO MODEL SELECTION DIAGNOSTICS
// ═══════════════════════════════════════════════════════════════════

static double estimate_invariant_proportion(const OneHotData& oh) {
    int n = oh.n, m = oh.m;
    int n_invariant = 0;
    for (int j = 0; j < m; ++j) {
        int first_nuc = -1;
        bool invariant = true;
        int count_valid = 0;
        for (int i = 0; i < n; ++i) {
            if (!oh.is_valid(i, j)) continue;
            ++count_valid;
            int nuc = oh.argmax(i, j);
            if (first_nuc < 0) first_nuc = nuc;
            else if (nuc != first_nuc) { invariant = false; break; }
        }
        if (count_valid >= 2 && invariant) ++n_invariant;
    }
    return (double)n_invariant / std::max(m, 1);
}

static double estimate_site_rate_cv(const OneHotData& oh) {
    int n = oh.n, m = oh.m;
    std::vector<double> site_rates(m, 0.0);
    for (int j = 0; j < m; ++j) {
        int counts[4] = {};
        int n_valid = 0;
        for (int i = 0; i < n; ++i) {
            if (!oh.is_valid(i, j)) continue;
            counts[oh.argmax(i, j)]++;
            ++n_valid;
        }
        if (n_valid < 2) continue;
        int mode_count = *std::max_element(counts, counts + 4);
        site_rates[j] = 1.0 - (double)mode_count / n_valid;
    }

    std::vector<double> variable;
    for (double r : site_rates)
        if (r > 0.0) variable.push_back(r);
    if (variable.size() < 10) return 0.0;

    double mean = std::accumulate(variable.begin(), variable.end(), 0.0) / variable.size();
    double sq_sum = 0.0;
    for (double v : variable) sq_sum += (v - mean) * (v - mean);
    double stddev = std::sqrt(sq_sum / variable.size());
    return stddev / std::max(mean, 1e-10);
}

static double estimate_gamma_alpha(double site_rate_cv) {
    if (site_rate_cv <= 0.01) return -1.0;
    double alpha = 1.0 / std::max(site_rate_cv * site_rate_cv, 0.01);
    alpha = std::clamp(alpha, 0.05, 100.0);
    return std::round(alpha * 100.0) / 100.0;
}

static double per_taxon_composition_test(const OneHotData& oh, int n_sample = 200, int seed = 42) {
    int n = oh.n, m = oh.m;

    double total_valid = 0;
    double global_counts[4] = {};
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (!oh.is_valid(i, j)) continue;
            total_valid += 1.0;
            global_counts[oh.argmax(i, j)] += 1.0;
        }
    }
    double global_freq[4];
    for (int b = 0; b < 4; ++b)
        global_freq[b] = global_counts[b] / std::max(total_valid, 1.0);

    std::mt19937 rng(seed);
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
    int to_sample = std::min(n_sample, n);

    std::vector<double> chi2_vals;
    for (int idx = 0; idx < to_sample; ++idx) {
        int i = indices[idx];
        int nv = 0;
        double obs[4] = {};
        for (int j = 0; j < m; ++j) {
            if (!oh.is_valid(i, j)) continue;
            ++nv;
            obs[oh.argmax(i, j)] += 1.0;
        }
        if (nv < 20) continue;

        double chi2 = 0.0;
        for (int b = 0; b < 4; ++b) {
            double expected = std::max(global_freq[b] * nv, 0.5);
            double diff = obs[b] - expected;
            chi2 += diff * diff / expected;
        }
        chi2_vals.push_back(chi2);
    }

    if (chi2_vals.empty()) return 1.0;
    double mean_chi2 = std::accumulate(chi2_vals.begin(), chi2_vals.end(), 0.0) / chi2_vals.size();
    double p_approx = (mean_chi2 < 50.0) ? std::exp(-mean_chi2 / 2.0) : 0.0;
    return p_approx;
}

ModelSelection select_model(const OneHotData& oh, bool verbose) {
    PairStats stats = compute_pair_counts(oh);
    double inv_prop = estimate_invariant_proportion(oh);
    double site_cv = estimate_site_rate_cv(oh);
    double comp_p = per_taxon_composition_test(oh);
    double gamma_est = estimate_gamma_alpha(site_cv);

    if (verbose) {
        std::cout << "  Model selection diagnostics (improved):\n"
                  << "    Base freqs:       A=" << stats.gA << "  C=" << stats.gC
                  << "  G=" << stats.gG << "  T=" << stats.gT << "\n"
                  << "    Ts/Tv ratio:      " << stats.ts_tv_ratio << "\n"
                  << "    Mean p-dist:      " << stats.mean_p << "\n"
                  << "    Freq uniform:     " << (stats.freq_uniform ? "True" : "False") << "\n"
                  << "    Invariant sites:  " << inv_prop << "\n"
                  << "    Site-rate CV:     " << site_cv << "\n"
                  << "    Composition p:    " << comp_p << "\n"
                  << "    Gamma alpha est:  " << (gamma_est > 0 ? std::to_string(gamma_est) : "None") << "\n";
    }

    bool recommend_gamma = site_cv > 1.0;
    ModelSelection result;
    result.gamma_alpha = recommend_gamma ? gamma_est : -1.0;

    if (comp_p < 0.01)
        result.model = "LOGDET";
    else if (!stats.freq_uniform)
        result.model = "TN93";
    else if (stats.ts_tv_ratio >= 1.2)
        result.model = "K2P";
    else
        result.model = "JC69";

    if (verbose)
        std::cout << "    -> Selected model: " << result.model << "\n";

    return result;
}

// ═══════════════════════════════════════════════════════════════════
//  FULL DISTANCE MATRIX
// ═══════════════════════════════════════════════════════════════════

MatrixXd compute_distance_matrix(
    const OneHotData& oh,
    DistModel model,
    double gamma_alpha,
    int min_shared_sites,
    bool verbose
) {
    int n = oh.n, m = oh.m;

    if (model == DistModel::AUTO) {
        auto sel = select_model(oh, verbose);
        model = parse_model(sel.model);
        if (sel.has_gamma() && gamma_alpha < 0.0) {
            gamma_alpha = sel.gamma_alpha;
            if (verbose)
                std::cout << "  AUTO -> applying gamma_alpha=" << gamma_alpha << "\n";
        }
    }

    // Compute global base frequencies for TN93
    double gA = 0.25, gC = 0.25, gG = 0.25, gT = 0.25;
    if (model == DistModel::TN93) {
        double tv_all = 0;
        double cA = 0, cC = 0, cG = 0, cT = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (!oh.is_valid(i, j)) continue;
                tv_all += 1.0;
                switch (oh.argmax(i, j)) {
                    case 0: cA++; break; case 1: cC++; break;
                    case 2: cG++; break; case 3: cT++; break;
                }
            }
        }
        if (tv_all > 0) {
            gA = cA / tv_all; gC = cC / tv_all;
            gG = cG / tv_all; gT = cT / tv_all;
        }
        if (verbose)
            std::cout << "  TN93 global freqs: A=" << gA << " C=" << gC
                      << " G=" << gG << " T=" << gT << "\n";
    }

    std::string gamma_label = (gamma_alpha > 0) ?
        "+G(a=" + std::to_string(gamma_alpha) + ")" : "";
    if (verbose)
        std::cout << "  Computing " << n << "x" << n << " distance matrix  "
                  << "[model=" << model_to_string(model) << gamma_label
                  << "]  min_shared_sites=" << min_shared_sites << "\n";

    double min_sites_threshold = std::max(min_shared_sites, 1);

    // Pre-compute nuc_idx and valid arrays for efficiency
    std::vector<std::vector<int>> nuc_idx(n, std::vector<int>(m));
    std::vector<std::vector<bool>> valid(n, std::vector<bool>(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            valid[i][j] = oh.is_valid(i, j);
            nuc_idx[i][j] = oh.argmax(i, j);
        }
    }

    MatrixXd D = MatrixXd::Zero(n, n);

    for (int i = 0; i < n; ++i) {
        if (verbose && n > 200 && i % 500 == 0 && i > 0)
            std::cout << "    row " << i << "/" << n << " ...\n";

        for (int j = i + 1; j < n; ++j) {
            int n_valid = 0;
            int n_diff = 0, n_ts_AG = 0, n_ts_CT = 0;

            for (int s = 0; s < m; ++s) {
                if (!valid[i][s] || !valid[j][s]) continue;
                ++n_valid;
                int ni = nuc_idx[i][s];
                int nj_v = nuc_idx[j][s];
                if (ni != nj_v) {
                    ++n_diff;
                    int sum = ni + nj_v;
                    if (sum == 2) ++n_ts_AG;
                    else if (sum == 4) ++n_ts_CT;
                }
            }

            double d = 0.0;

            if (n_valid < min_sites_threshold) {
                d = MAX_DIST;
            } else if (model == DistModel::LOGDET) {
                // Build 4x4 joint frequency matrix
                double F[4][4] = {};
                for (int s = 0; s < m; ++s) {
                    if (!valid[i][s] || !valid[j][s]) continue;
                    F[nuc_idx[i][s]][nuc_idx[j][s]] += 1.0;
                }
                double nv = std::max((double)n_valid, 1.0);
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                        F[a][b] /= nv;

                // Compute determinant of 4x4 matrix
                Eigen::Matrix4d Fmat;
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                        Fmat(a, b) = F[a][b];
                double det = Fmat.determinant();
                d = -std::log(std::max(std::abs(det), LOG_FLOOR)) / 4.0;
                if (!std::isfinite(d)) d = 0.0;
                d = std::max(d, 0.0);
            } else {
                double nv = (double)n_valid;
                int n_tv = n_diff - n_ts_AG - n_ts_CT;

                if (model == DistModel::JC69) {
                    double p = (double)n_diff / nv;
                    d = jc69_correct(p, gamma_alpha);
                } else if (model == DistModel::K2P) {
                    double P = (double)(n_ts_AG + n_ts_CT) / nv;
                    double Qv = (double)n_tv / nv;
                    d = k2p_correct(P, Qv, gamma_alpha);
                } else if (model == DistModel::TN93) {
                    double P1 = (double)n_ts_AG / nv;
                    double P2 = (double)n_ts_CT / nv;
                    double Qv = (double)n_tv / nv;
                    d = tn93_correct(P1, P2, Qv, gA, gC, gG, gT, gamma_alpha);
                }
            }

            d = std::max(d, 0.0);
            D(i, j) = d;
            D(j, i) = d;
        }
    }

    return D;
}

// ═══════════════════════════════════════════════════════════════════
//  Subset distance and encoding (for online insertion)
// ═══════════════════════════════════════════════════════════════════

static int char_to_nuc(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return -1;
    }
}

void encode_sequence_flat(
    const std::string& seq, int m,
    int8_t* out_nuc, uint8_t* out_valid
) {
    int slen = std::min((int)seq.size(), m);
    for (int s = 0; s < slen; ++s) {
        int8_t nuc = -1;
        char c = seq[s];
        if (c == 'A' || c == 'a') nuc = 0;
        else if (c == 'C' || c == 'c') nuc = 1;
        else if (c == 'G' || c == 'g') nuc = 2;
        else if (c == 'T' || c == 't' || c == 'U' || c == 'u') nuc = 3;
        out_nuc[s] = nuc;
        out_valid[s] = (nuc >= 0) ? 1 : 0;
    }
    for (int s = slen; s < m; ++s) {
        out_nuc[s] = -1;
        out_valid[s] = 0;
    }
}

std::unordered_map<int, double> compute_distance_to_subset(
    const OneHotData& oh_existing,
    const std::string& new_sequence,
    const std::vector<int>& target_indices,
    DistModel model,
    double gamma_alpha,
    double gA, double gC, double gG, double gT,
    int min_shared_sites
) {
    int m = oh_existing.m;
    int seq_len = (int)new_sequence.size();
    int sites = std::min(m, seq_len);

    std::vector<int> new_nuc(sites);
    std::vector<bool> new_valid(sites);
    for (int s = 0; s < sites; ++s) {
        int nuc = char_to_nuc(new_sequence[s]);
        new_nuc[s] = nuc;
        new_valid[s] = (nuc >= 0);
    }

    std::unordered_map<int, double> result;
    double min_sites_threshold = std::max(min_shared_sites, 1);

    for (int i : target_indices) {
        if (i < 0 || i >= oh_existing.n) continue;

        int n_valid = 0, n_diff = 0, n_ts_AG = 0, n_ts_CT = 0;

        for (int s = 0; s < sites; ++s) {
            if (!new_valid[s] || !oh_existing.is_valid(i, s)) continue;
            ++n_valid;
            int ni = oh_existing.argmax(i, s);
            int nj = new_nuc[s];
            if (ni != nj) {
                ++n_diff;
                int sum = ni + nj;
                if (sum == 2) ++n_ts_AG;
                else if (sum == 4) ++n_ts_CT;
            }
        }

        double d = 0.0;
        if (n_valid < min_sites_threshold) {
            d = MAX_DIST;
        } else if (model == DistModel::LOGDET) {
            double F[4][4] = {};
            for (int s = 0; s < sites; ++s) {
                if (!new_valid[s] || !oh_existing.is_valid(i, s)) continue;
                F[oh_existing.argmax(i, s)][new_nuc[s]] += 1.0;
            }
            double nv = std::max((double)n_valid, 1.0);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    F[a][b] /= nv;
            Eigen::Matrix4d Fmat;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    Fmat(a, b) = F[a][b];
            double det = Fmat.determinant();
            d = -std::log(std::max(std::abs(det), LOG_FLOOR)) / 4.0;
            if (!std::isfinite(d)) d = 0.0;
            d = std::max(d, 0.0);
        } else {
            double nv = (double)n_valid;
            int n_tv = n_diff - n_ts_AG - n_ts_CT;
            if (model == DistModel::JC69) {
                d = jc69_correct((double)n_diff / nv, gamma_alpha);
            } else if (model == DistModel::K2P) {
                d = k2p_correct((double)(n_ts_AG + n_ts_CT) / nv,
                                (double)n_tv / nv, gamma_alpha);
            } else if (model == DistModel::TN93) {
                d = tn93_correct((double)n_ts_AG / nv, (double)n_ts_CT / nv,
                                 (double)n_tv / nv, gA, gC, gG, gT, gamma_alpha);
            }
        }
        d = std::max(d, 0.0);
        result[i] = d;
    }
    return result;
}

std::unordered_map<int, double> compute_distance_to_subset_cached(
    const int8_t* nuc_idx,
    const uint8_t* valid,
    int n_existing, int m,
    const int8_t* new_nuc,
    const uint8_t* new_valid,
    const std::vector<int>& target_indices,
    DistModel model,
    double gamma_alpha,
    double gA, double gC, double gG, double gT,
    int min_shared_sites
) {
    std::unordered_map<int, double> result;
    double min_sites_threshold = std::max(min_shared_sites, 1);

    for (int i : target_indices) {
        if (i < 0 || i >= n_existing) continue;

        const int8_t* row_nuc = nuc_idx + (size_t)i * m;
        const uint8_t* row_val = valid + (size_t)i * m;

        int n_valid = 0, n_diff = 0, n_ts_AG = 0, n_ts_CT = 0;

        for (int s = 0; s < m; ++s) {
            if (!new_valid[s] || !row_val[s]) continue;
            ++n_valid;
            int ni = row_nuc[s];
            int nj = new_nuc[s];
            if (ni != nj) {
                ++n_diff;
                int sum = ni + nj;
                if (sum == 2) ++n_ts_AG;
                else if (sum == 4) ++n_ts_CT;
            }
        }

        double d = 0.0;
        if (n_valid < min_sites_threshold) {
            d = MAX_DIST;
        } else if (model == DistModel::LOGDET) {
            double F[4][4] = {};
            for (int s = 0; s < m; ++s) {
                if (!new_valid[s] || !row_val[s]) continue;
                F[row_nuc[s]][new_nuc[s]] += 1.0;
            }
            double nv = std::max((double)n_valid, 1.0);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    F[a][b] /= nv;
            Eigen::Matrix4d Fmat;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    Fmat(a, b) = F[a][b];
            double det = Fmat.determinant();
            d = -std::log(std::max(std::abs(det), LOG_FLOOR)) / 4.0;
            if (!std::isfinite(d)) d = 0.0;
            d = std::max(d, 0.0);
        } else {
            double nv = (double)n_valid;
            int n_tv = n_diff - n_ts_AG - n_ts_CT;
            if (model == DistModel::JC69) {
                d = jc69_correct((double)n_diff / nv, gamma_alpha);
            } else if (model == DistModel::K2P) {
                d = k2p_correct((double)(n_ts_AG + n_ts_CT) / nv,
                                (double)n_tv / nv, gamma_alpha);
            } else if (model == DistModel::TN93) {
                d = tn93_correct((double)n_ts_AG / nv, (double)n_ts_CT / nv,
                                 (double)n_tv / nv, gA, gC, gG, gT, gamma_alpha);
            }
        }
        d = std::max(d, 0.0);
        result[i] = d;
    }
    return result;
}

}  // namespace clnj
