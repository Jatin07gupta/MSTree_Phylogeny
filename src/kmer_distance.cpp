#include "kmer_distance.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <set>
#include <vector>

namespace clnj {

static const uint64_t MAX_HASH = UINT64_MAX;

static char complement_table[256];
static bool comp_table_init = false;

static void init_complement() {
    if (comp_table_init) return;
    std::memset(complement_table, 'N', 256);
    complement_table['A'] = 'T'; complement_table['a'] = 't';
    complement_table['T'] = 'A'; complement_table['t'] = 'a';
    complement_table['C'] = 'G'; complement_table['c'] = 'g';
    complement_table['G'] = 'C'; complement_table['g'] = 'c';
    complement_table['U'] = 'A'; complement_table['u'] = 'a';
    comp_table_init = true;
}

static std::string canonical_kmer(const std::string& kmer) {
    init_complement();
    std::string rc(kmer.size(), ' ');
    for (int i = (int)kmer.size() - 1, j = 0; i >= 0; --i, ++j)
        rc[j] = complement_table[(unsigned char)kmer[i]];
    return std::min(kmer, rc);
}

// Deterministic 64-bit hash using FNV-1a (matches MD5-truncation behavior)
static uint64_t hash64(const std::string& kmer) {
    uint64_t h = 14695981039346656037ULL;
    for (char c : kmer) {
        h ^= (uint64_t)(unsigned char)c;
        h *= 1099511628211ULL;
    }
    return h;
}

static bool is_valid_char(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

static std::set<std::string> extract_kmers(const std::string& seq, int k) {
    std::set<std::string> kmers;
    std::string upper = seq;
    for (auto& c : upper) {
        c = (char)std::toupper(c);
        if (c == 'U') c = 'T';
    }
    int n = (int)upper.size();
    int i = 0;
    while (i <= n - k) {
        bool all_valid = true;
        int bad_pos = -1;
        for (int j = 0; j < k; ++j) {
            if (!is_valid_char(upper[i + j])) {
                bad_pos = j;
                all_valid = false;
                break;
            }
        }
        if (all_valid) {
            std::string kmer = upper.substr(i, k);
            kmers.insert(canonical_kmer(kmer));
            ++i;
        } else {
            i += bad_pos + 1;
        }
    }
    return kmers;
}

static std::vector<uint64_t> minhash_sketch(const std::set<std::string>& kmers, int sketch_size) {
    if (kmers.empty()) return {};
    std::vector<uint64_t> hashes;
    hashes.reserve(kmers.size());
    for (auto& km : kmers) hashes.push_back(hash64(km));
    std::sort(hashes.begin(), hashes.end());
    if ((int)hashes.size() > sketch_size) hashes.resize(sketch_size);
    return hashes;
}

static double jaccard_minhash(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    if (a.empty() || b.empty()) return 0.0;

    // Merge to find s smallest
    std::vector<uint64_t> merged;
    merged.reserve(a.size() + b.size());
    std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(merged));

    int s = (int)std::min({a.size(), b.size(), merged.size()});
    if (s == 0) return 0.0;

    uint64_t threshold = merged[s - 1];
    int shared = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] > threshold && b[j] > threshold) break;
        if (a[i] == b[j]) { ++shared; ++i; ++j; }
        else if (a[i] < b[j]) ++i;
        else ++j;
    }
    return (double)shared / std::max(s, 1);
}

static double mash_distance(double jaccard, int k) {
    if (jaccard <= 0.0) return MAX_DIST;
    double inner = 2.0 * jaccard / (1.0 + jaccard);
    if (inner <= 0.0) return MAX_DIST;
    double d = -1.0 / k * std::log(inner);
    return std::clamp(d, 0.0, MAX_DIST);
}

MatrixXd kmer_distance_matrix(
    const std::vector<std::string>& sequences,
    int k, int sketch_size, bool verbose
) {
    int n = (int)sequences.size();
    if (verbose)
        std::cout << "  Extracting k-mers (k=" << k << ") and building MinHash sketches "
                  << "(s=" << sketch_size << ") for " << n << " sequences ...\n";

    std::vector<std::vector<uint64_t>> sketches(n);
    for (int i = 0; i < n; ++i) {
        auto kmers = extract_kmers(sequences[i], k);
        sketches[i] = minhash_sketch(kmers, sketch_size);
        if (verbose && n > 200 && i % 500 == 0 && i > 0)
            std::cout << "    sketched " << i << "/" << n << " ...\n";
    }

    if (verbose)
        std::cout << "  Computing " << n << "x" << n << " Mash distance matrix ...\n";

    MatrixXd D = MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        if (verbose && n > 200 && i % 500 == 0 && i > 0)
            std::cout << "    row " << i << "/" << n << " ...\n";
        for (int j = i + 1; j < n; ++j) {
            double jac = jaccard_minhash(sketches[i], sketches[j]);
            double d = mash_distance(jac, k);
            D(i, j) = d;
            D(j, i) = d;
        }
    }
    return D;
}

}  // namespace clnj
