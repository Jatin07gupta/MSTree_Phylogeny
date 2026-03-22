#include "fasta_parser.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_map>

namespace clnj {

static const std::set<char> VALID_NUC = {
    'A','C','G','T','U','a','c','g','t','u'
};
static const std::set<char> AMBIGUITY = {
    'R','Y','S','W','K','M','B','D','H','V','N',
    'r','y','s','w','k','m','b','d','h','v','n'
};

std::vector<FastaRecord> parse_fasta(const std::string& path) {
    std::vector<FastaRecord> records;
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "  ERROR: Cannot open " << path << "\n";
        return records;
    }

    std::string line;
    FastaRecord current;
    bool has_record = false;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (has_record) {
                records.push_back(std::move(current));
                current = FastaRecord{};
            }
            current.name = line.substr(1);
            // Trim trailing whitespace
            while (!current.name.empty() && std::isspace(current.name.back()))
                current.name.pop_back();
            has_record = true;
        } else {
            // Remove trailing whitespace/CR
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' '))
                line.pop_back();
            current.sequence += line;
        }
    }
    if (has_record) records.push_back(std::move(current));

    return records;
}

CleanResult load_clean_fasta(
    const std::string& fasta_path,
    int min_non_gap,
    double max_ambiguity,
    double length_percentile_cutoff,
    bool skip_length_check
) {
    auto records = parse_fasta(fasta_path);
    std::cout << "  Raw sequences loaded: " << records.size() << "\n";

    int aln_len = 0;
    for (auto& r : records)
        aln_len = std::max(aln_len, (int)r.sequence.size());

    if (skip_length_check)
        std::cout << "  Max length: " << aln_len << "\n";
    else
        std::cout << "  Alignment length: " << aln_len << "\n";

    struct KeptEntry {
        FastaRecord record;
        int non_gap;
    };

    std::vector<KeptEntry> kept;
    int rm_divider = 0, rm_short = 0, rm_low = 0, rm_ambig = 0, rm_outlier = 0;

    for (auto& r : records) {
        if (r.name.find("DIVIDER") != std::string::npos) { ++rm_divider; continue; }
        if (!skip_length_check && (int)r.sequence.size() < aln_len) { ++rm_short; continue; }

        int non_gap = 0, n_ambig = 0;
        for (char c : r.sequence) {
            if (VALID_NUC.count(c)) ++non_gap;
            if (AMBIGUITY.count(c)) ++n_ambig;
        }
        if (non_gap < min_non_gap) { ++rm_low; continue; }

        double frac_ambig = (double)n_ambig / std::max(non_gap + n_ambig, 1);
        if (frac_ambig > max_ambiguity) { ++rm_ambig; continue; }

        kept.push_back({std::move(r), non_gap});
    }

    // Length-percentile outlier removal
    if (!kept.empty() && length_percentile_cutoff > 0.0) {
        std::vector<int> valid_counts;
        valid_counts.reserve(kept.size());
        for (auto& k : kept) valid_counts.push_back(k.non_gap);
        std::sort(valid_counts.begin(), valid_counts.end());

        double idx = (length_percentile_cutoff / 100.0) * (valid_counts.size() - 1);
        int lo = (int)std::floor(idx);
        int hi = std::min(lo + 1, (int)valid_counts.size() - 1);
        double frac = idx - lo;
        double cutoff = valid_counts[lo] * (1.0 - frac) + valid_counts[hi] * frac;

        std::vector<KeptEntry> filtered;
        for (auto& k : kept) {
            if (k.non_gap >= cutoff)
                filtered.push_back(std::move(k));
            else
                ++rm_outlier;
        }
        kept = std::move(filtered);
        std::cout << "  Length-outlier cutoff (" << length_percentile_cutoff
                  << "th pctile): " << (int)cutoff << " valid sites\n";
    }

    std::cout << "  Removed: {'divider': " << rm_divider
              << ", 'short': " << rm_short
              << ", 'low_content': " << rm_low
              << ", 'high_ambiguity': " << rm_ambig
              << ", 'length_outlier': " << rm_outlier << "}\n";
    std::cout << "  Usable:  " << kept.size() << "\n";

    // Unique naming (matches Python: uses :: splitting or record id)
    std::unordered_map<std::string, int> name_counts;
    for (auto& k : kept) {
        auto& desc = k.record.name;
        std::string base;
        auto pos = desc.find("::");
        if (pos != std::string::npos) {
            auto pos2 = desc.find("::", pos + 2);
            if (pos2 != std::string::npos)
                base = desc.substr(pos + 2, pos2 - pos - 2);
            else
                base = desc.substr(pos + 2);
        }
        // Trim trailing dots
        while (!base.empty() && base.back() == '.') base.pop_back();
        if (base.empty() || base.find_first_not_of(' ') == std::string::npos) {
            // Use first word as ID
            std::istringstream iss(desc);
            iss >> base;
        }
        name_counts[base]++;
    }

    CleanResult result;
    result.aln_len = aln_len;
    std::unordered_map<std::string, int> seen;

    for (auto& k : kept) {
        auto& desc = k.record.name;
        std::string base;
        auto pos = desc.find("::");
        if (pos != std::string::npos) {
            auto pos2 = desc.find("::", pos + 2);
            if (pos2 != std::string::npos)
                base = desc.substr(pos + 2, pos2 - pos - 2);
            else
                base = desc.substr(pos + 2);
        }
        while (!base.empty() && base.back() == '.') base.pop_back();
        if (base.empty() || base.find_first_not_of(' ') == std::string::npos) {
            std::istringstream iss(desc);
            iss >> base;
        }

        if (name_counts[base] > 1) {
            seen[base]++;
            result.names.push_back(base + "__v" + std::to_string(seen[base]));
        } else {
            result.names.push_back(base);
        }

        // Uppercase the sequence
        std::string seq = k.record.sequence;
        for (auto& c : seq) c = static_cast<char>(std::toupper(c));
        result.sequences.push_back(std::move(seq));
    }

    return result;
}

}  // namespace clnj
