#pragma once

#include <string>
#include <vector>

namespace clnj {

// One-hot encoding stored as flat array: shape (n, m, 4), row-major.
// Access: data[i * m * 4 + j * 4 + base]
struct OneHotData {
    std::vector<float> data;
    int n = 0;  // number of sequences
    int m = 0;  // alignment length

    float& operator()(int seq, int site, int base) {
        return data[seq * m * 4 + site * 4 + base];
    }
    const float& operator()(int seq, int site, int base) const {
        return data[seq * m * 4 + site * 4 + base];
    }
    bool is_valid(int seq, int site) const {
        int off = seq * m * 4 + site * 4;
        return data[off] + data[off+1] + data[off+2] + data[off+3] > 0.0f;
    }
    int argmax(int seq, int site) const {
        int off = seq * m * 4 + site * 4;
        int best = 0;
        float best_val = data[off];
        for (int b = 1; b < 4; ++b) {
            if (data[off + b] > best_val) {
                best_val = data[off + b];
                best = b;
            }
        }
        return best;
    }
    size_t nbytes() const { return data.size() * sizeof(float); }
};

OneHotData one_hot_encode(const std::vector<std::string>& sequences);

}  // namespace clnj
