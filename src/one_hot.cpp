#include "one_hot.h"
#include <cstring>

namespace clnj {

// Lookup table: char -> one-hot [A, C, G, T]
static float ONEHOT_TABLE[256][4];
static bool table_initialized = false;

static void init_table() {
    if (table_initialized) return;
    std::memset(ONEHOT_TABLE, 0, sizeof(ONEHOT_TABLE));

    ONEHOT_TABLE['A'][0] = 1.0f; ONEHOT_TABLE['a'][0] = 1.0f;
    ONEHOT_TABLE['C'][1] = 1.0f; ONEHOT_TABLE['c'][1] = 1.0f;
    ONEHOT_TABLE['G'][2] = 1.0f; ONEHOT_TABLE['g'][2] = 1.0f;
    ONEHOT_TABLE['T'][3] = 1.0f; ONEHOT_TABLE['t'][3] = 1.0f;
    ONEHOT_TABLE['U'][3] = 1.0f; ONEHOT_TABLE['u'][3] = 1.0f;

    table_initialized = true;
}

OneHotData one_hot_encode(const std::vector<std::string>& sequences) {
    init_table();

    OneHotData result;
    if (sequences.empty()) return result;

    result.n = (int)sequences.size();
    result.m = (int)sequences[0].size();
    result.data.resize((size_t)result.n * result.m * 4, 0.0f);

    for (int i = 0; i < result.n; ++i) {
        const auto& seq = sequences[i];
        int len = std::min((int)seq.size(), result.m);
        for (int j = 0; j < len; ++j) {
            unsigned char c = static_cast<unsigned char>(seq[j]);
            int off = i * result.m * 4 + j * 4;
            result.data[off + 0] = ONEHOT_TABLE[c][0];
            result.data[off + 1] = ONEHOT_TABLE[c][1];
            result.data[off + 2] = ONEHOT_TABLE[c][2];
            result.data[off + 3] = ONEHOT_TABLE[c][3];
        }
    }
    return result;
}

}  // namespace clnj
