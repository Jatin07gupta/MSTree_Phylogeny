#include "tree_state.h"

#include <cstdint>
#include <fstream>
#include <iostream>

namespace clnj {

static const uint32_t MAGIC = 0x434C4E4A;  // "CLNJ"
static const uint32_t VERSION = 2;

template<typename T>
static void write_val(std::ofstream& f, const T& v) { f.write(reinterpret_cast<const char*>(&v), sizeof(T)); }
template<typename T>
static void read_val(std::ifstream& f, T& v) { f.read(reinterpret_cast<char*>(&v), sizeof(T)); }

static void write_string(std::ofstream& f, const std::string& s) {
    uint32_t len = (uint32_t)s.size();
    write_val(f, len);
    f.write(s.data(), len);
}
static void read_string(std::ifstream& f, std::string& s) {
    uint32_t len;
    read_val(f, len);
    s.resize(len);
    f.read(s.data(), len);
}

bool save_tree_state(const TreeState& state, const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    if (!f) { std::cerr << "Cannot open " << path << " for writing\n"; return false; }

    write_val(f, MAGIC);
    write_val(f, VERSION);

    write_val(f, (int32_t)state.n_observed);
    write_val(f, (int32_t)state.next_hidden_id);
    write_val(f, (int32_t)state.aln_len);
    write_val(f, (int32_t)state.model);
    write_val(f, state.gamma_alpha);
    write_val(f, state.gA); write_val(f, state.gC);
    write_val(f, state.gG); write_val(f, state.gT);
    write_val(f, (int32_t)state.min_shared_sites);
    write_val(f, (int32_t)state.tree_algo);
    write_val(f, (int32_t)state.kmer_size);
    write_val(f, (int32_t)state.sketch_size);
    write_string(f, state.distance_method);

    uint32_t n_names = (uint32_t)state.names.size();
    write_val(f, n_names);
    for (const auto& name : state.names) write_string(f, name);
    for (const auto& seq : state.sequences) write_string(f, seq);

    int rows = (int)state.D.rows(), cols = (int)state.D.cols();
    write_val(f, (int32_t)rows);
    write_val(f, (int32_t)cols);
    f.write(reinterpret_cast<const char*>(state.D.data()), rows * cols * sizeof(double));

    uint32_t n_adj_nodes = (uint32_t)state.adjacency.size();
    write_val(f, n_adj_nodes);
    for (const auto& [node, neighbors] : state.adjacency) {
        write_val(f, (int32_t)node);
        uint32_t n_nb = (uint32_t)neighbors.size();
        write_val(f, n_nb);
        for (int nb : neighbors) write_val(f, (int32_t)nb);
    }

    uint32_t n_edges = (uint32_t)state.edge_weights.size();
    write_val(f, n_edges);
    for (const auto& [ek, w] : state.edge_weights) {
        write_val(f, (int32_t)ek.first);
        write_val(f, (int32_t)ek.second);
        write_val(f, w);
    }

    uint32_t n_hidden = (uint32_t)state.hidden_info.size();
    write_val(f, n_hidden);
    for (const auto& [hid, val] : state.hidden_info) {
        write_val(f, (int32_t)hid);
        write_val(f, (int32_t)val);
    }

    write_val(f, state.cluster_tree.mj);
    uint32_t n_clusters = (uint32_t)state.cluster_tree.clusters.size();
    write_val(f, n_clusters);
    for (const auto& c : state.cluster_tree.clusters) {
        write_val(f, (int32_t)c.id);
        write_val(f, (int32_t)c.center_node);
        write_val(f, c.description_length);
        write_val(f, (int32_t)c.parent_cluster_id);
        uint32_t nm = (uint32_t)c.observed_members.size();
        write_val(f, nm);
        for (int m : c.observed_members) write_val(f, (int32_t)m);
        uint32_t nc = (uint32_t)c.child_cluster_ids.size();
        write_val(f, nc);
        for (int ch : c.child_cluster_ids) write_val(f, (int32_t)ch);
    }

    uint32_t n_map = (uint32_t)state.cluster_tree.node_to_cluster.size();
    write_val(f, n_map);
    for (const auto& [node, cid] : state.cluster_tree.node_to_cluster) {
        write_val(f, (int32_t)node);
        write_val(f, (int32_t)cid);
    }

    uint32_t cache_n = state.has_cached_encoding() ? (uint32_t)state.n_observed : 0;
    uint32_t cache_m = state.has_cached_encoding() ? (uint32_t)state.aln_len : 0;
    write_val(f, cache_n);
    write_val(f, cache_m);
    if (cache_n > 0 && cache_m > 0) {
        f.write(reinterpret_cast<const char*>(state.cached_nuc_idx.data()),
                cache_n * cache_m * sizeof(int8_t));
        f.write(reinterpret_cast<const char*>(state.cached_valid.data()),
                cache_n * cache_m * sizeof(uint8_t));
    }

    return f.good();
}

bool load_tree_state(TreeState& state, const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { std::cerr << "Cannot open " << path << " for reading\n"; return false; }

    uint32_t magic, version;
    read_val(f, magic);
    read_val(f, version);
    if (magic != MAGIC || (version != 1 && version != VERSION)) {
        std::cerr << "Invalid tree state file or version mismatch (got v"
                  << version << ", support v1-v" << VERSION << ")\n";
        return false;
    }

    int32_t i32;
    read_val(f, i32); state.n_observed = i32;
    read_val(f, i32); state.next_hidden_id = i32;
    read_val(f, i32); state.aln_len = i32;
    read_val(f, i32); state.model = static_cast<DistModel>(i32);
    read_val(f, state.gamma_alpha);
    read_val(f, state.gA); read_val(f, state.gC);
    read_val(f, state.gG); read_val(f, state.gT);
    read_val(f, i32); state.min_shared_sites = i32;
    read_val(f, i32); state.tree_algo = static_cast<TreeAlgo>(i32);
    read_val(f, i32); state.kmer_size = i32;
    read_val(f, i32); state.sketch_size = i32;
    read_string(f, state.distance_method);

    uint32_t n_names;
    read_val(f, n_names);
    state.names.resize(n_names);
    state.sequences.resize(n_names);
    for (uint32_t i = 0; i < n_names; ++i) read_string(f, state.names[i]);
    for (uint32_t i = 0; i < n_names; ++i) read_string(f, state.sequences[i]);

    int32_t rows, cols;
    read_val(f, rows); read_val(f, cols);
    state.D.resize(rows, cols);
    f.read(reinterpret_cast<char*>(state.D.data()), rows * cols * sizeof(double));

    uint32_t n_adj_nodes;
    read_val(f, n_adj_nodes);
    state.adjacency.clear();
    for (uint32_t i = 0; i < n_adj_nodes; ++i) {
        int32_t node; read_val(f, node);
        uint32_t n_nb; read_val(f, n_nb);
        auto& nbs = state.adjacency[node];
        for (uint32_t j = 0; j < n_nb; ++j) {
            int32_t nb; read_val(f, nb);
            nbs.insert(nb);
        }
    }

    uint32_t n_edges;
    read_val(f, n_edges);
    state.edge_weights.clear();
    for (uint32_t i = 0; i < n_edges; ++i) {
        int32_t a, b; double w;
        read_val(f, a); read_val(f, b); read_val(f, w);
        state.edge_weights[{a, b}] = w;
    }

    uint32_t n_hidden;
    read_val(f, n_hidden);
    state.hidden_info.clear();
    for (uint32_t i = 0; i < n_hidden; ++i) {
        int32_t hid, val;
        read_val(f, hid); read_val(f, val);
        state.hidden_info[hid] = val;
    }

    read_val(f, state.cluster_tree.mj);
    uint32_t n_clusters;
    read_val(f, n_clusters);
    state.cluster_tree.clusters.resize(n_clusters);
    for (uint32_t i = 0; i < n_clusters; ++i) {
        auto& c = state.cluster_tree.clusters[i];
        read_val(f, i32); c.id = i32;
        read_val(f, i32); c.center_node = i32;
        read_val(f, c.description_length);
        read_val(f, i32); c.parent_cluster_id = i32;
        uint32_t nm; read_val(f, nm);
        for (uint32_t j = 0; j < nm; ++j) {
            read_val(f, i32); c.observed_members.insert(i32);
        }
        uint32_t nc; read_val(f, nc);
        c.child_cluster_ids.resize(nc);
        for (uint32_t j = 0; j < nc; ++j) {
            read_val(f, i32); c.child_cluster_ids[j] = i32;
        }
    }

    uint32_t n_map;
    read_val(f, n_map);
    state.cluster_tree.node_to_cluster.clear();
    for (uint32_t i = 0; i < n_map; ++i) {
        int32_t node, cid;
        read_val(f, node); read_val(f, cid);
        state.cluster_tree.node_to_cluster[node] = cid;
    }

    state.cached_nuc_idx.clear();
    state.cached_valid.clear();
    if (version >= 2) {
        uint32_t cache_n, cache_m;
        read_val(f, cache_n);
        read_val(f, cache_m);
        if (cache_n > 0 && cache_m > 0) {
            size_t total = (size_t)cache_n * cache_m;
            state.cached_nuc_idx.resize(total);
            state.cached_valid.resize(total);
            f.read(reinterpret_cast<char*>(state.cached_nuc_idx.data()),
                   total * sizeof(int8_t));
            f.read(reinterpret_cast<char*>(state.cached_valid.data()),
                   total * sizeof(uint8_t));
        }
    }

    return f.good();
}

}  // namespace clnj
