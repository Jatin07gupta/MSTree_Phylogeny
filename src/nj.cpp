#include "nj.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace clnj {

LocalResult nj_local(
    const std::vector<int>& nodes,
    DistanceOracle& oracle,
    int m,
    int next_hidden_id,
    double min_branch_length,
    NegCounts* neg_cnt,
    std::vector<ZeroEdgeEntry>* zero_edge_log
) {
    LocalResult result;
    int n_local = (int)nodes.size();
    double _floor = std::max(0.0, min_branch_length);
    bool has_hidden = false;
    for (int nd : nodes) if (nd >= m) { has_hidden = true; break; }

    auto _dist = [&](int i, int j) {
        return oracle.dist_for_neighbourhood(i, j, has_hidden);
    };

    auto _count = [&](const char* key, double raw_val) {
        if (neg_cnt && raw_val < -1e-12)
            (*neg_cnt)[key] += 1;
    };

    if (n_local < 2) return result;

    if (n_local == 2) {
        double d_raw = _dist(nodes[0], nodes[1]);
        _count("initial_pair", d_raw);
        double d = std::max(_floor, d_raw);
        if (zero_edge_log && d < 1e-9) {
            ZeroEdgeEntry e;
            e.case_type = "initial_pair";
            e.a = nodes[0]; e.b = nodes[1]; e.w = d;
            e.d_raw = d_raw;
            e.reason = "only 2 nodes, oracle returned d~0";
            e.nodes.assign(nodes.begin(), nodes.end());
            zero_edge_log->push_back(std::move(e));
        }
        result.edges.emplace_back(nodes[0], nodes[1], d);
        return result;
    }

    int sz = n_local;
    int cap = n_local + n_local;
    std::vector<double> D(cap * cap, 0.0);
    auto Dref = [&](int i, int j) -> double& { return D[i * cap + j]; };

    for (int i = 0; i < n_local; ++i)
        for (int j = 0; j < n_local; ++j)
            if (i != j) {
                double d_raw = _dist(nodes[i], nodes[j]);
                _count("oracle_dist", d_raw);
                Dref(i, j) = std::max(0.0, d_raw);
            }

    std::vector<int> active(n_local);
    std::iota(active.begin(), active.end(), 0);
    std::vector<int> node_map(cap);
    for (int i = 0; i < n_local; ++i) node_map[i] = nodes[i];

    while ((int)active.size() > 2) {
        int n_active = (int)active.size();

        std::vector<double> row_sums(n_active, 0.0);
        for (int ai = 0; ai < n_active; ++ai)
            for (int aj = 0; aj < n_active; ++aj)
                if (active[ai] != active[aj])
                    row_sums[ai] += Dref(active[ai], active[aj]);

        double q_min = std::numeric_limits<double>::max();
        int best_ai = 0, best_aj = 1;
        for (int ai = 0; ai < n_active; ++ai)
            for (int aj = ai + 1; aj < n_active; ++aj) {
                double q = (n_active - 2) * Dref(active[ai], active[aj])
                           - row_sums[ai] - row_sums[aj];
                if (q < q_min) { q_min = q; best_ai = ai; best_aj = aj; }
            }

        int min_i = active[best_ai], min_j = active[best_aj];
        double d_ij = Dref(min_i, min_j);

        double sum_i = 0, sum_j = 0;
        for (int k : active) {
            sum_i += Dref(min_i, k);
            sum_j += Dref(min_j, k);
        }

        double delta_i_raw, delta_j_raw;
        if (n_active > 2) {
            delta_i_raw = 0.5 * d_ij + (sum_i - sum_j) / (2.0 * (n_active - 2));
        } else {
            delta_i_raw = 0.5 * d_ij;
        }
        delta_j_raw = d_ij - delta_i_raw;

        _count("delta_i", delta_i_raw);
        _count("delta_j", delta_j_raw);

        double delta_i, delta_j;
        if (delta_i_raw < _floor) {
            delta_i = _floor; delta_j = d_ij - _floor;
        } else if (delta_j_raw < _floor) {
            delta_j = _floor; delta_i = d_ij - _floor;
        } else {
            delta_i = delta_i_raw; delta_j = delta_j_raw;
        }

        int u_id = next_hidden_id++;
        int u_local = sz++;
        if (sz > cap) {
            int new_cap = cap * 2;
            std::vector<double> new_D(new_cap * new_cap, 0.0);
            for (int r = 0; r < cap; ++r)
                for (int c = 0; c < cap; ++c)
                    new_D[r * new_cap + c] = D[r * cap + c];
            D = std::move(new_D);
            cap = new_cap;
            node_map.resize(cap);
        }
        node_map[u_local] = u_id;

        int f_id = node_map[min_i], g_id = node_map[min_j];
        result.edges.emplace_back(f_id, u_id, delta_i);
        result.edges.emplace_back(g_id, u_id, delta_j);
        result.new_hidden_ids.push_back(u_id);

        if (zero_edge_log) {
            if (delta_i < 1e-9) {
                ZeroEdgeEntry e;
                e.case_type = "create_hidden";
                e.a = f_id; e.b = u_id; e.w = delta_i;
                e.d_ij = d_ij; e.delta_i_raw = delta_i_raw; e.delta_j_raw = delta_j_raw;
                e.reason = "delta_i=0 after redistribution";
                e.nodes.assign(nodes.begin(), nodes.end());
                zero_edge_log->push_back(std::move(e));
            }
            if (delta_j < 1e-9) {
                ZeroEdgeEntry e;
                e.case_type = "create_hidden";
                e.a = g_id; e.b = u_id; e.w = delta_j;
                e.d_ij = d_ij; e.delta_i_raw = delta_i_raw; e.delta_j_raw = delta_j_raw;
                e.reason = "delta_j=0 after redistribution";
                e.nodes.assign(nodes.begin(), nodes.end());
                zero_edge_log->push_back(std::move(e));
            }
        }

        for (int k : active) {
            if (k != min_i && k != min_j) {
                double d_uk_raw = 0.5 * (Dref(min_i, k) + Dref(min_j, k) - d_ij);
                _count("NJ_update_d_uk", d_uk_raw);
                double d_uk = std::max(0.0, d_uk_raw);
                Dref(u_local, k) = d_uk;
                Dref(k, u_local) = d_uk;
            }
        }

        active.erase(std::remove(active.begin(), active.end(), min_i), active.end());
        active.erase(std::remove(active.begin(), active.end(), min_j), active.end());
        active.push_back(u_local);
    }

    if ((int)active.size() == 2) {
        double d_raw = Dref(active[0], active[1]);
        _count("final_edge", d_raw);
        double d = std::max(_floor, d_raw);
        result.edges.emplace_back(node_map[active[0]], node_map[active[1]], d);
        if (zero_edge_log && d < 1e-9) {
            ZeroEdgeEntry e;
            e.case_type = "final_edge";
            e.a = node_map[active[0]]; e.b = node_map[active[1]]; e.w = d;
            e.d_raw = d_raw;
            e.reason = "final 2 nodes had d=0";
            e.nodes.assign(nodes.begin(), nodes.end());
            zero_edge_log->push_back(std::move(e));
        }
    }

    return result;
}

}  // namespace clnj
