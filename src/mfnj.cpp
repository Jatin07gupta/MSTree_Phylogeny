#include "mfnj.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>

namespace clnj {

static std::vector<std::set<int>> find_tied_groups(
    const std::vector<std::tuple<double, int, int>>& q_list,
    double tol = TIE_TOL
) {
    if (q_list.empty()) return {};

    double q_min = std::get<0>(q_list[0]);
    std::vector<std::pair<int,int>> tied_pairs;
    for (auto& [qv, ai, aj] : q_list) {
        if (qv > q_min + tol) break;
        tied_pairs.emplace_back(ai, aj);
    }

    if (tied_pairs.size() == 1)
        return {{tied_pairs[0].first, tied_pairs[0].second}};

    std::map<int, std::set<int>> adj;
    for (auto& [a, b] : tied_pairs) {
        adj[a].insert(b);
        adj[b].insert(a);
    }

    std::set<int> visited;
    std::vector<std::set<int>> groups;
    for (auto& [start, _] : adj) {
        if (visited.count(start)) continue;
        std::set<int> component;
        std::vector<int> stack = {start};
        while (!stack.empty()) {
            int node = stack.back(); stack.pop_back();
            if (visited.count(node)) continue;
            visited.insert(node);
            component.insert(node);
            for (int nb : adj[node])
                if (!visited.count(nb)) stack.push_back(nb);
        }
        groups.push_back(std::move(component));
    }
    return groups;
}

LocalResult mfnj_local(
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

    int cap = n_local + n_local;
    int sz = n_local;
    std::vector<double> D(cap * cap, 0.0);
    auto Dref = [&](int i, int j) -> double& { return D[i * cap + j]; };

    auto grow_D = [&]() {
        int new_cap = cap * 2;
        std::vector<double> new_D(new_cap * new_cap, 0.0);
        for (int r = 0; r < cap; ++r)
            for (int c = 0; c < cap; ++c)
                new_D[r * new_cap + c] = D[r * cap + c];
        D = std::move(new_D);
        cap = new_cap;
    };

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

    while ((int)active.size() > 1) {
        int n_act = (int)active.size();

        if (n_act == 2) {
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
            break;
        }

        if (n_act == 3) {
            int i0 = active[0], i1 = active[1], i2 = active[2];
            double L0 = std::max(_floor, (Dref(i0,i1) + Dref(i0,i2) - Dref(i1,i2)) / 2.0);
            double L1 = std::max(_floor, (Dref(i0,i1) + Dref(i1,i2) - Dref(i0,i2)) / 2.0);
            double L2 = std::max(_floor, (Dref(i0,i2) + Dref(i1,i2) - Dref(i0,i1)) / 2.0);

            int u_id = next_hidden_id++;
            result.new_hidden_ids.push_back(u_id);
            result.edges.emplace_back(node_map[i0], u_id, L0);
            result.edges.emplace_back(node_map[i1], u_id, L1);
            result.edges.emplace_back(node_map[i2], u_id, L2);

            if (zero_edge_log) {
                struct { const char* lbl; int lid; double lval; } arms[] = {
                    {"i0", i0, L0}, {"i1", i1, L1}, {"i2", i2, L2}
                };
                for (auto& arm : arms) {
                    if (arm.lval < 1e-9) {
                        ZeroEdgeEntry e;
                        e.case_type = "final_star";
                        e.a = node_map[arm.lid]; e.b = u_id; e.w = arm.lval;
                        e.reason = std::string("final 3 nodes, branch ") + arm.lbl + "=0";
                        e.nodes.assign(nodes.begin(), nodes.end());
                        zero_edge_log->push_back(std::move(e));
                    }
                }
            }
            break;
        }

        // Q-matrix
        std::vector<double> row_sums(n_act, 0.0);
        for (int pi = 0; pi < n_act; ++pi)
            for (int pj = 0; pj < n_act; ++pj)
                if (active[pi] != active[pj])
                    row_sums[pi] += Dref(active[pi], active[pj]);

        std::vector<std::tuple<double, int, int>> q_list;
        for (int pi = 0; pi < n_act; ++pi)
            for (int pj = pi + 1; pj < n_act; ++pj) {
                double q = (n_act - 2) * Dref(active[pi], active[pj])
                           - row_sums[pi] - row_sums[pj];
                q_list.emplace_back(q, active[pi], active[pj]);
            }
        std::sort(q_list.begin(), q_list.end());

        auto groups = find_tied_groups(q_list);

        std::vector<std::pair<std::set<int>, int>> processed_groups;
        for (auto& group : groups) {
            std::vector<int> I(group.begin(), group.end());
            std::sort(I.begin(), I.end());
            int P = (int)I.size();
            if (P < 2) continue;

            std::vector<int> complement;
            for (int k : active)
                if (!group.count(k)) complement.push_back(k);

            double R_II = 0.0;
            for (int ai = 0; ai < P; ++ai)
                for (int bi = ai + 1; bi < P; ++bi)
                    R_II += Dref(I[ai], I[bi]);

            int u_id = next_hidden_id++;
            result.new_hidden_ids.push_back(u_id);

            if (sz >= cap) { grow_D(); node_map.resize(cap); }
            int u_local = sz++;
            node_map[u_local] = u_id;

            std::map<int, double> branch_lengths;
            if (!complement.empty()) {
                double R_IIc = 0.0;
                for (int i : I)
                    for (int k : complement)
                        R_IIc += Dref(i, k);

                for (int i : I) {
                    double R_iIc = 0.0;
                    for (int k : complement) R_iIc += Dref(i, k);
                    double L_raw = R_II / (P * (P - 1))
                                 + R_iIc / (n_act - P)
                                 - R_IIc / (P * (n_act - P));
                    _count("mfnj_branch", L_raw);
                    branch_lengths[i] = std::max(_floor, L_raw);
                }
            } else {
                if (P == 2) {
                    double d_val = Dref(I[0], I[1]);
                    branch_lengths[I[0]] = std::max(_floor, d_val / 2.0);
                    branch_lengths[I[1]] = std::max(_floor, d_val / 2.0);
                } else {
                    for (int i : I) {
                        double R_iI = 0.0;
                        for (int ip : I)
                            if (ip != i) R_iI += Dref(i, ip);
                        double L_raw = R_iI / (P - 2) - R_II / ((P - 1) * (P - 2));
                        _count("mfnj_branch_final", L_raw);
                        branch_lengths[i] = std::max(_floor, L_raw);
                    }
                }
            }

            for (int i : I) {
                result.edges.emplace_back(node_map[i], u_id, branch_lengths[i]);
                if (zero_edge_log && branch_lengths[i] < 1e-9) {
                    ZeroEdgeEntry e;
                    e.case_type = (P > 2) ? "mfnj_polytomy" : "create_hidden";
                    e.a = node_map[i]; e.b = u_id; e.w = branch_lengths[i];
                    e.reason = "MFNJ group size " + std::to_string(P)
                             + ", branch for local idx " + std::to_string(i);
                    e.nodes.assign(nodes.begin(), nodes.end());
                    zero_edge_log->push_back(std::move(e));
                }
            }

            for (int k : complement) {
                double R_Ik = 0.0;
                for (int i : I) R_Ik += Dref(i, k);
                double D_uk_raw = R_Ik / P - R_II / (P * (P - 1));
                _count("mfnj_update_d_uk", D_uk_raw);
                double D_uk = std::max(0.0, D_uk_raw);
                Dref(u_local, k) = D_uk;
                Dref(k, u_local) = D_uk;
            }

            for (int i : I)
                active.erase(std::remove(active.begin(), active.end(), i), active.end());
            active.push_back(u_local);
            processed_groups.push_back({group, u_local});
        }

        if (processed_groups.size() > 1) {
            for (size_t ai = 0; ai < processed_groups.size(); ++ai) {
                for (size_t bi = ai + 1; bi < processed_groups.size(); ++bi) {
                    int ua = processed_groups[ai].second;
                    int ub = processed_groups[bi].second;
                    auto& ga = processed_groups[ai].first;
                    auto& gb = processed_groups[bi].first;
                    int Pa = (int)ga.size(), Pb = (int)gb.size();

                    double R_IJ = 0.0;
                    for (int i : ga) for (int j : gb) R_IJ += Dref(i, j);

                    double term_II = 0.0;
                    if (Pa > 1) {
                        std::vector<int> ga_v(ga.begin(), ga.end());
                        for (int x = 0; x < Pa; ++x)
                            for (int y = x + 1; y < Pa; ++y)
                                term_II += Dref(ga_v[x], ga_v[y]);
                        term_II /= Pa * (Pa - 1);
                    }
                    double term_JJ = 0.0;
                    if (Pb > 1) {
                        std::vector<int> gb_v(gb.begin(), gb.end());
                        for (int x = 0; x < Pb; ++x)
                            for (int y = x + 1; y < Pb; ++y)
                                term_JJ += Dref(gb_v[x], gb_v[y]);
                        term_JJ /= Pb * (Pb - 1);
                    }
                    double D_uv = std::max(0.0, R_IJ / (Pa * Pb) - term_II - term_JJ);
                    Dref(ua, ub) = D_uv;
                    Dref(ub, ua) = D_uv;
                }
            }
        }
    }

    return result;
}

}  // namespace clnj
