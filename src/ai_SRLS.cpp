/*
 * SRLS - Swap Relaxation-based Local Search for Latin Square Completion
 * Based on Xie et al., IJCAI-24
 * Usage: srls <seed> <instance> <output> [time_limit]
 *
 * Build: g++ -O3 -march=native -std=c++23 -o srls.exe srls_new.cpp
 */
 #include "utils.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cstring>

using namespace std;
using namespace std::chrono;

// ============== Minimal mdspan-like 2D view (std::mdspan unavailable in this toolchain)
// ==============
template <typename T> class mdspan2d {
    T *ptr_;
    size_t rows_, cols_;

  public:
    mdspan2d() : ptr_(nullptr), rows_(0), cols_(0) {}
    mdspan2d(T *p, size_t r, size_t c) : ptr_(p), rows_(r), cols_(c) {}
    T &operator[](size_t i, size_t j) const { return ptr_[i * cols_ + j]; }
    size_t extent(size_t d) const { return d == 0 ? rows_ : cols_; }
    T *data_handle() const { return ptr_; }
};

// ============== Owning 2D matrix: vector storage + mdspan2d view ==============
template <typename T> class Mat {
    vector<T> buf_;
    size_t rows_ = 0, cols_ = 0;

  public:
    Mat() = default;
    void resize(size_t r, size_t c, T v = {}) {
        rows_ = r;
        cols_ = c;
        buf_.assign(r * c, v);
    }
    mdspan2d<T> span() { return {buf_.data(), rows_, cols_}; }
    T &operator[](size_t i, size_t j) { return buf_[i * cols_ + j]; }
    const T &operator[](size_t i, size_t j) const { return buf_[i * cols_ + j]; }
    void fill(T v = {}) { std::fill(buf_.begin(), buf_.end(), v); }
    T *data() { return buf_.data(); }
    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
};

static constexpr int ALPHA_NUM = 2, ALPHA_DEN = 5;
static constexpr int RT0 = 10, RT_UB = 15, ACCU_UB = 1000, TABU_P = 10;
static constexpr int CHECK_INTERVAL = 1024;

// ============== Instance ==============
static int n, ncells;
static Mat<uint8_t> dom;         // dom[v, k]: whether color k is in domain of cell v
static vector<int> dsz;          // dsz[v]: domain size of cell v
static vector<uint8_t> is_fixed; // is_fixed[v]: whether cell v is fixed
static vector<int> crow, ccol;   // crow[v], ccol[v]: row and column of cell v

// ============== Solution ==============
static vector<int> clr, bclr; // current / best color assignment
static Mat<int> cc;           // cc[col, color]: count of color in column
static Mat<int> crx;          // crx[col, color]: XOR of row indices
static int CE, DV, bCE, bDV;
static Mat<int> tabu;        // tabu[cell, color]: iteration expiry
static Mat<int> ruf;         // ruf[row, idx]: unfixed cells per row
static vector<int> rufc;     // rufc[row]: count of unfixed cells in row
static vector<int8_t> dv_of; // dv_of[v]: 1 if clr[v] not in dom[v]

// ============== Incremental conflict set ==============
static vector<int> clist, cpos;
static int csize;
static vector<uint8_t> cin_f;

static inline void conf_add(int v) {
    if (cin_f[v])
        return;
    cin_f[v] = 1;
    cpos[v] = csize;
    clist[csize++] = v;
}
static inline void conf_rem(int v) {
    if (!cin_f[v])
        return;
    cin_f[v] = 0;
    int p = cpos[v];
    int last = clist[--csize];
    clist[p] = last;
    cpos[last] = p;
}
static inline void conf_update(int v) {
    if (is_fixed[v])
        return;
    if (cc[ccol[v], clr[v]] >= 2 || !dom[v, clr[v]])
        conf_add(v);
    else
        conf_rem(v);
}
static void build_conf() {
    csize = 0;
    fill(cin_f.begin(), cin_f.end(), (uint8_t)0);
    for (int v = 0; v < ncells; v++) {
        if (is_fixed[v])
            continue;
        if (cc[ccol[v], clr[v]] >= 2 || !dom[v, clr[v]]) {
            cin_f[v] = 1;
            cpos[v] = csize;
            clist[csize++] = v;
        }
    }
}

// ============== Timer ==============
static high_resolution_clock::time_point T0;
static double TL = 1000.0;
static inline double timer() { return duration<double>(high_resolution_clock::now() - T0).count(); }

// ============== Allocate all dynamic arrays ==============
static void allocate_arrays() {
    dom.resize(ncells, n + 1, 0);
    dsz.assign(ncells, 0);
    is_fixed.assign(ncells, 0);
    clr.assign(ncells, 0);
    bclr.resize(ncells);
    crow.resize(ncells);
    ccol.resize(ncells);
    cc.resize(n, n + 1, 0);
    crx.resize(n, n + 1, 0);
    tabu.resize(ncells, n + 1, 0);
    ruf.resize(n, n, 0);
    rufc.assign(n, 0);
    dv_of.assign(ncells, 0);
    clist.resize(ncells);
    cpos.resize(ncells);
    cin_f.assign(ncells, 0);
}

// ============== Input ==============
static void read_instance(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "open fail\n");
        exit(1);
    }
    n = 0;
    ncells = 0;
    char buf[1 << 16];
    int fv = -1;
    vector<int> fc;
    bool allocated = false;

    auto ensure_alloc = [&]() {
        if (!allocated && n > 0) {
            ncells = n * n;
            allocate_arrays();
            allocated = true;
        }
    };

    auto flush = [&]() {
        ensure_alloc();
        if (n > 0 && fv >= 0 && fv < ncells) {
            for (int k : fc)
                if (k >= 1 && k <= n)
                    dom[fv, k] = 1;
            int c = 0;
            for (int k = 1; k <= n; k++)
                c += dom[fv, k];
            dsz[fv] = c;
            // Don't set is_fixed here; let reduce() handle propagation
            if (c == 1) {
                for (int k = 1; k <= n; k++)
                    if (dom[fv, k]) {
                        clr[fv] = k;
                        break;
                    }
            }
        }
        fv = -1;
        fc.clear();
    };

    auto si = [](char *p, vector<int> &o) {
        while (*p) {
            while (*p && (*p < '0' || *p > '9'))
                ++p;
            if (*p >= '0' && *p <= '9') {
                o.push_back(atoi(p));
                while (*p >= '0' && *p <= '9')
                    ++p;
            }
        }
    };

    while (fgets(buf, sizeof(buf), fp)) {
        char *p = buf;
        while (*p == ' ' || *p == '\t')
            ++p;
        if (!*p || *p == '\n' || *p == '\r')
            continue;
        if (*p == 'p') {
            flush();
            vector<int> v;
            si(p + 1, v);
            if (!v.empty()) {
                n = (int)round(sqrt((double)v[0]));
                ncells = n * n;
            }
        } else if (*p == 'e' || *p == 'c') {
            flush();
        } else if (*p == 'f') {
            flush();
            vector<int> v;
            si(p + 1, v);
            if (!v.empty()) {
                fv = v[0] - 1;
                for (int i = 1; i < (int)v.size(); i++)
                    fc.push_back(v[i]);
            }
        } else if (*p >= '0' && *p <= '9' && fv >= 0) {
            vector<int> v;
            si(p, v);
            for (int x : v)
                fc.push_back(x);
        }
    }
    flush();
    fclose(fp);
    ensure_alloc();
    if (n <= 0) {
        fprintf(stderr, "parse fail\n");
        exit(1);
    }
    for (int v = 0; v < ncells; v++) {
        if (dsz[v] == 0) {
            for (int k = 1; k <= n; k++)
                dom[v, k] = 1;
            dsz[v] = n;
        }
        crow[v] = v / n;
        ccol[v] = v % n;
    }
}

// ============== Reduction ==============
struct BS {
    unsigned long long w[2] = {};
    void set(int i) { w[i >> 6] |= 1ULL << (i & 63); }
    bool get(int i) const { return (w[i >> 6] >> (i & 63)) & 1; }
    BS operator|(const BS &o) const {
        BS r;
        r.w[0] = w[0] | o.w[0];
        r.w[1] = w[1] | o.w[1];
        return r;
    }
    bool operator==(const BS &o) const { return w[0] == o.w[0] && w[1] == o.w[1]; }
    int count() const { return __builtin_popcountll(w[0]) + __builtin_popcountll(w[1]); }
};
static bool reduce_group(vector<int> &cells) {
    int m = (int)cells.size();
    if (m <= 1)
        return false;
    if (m > 16)
        return false; // cap exponential naked tuple search
    bool changed = false;
    for (int i = 0; i < m; i++) {
        int v = cells[i];
        if (dsz[v] != 1)
            continue;
        int k = 0;
        for (int c = 1; c <= n; c++)
            if (dom[v, c]) {
                k = c;
                break;
            }
        for (int j = 0; j < m; j++) {
            if (j == i)
                continue;
            int u = cells[j];
            if (dom[u, k]) {
                dom[u, k] = 0;
                dsz[u]--;
                changed = true;
            }
        }
        cells.erase(cells.begin() + i);
        return changed | reduce_group(cells);
    }
    struct E {
        BS b;
        vector<int> idx;
    };
    vector<vector<E>> f(m + 1);
    for (int i = 0; i < m; i++) {
        BS vi;
        for (int c = 1; c <= n; c++)
            if (dom[cells[i], c])
                vi.set(c);
        for (int sz = min(i, m - 2); sz >= 1; sz--) {
            for (auto &e : f[sz]) {
                BS u = e.b | vi;
                if (u.count() == sz + 1) {
                    vector<bool> ins(m, false);
                    for (int x : e.idx)
                        ins[x] = true;
                    ins[i] = true;
                    for (int j = 0; j < m; j++) {
                        if (ins[j])
                            continue;
                        int cv = cells[j];
                        for (int c = 1; c <= n; c++)
                            if (u.get(c) && dom[cv, c]) {
                                dom[cv, c] = 0;
                                dsz[cv]--;
                                changed = true;
                            }
                    }
                    if (changed) {
                        vector<int> rest, hall;
                        for (int j = 0; j < m; j++) {
                            if (ins[j])
                                hall.push_back(cells[j]);
                            else
                                rest.push_back(cells[j]);
                        }
                        return true | reduce_group(rest) | reduce_group(hall);
                    }
                }
                bool found = false;
                for (auto &e2 : f[sz + 1])
                    if (e2.b == u) {
                        found = true;
                        break;
                    }
                if (!found) {
                    E ne;
                    ne.b = u;
                    ne.idx = e.idx;
                    ne.idx.push_back(i);
                    f[sz + 1].push_back(ne);
                }
            }
        }
        E se;
        se.b = vi;
        se.idx = {i};
        f[1].push_back(se);
    }
    return changed;
}
static bool reduce() {
    bool any = true;
    while (any) {
        any = false;
        for (int v = 0; v < ncells; v++) {
            if (is_fixed[v] || dsz[v] != 1)
                continue;
            is_fixed[v] = 1;
            int k = 0;
            for (int c = 1; c <= n; c++)
                if (dom[v, c]) {
                    k = c;
                    break;
                }
            clr[v] = k;
            int r = v / n, col = v % n;
            for (int j = 0; j < n; j++) {
                int u = r * n + j;
                if (u != v && dom[u, k]) {
                    dom[u, k] = 0;
                    if (--dsz[u] == 0)
                        return false;
                }
            }
            for (int i = 0; i < n; i++) {
                int u = i * n + col;
                if (u != v && dom[u, k]) {
                    dom[u, k] = 0;
                    if (--dsz[u] == 0)
                        return false;
                }
            }
            any = true;
        }
        for (int r = 0; r < n; r++) {
            for (int k = 1; k <= n; k++) {
                int cand = -1, cnt = 0;
                bool done = false;
                for (int j = 0; j < n; j++) {
                    int v = r * n + j;
                    if (is_fixed[v] && clr[v] == k) {
                        done = true;
                        break;
                    }
                    if (!is_fixed[v] && dom[v, k]) {
                        cand = v;
                        cnt++;
                    }
                }
                if (done || cnt != 1)
                    continue;
                for (int c = 1; c <= n; c++)
                    dom[cand, c] = 0;
                dom[cand, k] = 1;
                dsz[cand] = 1;
                any = true;
            }
        }
        for (int col = 0; col < n; col++) {
            for (int k = 1; k <= n; k++) {
                int cand = -1, cnt = 0;
                bool done = false;
                for (int i = 0; i < n; i++) {
                    int v = i * n + col;
                    if (is_fixed[v] && clr[v] == k) {
                        done = true;
                        break;
                    }
                    if (!is_fixed[v] && dom[v, k]) {
                        cand = v;
                        cnt++;
                    }
                }
                if (done || cnt != 1)
                    continue;
                for (int c = 1; c <= n; c++)
                    dom[cand, c] = 0;
                dom[cand, k] = 1;
                dsz[cand] = 1;
                any = true;
            }
        }
        for (int r = 0; r < n; r++) {
            vector<int> cells;
            for (int j = 0; j < n; j++) {
                int v = r * n + j;
                if (!is_fixed[v])
                    cells.push_back(v);
            }
            if (reduce_group(cells))
                any = true;
        }
        for (int col = 0; col < n; col++) {
            vector<int> cells;
            for (int i = 0; i < n; i++) {
                int v = i * n + col;
                if (!is_fixed[v])
                    cells.push_back(v);
            }
            if (reduce_group(cells))
                any = true;
        }
    }
    return true;
}

// ============== Init ==============
static void init_solution() {
    for (int r = 0; r < n; r++) {
        rufc[r] = 0;
        vector<uint8_t> used(n + 1, 0);
        for (int j = 0; j < n; j++) {
            int v = r * n + j;
            if (is_fixed[v])
                used[clr[v]] = 1;
            else
                ruf[r, rufc[r]++] = v;
        }
        // Domain-aware init (Algorithm 2): sort by dsz, assign from D(v)∩Available
        int uc = rufc[r];
        vector<int> sorted_uf(uc);
        for (int i = 0; i < uc; i++)
            sorted_uf[i] = ruf[r, i];
        for (int i = 0; i < uc - 1; i++) {
            int mi = i;
            for (int j = i + 1; j < uc; j++)
                if (dsz[sorted_uf[j]] < dsz[sorted_uf[mi]])
                    mi = j;
            if (mi != i)
                swap(sorted_uf[i], sorted_uf[mi]);
        }
        vector<uint8_t> avail(n + 1, 0);
        for (int k = 1; k <= n; k++)
            avail[k] = !used[k];
        for (int i = 0; i < uc; i++) {
            int v = sorted_uf[i];
            vector<int> cands, others;
            cands.reserve(n);
            others.reserve(n);
            for (int k = 1; k <= n; k++)
                if (avail[k]) {
                    if (dom[v, k])
                        cands.push_back(k);
                    else
                        others.push_back(k);
                }
            int chosen;
            if (!cands.empty())
                chosen = cands[rng32() % (unsigned)cands.size()];
            else
                chosen = others[rng32() % (unsigned)others.size()];
            clr[v] = chosen;
            avail[chosen] = 0;
        }
    }
    cc.fill(0);
    crx.fill(0);
    for (int v = 0; v < ncells; v++) {
        int r = crow[v], j = ccol[v], c = clr[v];
        cc[j, c]++;
        crx[j, c] ^= r;
    }
    CE = 0;
    for (int j = 0; j < n; j++)
        for (int k = 1; k <= n; k++)
            CE += cc[j, k] * (cc[j, k] - 1) / 2;
    DV = 0;
    for (int v = 0; v < ncells; v++) {
        dv_of[v] = (!is_fixed[v] && !dom[v, clr[v]]) ? 1 : 0;
        DV += dv_of[v];
    }
    build_conf();
}

// ============== Solve ==============
static bool solve() {
    bCE = CE;
    bDV = DV;
    bclr = clr;
    if (CE == 0 && DV == 0)
        return true;

    int rt = RT0, accu = 0, iter = 0;
    tabu.fill(0);

    while (true) {
        if ((iter & (CHECK_INTERVAL - 1)) == 0 && timer() >= TL)
            break;
        if (csize == 0)
            return CE == 0 && DV == 0;
        ++iter;

        // ---- Move evaluation ----
        int bn_u = -1, bn_v = -1, bn_f1 = 1 << 29, bn_f2 = 1 << 29, bn_c = 0;
        int bt_u = -1, bt_v = -1, bt_f1 = 1 << 29, bt_f2 = 1 << 29, bt_c = 0;
        int nr_u = -1, nr_v = -1, nr_f1 = 1 << 29, nr_f2 = 1 << 29;
        int rd_u = -1, rd_v = -1, rd_f1 = 1 << 29, rd_f2 = 1 << 29;
        int nr_cnt = 1, rd_cnt = 1; // reservoir sampling counters

        for (int ci = 0; ci < csize; ci++) {
            const int u = clist[ci];
            const int r = crow[u], ju = ccol[u], cu = clr[u];
            int pbn_u = -1, pbn_v = -1, pbn_f1 = 1 << 29, pbn_f2 = 1 << 29, pbn_c = 0;
            int pbt_u = -1, pbt_v = -1, pbt_f1 = 1 << 29, pbt_f2 = 1 << 29, pbt_c = 0;
            const int cc_ju_cu = cc[ju, cu];
            const int du_cu = dv_of[u];

            const int rc = rufc[r];
            for (int vi = 0; vi < rc; vi++) {
                const int v = ruf[r, vi];
                if (v == u)
                    continue;
                const int jv = ccol[v], cv = clr[v];
                if (cu == cv)
                    continue;

                const int f1 = (cc[ju, cv] - cc_ju_cu + 1) + (cc[jv, cu] - cc[jv, cv] + 1);
                const int f2 = (!dom[u, cv] - du_cu) + (!dom[v, cu] - dv_of[v]);
                const bool tb = (tabu[u, cv] > iter) || (tabu[v, cu] > iter);

                if (tb) {
                    if (f1 < bt_f1 || (f1 == bt_f1 && f2 < bt_f2)) {
                        bt_u = u;
                        bt_v = v;
                        bt_f1 = f1;
                        bt_f2 = f2;
                        bt_c = 1;
                    } else if (f1 == bt_f1 && f2 == bt_f2) {
                        if (rng32() % (unsigned)(++bt_c) == 0) {
                            bt_u = u;
                            bt_v = v;
                        }
                    }
                    if (f1 < pbt_f1 || (f1 == pbt_f1 && f2 < pbt_f2)) {
                        pbt_u = u;
                        pbt_v = v;
                        pbt_f1 = f1;
                        pbt_f2 = f2;
                        pbt_c = 1;
                    } else if (f1 == pbt_f1 && f2 == pbt_f2) {
                        if (rng32() % (unsigned)(++pbt_c) == 0) {
                            pbt_u = u;
                            pbt_v = v;
                        }
                    }
                } else {
                    if (f1 < bn_f1 || (f1 == bn_f1 && f2 < bn_f2)) {
                        bn_u = u;
                        bn_v = v;
                        bn_f1 = f1;
                        bn_f2 = f2;
                        bn_c = 1;
                    } else if (f1 == bn_f1 && f2 == bn_f2) {
                        if (rng32() % (unsigned)(++bn_c) == 0) {
                            bn_u = u;
                            bn_v = v;
                        }
                    }
                    if (f1 < pbn_f1 || (f1 == pbn_f1 && f2 < pbn_f2)) {
                        pbn_u = u;
                        pbn_v = v;
                        pbn_f1 = f1;
                        pbn_f2 = f2;
                        pbn_c = 1;
                    } else if (f1 == pbn_f1 && f2 == pbn_f2) {
                        if (rng32() % (unsigned)(++pbn_c) == 0) {
                            pbn_u = u;
                            pbn_v = v;
                        }
                    }
                }
            }
            if (pbt_u >= 0 && rng32() % (unsigned)(rd_cnt++) == 0) {
                rd_u = pbt_u;
                rd_v = pbt_v;
                rd_f1 = pbt_f1;
                rd_f2 = pbt_f2;
            }
            if (pbn_u >= 0 && rng32() % (unsigned)(nr_cnt++) == 0) {
                nr_u = pbn_u;
                nr_v = pbn_v;
                nr_f1 = pbn_f1;
                nr_f2 = pbn_f2;
            }
        }

        // Exploration fallback
        if (bt_f1 > 0 && rd_u >= 0) {
            bt_u = rd_u;
            bt_v = rd_v;
            bt_f1 = rd_f1;
            bt_f2 = rd_f2;
        }
        if (bn_f1 > 0 && nr_u >= 0) {
            bn_u = nr_u;
            bn_v = nr_v;
            bn_f1 = nr_f1;
            bn_f2 = nr_f2;
        }

        // Aspiration + selection (lexicographic: (f1, f2) comparison)
        int su, sv, sf1, sf2;
        bool tb_better = (bt_f1 < bn_f1) || (bt_f1 == bn_f1 && bt_f2 < bn_f2);
        bool aspiration = (CE + bt_f1 < bCE) || (CE + bt_f1 == bCE && DV + bt_f2 < bDV);
        bool use_tb = bt_u >= 0 && tb_better && aspiration;
        if (use_tb) {
            su = bt_u;
            sv = bt_v;
            sf1 = bt_f1;
            sf2 = bt_f2;
        } else if (bn_u >= 0) {
            su = bn_u;
            sv = bn_v;
            sf1 = bn_f1;
            sf2 = bn_f2;
        } else
            continue;

        // ---- Execute swap ----
        const int u = su, v = sv, cu = clr[u], cv = clr[v];
        const int ju = ccol[u], jv = ccol[v], ru = crow[u];

        // Tabu (OLD CE for tenure)
        int ten1 = CE * ALPHA_NUM / ALPHA_DEN + 1 + (int)(rng32() % TABU_P);
        int ten2 = CE * ALPHA_NUM / ALPHA_DEN + 1 + (int)(rng32() % TABU_P);
        tabu[u, cu] = iter + ten1;
        tabu[v, cv] = iter + ten2;

        // Update column counts and XOR
        cc[ju, cu]--;
        crx[ju, cu] ^= ru;
        cc[ju, cv]++;
        crx[ju, cv] ^= ru;
        cc[jv, cv]--;
        crx[jv, cv] ^= ru;
        cc[jv, cu]++;
        crx[jv, cu] ^= ru;
        CE += sf1;
        DV += sf2;
        clr[u] = cv;
        clr[v] = cu;
        dv_of[u] = !dom[u, cv];
        dv_of[v] = !dom[v, cu];

        // ---- Incremental conflict set update ----
        // (ju, cu) decreased: remaining cell might leave conf
        if (cc[ju, cu] == 1) {
            int w = crx[ju, cu] * n + ju;
            if (!is_fixed[w])
                conf_update(w);
        }
        // (ju, cv) increased: other cell might join conf
        if (cc[ju, cv] == 2) {
            int w = (crx[ju, cv] ^ ru) * n + ju;
            if (!is_fixed[w])
                conf_update(w);
        }
        // (jv, cv) decreased
        if (cc[jv, cv] == 1) {
            int w = crx[jv, cv] * n + jv;
            if (!is_fixed[w])
                conf_update(w);
        }
        // (jv, cu) increased
        if (cc[jv, cu] == 2) {
            int w = (crx[jv, cu] ^ ru) * n + jv;
            if (!is_fixed[w])
                conf_update(w);
        }
        // Swapped cells
        conf_update(u);
        conf_update(v);

        if (CE == 0 && DV == 0)
            return true;

        if (CE < bCE || (CE == bCE && DV <= bDV)) {
            bCE = CE;
            bDV = DV;
            bclr = clr;
        } else if (CE - bCE > rt) {
            // Adaptive restart from best
            clr = bclr;
            cc.fill(0);
            crx.fill(0);
            for (int vv = 0; vv < ncells; vv++) {
                cc[ccol[vv], clr[vv]]++;
                crx[ccol[vv], clr[vv]] ^= crow[vv];
            }
            CE = bCE;
            DV = bDV;
            for (int vv = 0; vv < ncells; vv++)
                dv_of[vv] = (!is_fixed[vv] && !dom[vv, clr[vv]]) ? 1 : 0;
            tabu.fill(0);
            iter = 0;
            build_conf();
            if (rt < RT_UB) {
                if (++accu == ACCU_UB) {
                    accu = 0;
                    rt++;
                }
            }
        }
    }
    return false;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <seed> <inst> <out> [tl]\n", argv[0]);
        return 1;
    }
    rng_seed((unsigned)atoi(argv[1]));
    if (argc >= 5)
        TL = atof(argv[4]);
    T0 = high_resolution_clock::now();
    read_instance(argv[2]);
    reduce();

    int nf = 0;
    for (int v = 0; v < ncells; v++)
        nf += is_fixed[v];
    int sum_dsz = 0, min_dsz = n, max_dsz = 0, uf_cnt = 0;
    for (int v = 0; v < ncells; v++)
        if (!is_fixed[v]) {
            sum_dsz += dsz[v];
            if (dsz[v] < min_dsz)
                min_dsz = dsz[v];
            if (dsz[v] > max_dsz)
                max_dsz = dsz[v];
            uf_cnt++;
        }
    fprintf(stderr, "n=%d fixed=%d unfixed=%d dsz[min=%d avg=%.1f max=%d]\n", n, nf, uf_cnt,
            min_dsz, uf_cnt ? (double)sum_dsz / uf_cnt : 0, max_dsz);

    init_solution();
    fprintf(stderr, "init CE=%d DV=%d\n", CE, DV);
    bool ok = solve();
    if (!ok)
        clr = bclr;
    FILE *out = fopen(argv[3], "w");
    for (int v = 0; v < ncells; v++)
        fprintf(out, "%d\n", clr[v]);
    fclose(out);
    double t = timer();
    if (ok)
        fprintf(stderr, "SOLVED seed=%s time=%.3fs\n", argv[1], t);
    else
        fprintf(stderr, "TIMEOUT seed=%s time=%.3fs best_CE=%d best_DV=%d\n", argv[1], t, bCE, bDV);
	    const string logFile = "results.csv";
    size_t len = strlen(argv[2]);
    while (len > 0 && (argv[2][len - 1] != '/' && argv[2][len - 1] != '\\')) {
        len--;
    }
    ofstream csvFile(logFile, ios::app);
    csvFile << now_str() << ", " << (argv[2] + len) << ", "
            << "ai_SRLS, " << argv[1] << ", " << t << ", " << (ok ? 0 : bCE) << "\n";
    return ok ? 0 : 1;
}
