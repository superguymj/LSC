#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#define rep(i, x, y) for (int i = x; i <= y; i++)
#define repd(i, x, y) for (int i = x; i >= y; i--)
#define rep(i, x, y) for (int i = x; i <= y; i++)
#define repd(i, x, y) for (int i = x; i >= y; i--)

#define MESSAGE

using namespace std;
using namespace __gnu_pbds;

namespace Constants {

using uint = unsigned int;
using i64 = long long;
using ui64 = unsigned long long;

const int inf = 0x3f3f3f3f;
const i64 inf64 = 1e18;
const float infF32 = 1e100, eps32 = 1e-4;
const double infF64 = 1e100, eps64 = 1e-8;
template <typename T>
bool eqls(T a, T b) {
    return a == b;
}
bool eqls(float a, float b) {
    return std::abs(a - b) <= eps32;
}
bool eqls(double a, double b) {
    return std::abs(a - b) <= eps64;
}
int getinf(int z) {
    return inf;
}
i64 getinf(i64 z) {
    return inf64;
}
float getinf(float z) {
    return infF32;
}
double getinf(double z) {
    return infF64;
}
}  // namespace Constants

namespace MaxFlow {
using namespace Constants;

template <typename FlowType>
struct FlowGraph {
    int N, S, T;
    struct Edge {
        int to;
        FlowType f;
    };
    std::vector<Edge> e;
    std::vector<std::vector<int>> eid;
    FlowGraph(int n, int s, int t) : N(n), S(s), T(t), eid(N) {}
    std::vector<std::tuple<int, int, FlowType>> edgeList;
    int addpt() {
        eid.emplace_back();
        return N++;
    }
    void preAddEdge(int u, int v, FlowType f) {
        edgeList.emplace_back(u, v, f);
    }
    void addedge(int u, int v, FlowType f) {
        eid[u].emplace_back(e.size()), e.push_back({v, f});
        eid[v].emplace_back(e.size()), e.push_back({u, 0});
    }
    void build() {
        sort(edgeList.begin(), edgeList.end());
        FlowType f = 0;
        std::vector<int> degree(N);
        for (int i = 0; i < (int)edgeList.size(); i++) {
            auto [x, y, z] = edgeList[i];
            ++degree[x], ++degree[y];
        }
        for (int i = 0; i < N; i++) {
            eid[i].reserve(degree[i]);
        }
        for (int i = 0; i < (int)edgeList.size(); i++) {
            auto [x, y, z] = edgeList[i];
            f += z;
            if (i + 1 < (int)edgeList.size()) {
                auto [nx, ny, nz] = edgeList[i + 1];
                if (nx == x && ny == y) {
                    continue;
                }
            }
            addedge(x, y, f);
            f = 0;
        }
    }
    void clearFlow() {
        for (int i = 0; i < (int)e.size(); i += 2) {
            e[i].f += e[i + 1].f;
            e[i + 1].f = 0;
        }
    }
    void map(int mapN, const std::vector<int>& id) {
        std::vector<std::vector<int>> tmp(mapN);
        for (int i = 0; i < N; i++) {
            if (tmp[id[i]].size() < eid[i].size()) {
                tmp[id[i]].swap(eid[i]);
            }
            for (auto v : eid[i]) {
                tmp[id[i]].emplace_back(v);
            }
        }
        swap(eid, tmp);
        for (auto& [to, f] : e) {
            to = id[to];
        }
        S = id[S], T = id[T], N = mapN;
    }
    std::vector<int> h, gap, cur;
    std::vector<FlowType> flow;
    std::vector<std::vector<int>> pts;
    std::vector<std::list<int>> hp;
    std::vector<std::list<int>::iterator> hpi;
    int top;
    void hlpp_append() {
        h.emplace_back();
        gap.emplace_back();
        cur.emplace_back();
        flow.emplace_back();
        pts.emplace_back();
        hp.emplace_back();
        hpi.emplace_back();
    }
    void bfs() {
        h.assign(N, inf);
        cur.assign(N, 0);
        std::queue<int> que;
        que.push(T), h[T] = 0;
        while (!que.empty()) {
            int v = que.front();
            que.pop();
            for (auto ei : eid[v]) {
                int u = e[ei].to;
                if (!eqls(e[ei ^ 1].f, FlowType(0)) && h[u] == inf) {
                    h[u] = h[v] + 1;
                    que.push(u);
                }
            }
        }
    }
    void hlpp_init() {
        bfs();
        hp.assign(N, std::list<int>{});
        hpi.assign(N, std::list<int>::iterator());
        gap.assign(N, 0);
        h[S] = std::max(h[S], inf - 1);
        for (int i = 0; i < N; i++) {
            if (h[i] < N) {
                ++gap[h[i]];
                hpi[i] = hp[h[i]].insert(hp[h[i]].end(), i);
            }
        }
    }
    void update(int v) {
        if (h[v] < N && v != T) {
            top = std::max(top, h[v]);
            pts[h[v]].push_back(v);
        }
    }
    void pushflow(int v, int u, int ei, FlowType f) {
        e[ei].f -= f, e[ei ^ 1].f += f;
        flow[v] -= f, flow[u] += f;
        if (eqls(flow[u], f)) {
            update(u);
        }
    }
    bool init() {
        hlpp_init();
        flow.assign(N, 0);
        pts.assign(N, std::vector<int>{});
        top = 0;
        if (h[S] == inf) {
            return false;
        }
        for (auto ei : eid[S]) {
            if (!eqls(e[ei].f, FlowType(0))) {
                pushflow(S, e[ei].to, ei, e[ei].f);
            }
        }
        return true;
    }
    void push(int v) {
        for (int& eii = cur[v]; eii < (int)eid[v].size(); eii++) {
            int ei = eid[v][eii];
            int u = e[ei].to;
            if (!eqls(e[ei].f, FlowType(0)) && h[u] == h[v] - 1) {
                pushflow(v, u, ei, std::min(flow[v], (FlowType)e[ei].f));
            }
            if (eqls(flow[v], FlowType(0))) {
                break;
            }
        }
    }
    void relabel(int v) {
        cur[v] = 0;
        int val = inf;
        for (auto ei : eid[v]) {
            if (!eqls(e[ei].f, FlowType(0))) {
                val = std::min(val, h[e[ei].to] + 1);
            }
        }
        if (val >= N) {
            val = inf;
        }
        hp[h[v]].erase(hpi[v]);
        if (!--gap[h[v]]) {
            for (int i = h[v] + 1; i < N; i++) {
                pts[i].clear();
                gap[i] = 0;
                for (auto it = hp[i].begin(); it != hp[i].end();
                     it = hp[i].erase(it)) {
                    h[*it] = inf;
                }
            }
            h[v] = inf;
        } else {
            h[v] = val;
            if (val < N) {
                update(v);
                ++gap[val];
                hpi[v] = hp[h[v]].insert(hp[h[v]].end(), v);
            }
        }
    }
    int select() {
        while (top >= 0 && pts[top].empty()) {
            --top;
        }
        if (top < 0) {
            return -1;
        }
        int res = pts[top].back();
        pts[top].pop_back();
        return res;
    }
    void rebuild() {
        hlpp_init();
        pts.assign(N, std::vector<int>{});
        top = 0;
        for (int i = 0; i < N; i++) {
            if (!eqls(flow[i], FlowType(0))) {
                update(i);
            }
        }
    }
    FlowType hlpp() {
        if (!init()) {
            return 0;
        }
        int cnt = 0, tag = (1ll * N * N / sqrt(e.size()) + 1) * 4;
        for (int v = select(); v != -1; v = select()) {
            push(v);
            if (!eqls(flow[v], FlowType(0))) {
                relabel(v);
            }
            cnt = (cnt + 1) % tag;
            if (cnt == 0) {
                rebuild();
            }
        }
        return flow[T];
    }
    FlowType dynamic_hlpp(const std::function<void(int, int)>& update) {
        if (!init()) {
            return 0;
        }
        int cnt = 0, tag = (1ll * N * N / sqrt(e.size()) + 1) * 4;
        for (int v = select(); v != -1; v = select()) {
            push(v);
            if (!eqls(flow[v], FlowType(0))) {
                update(v, h[v]);
                relabel(v);
            }
            cnt = (cnt + 1) % tag;
            if (cnt == 0) {
                rebuild();
            }
        }
        return flow[T];
    }
    bool dinic_init() {
        bfs();
        cur.assign(N, 0);
        return h[S] < N;
    }
    FlowType dinic_search(int v, FlowType flow) {
        if (v == T) {
            return flow;
        }
        FlowType ref = flow;
        for (int& eii = cur[v];
             eii < (int)eid[v].size() && !eqls(ref, FlowType(0)); eii++) {
            int ei = eid[v][eii], u = e[ei].to;
            if (h[u] == h[v] - 1 && !eqls(e[ei].f, FlowType(0))) {
                FlowType fl = dinic_search(u, std::min(ref, (FlowType)e[ei].f));
                ref -= fl, e[ei].f -= fl, e[ei ^ 1].f += fl;
            }
        }
        return flow - ref;
    }
    FlowType dinic() {
        FlowType flow = 0, total = 0;
        for (auto ei : eid[S]) {
            total = std::min(getinf(FlowType(0)), total + e[ei].f);
        }
        while (dinic_init()) {
            flow += dinic_search(S, total);
        }
        return flow;
    }

    // 执行前需要先执行hlpp或dinic获取最大流
    std::vector<int> getTCut() {
        std::vector<int> res;
        for (int i = 0; i < N; i++) {
            if (h[i] < N) {
                res.emplace_back(i);
            }
        }
        return res;
    }

    // 执行前需要先执行hlpp或dinic获取最大流
    std::vector<int> getSCut() {
        std::vector<int> res;
        for (int i = 0; i < N; i++) {
            if (h[i] >= N) {
                res.emplace_back(i);
            }
        }
        return res;
    }
};
}  // namespace MaxFlow

using u32 = unsigned int;
using i64 = long long;
using u64 = unsigned long long;

const string logFile = "results.csv";

mt19937 rnd;

struct RandSelect {
    int count = 2;

    RandSelect() {}
    RandSelect(int c) : count(c) {}

    void reset() { count = 2; }
    bool isSelect() { return rnd() % (count++) == 0; }
};

struct Bitset {
    vector<u32> b;
    Bitset(int n = 0) : b(((n - 1) >> 5) + 1) {}
    Bitset(vector<bool> p) : b(((p.size() - 1) >> 5) + 1) {
        for (int i = 0; i < p.size(); i++) {
            if (p[i]) {
                set(i);
            }
        }
    }

    void expand(int x) {
        if ((int)b.size() <= (x >> 5)) {
            b.resize((x >> 5) + 1);
        }
    }
    void insert(int x) {
        expand(x);
        set(x);
    }
    void erase(int x) {
        expand(x);
        reset(x);
    }
    bool empty() const {
        for (auto v : b) {
            if (v != 0u) {
                return false;
            }
        }
        return true;
    }

    void set(int x) { b[x >> 5] |= 1 << (x & 31); }
    void reset(int x) { b[x >> 5] &= ~(1 << (x & 31)); }
    void clear() { b.assign(b.size(), 0); }

    friend Bitset operator|(Bitset lhs, const Bitset& rhs) {
        for (int i = 0; i < lhs.b.size(); i++) {
            lhs.b[i] |= rhs.b[i];
        }
        return lhs;
    }
    void operator|=(const Bitset& rhs) {
        for (int i = 0; i < b.size(); i++) {
            b[i] |= rhs.b[i];
        }
    }

    int count() {
        int res = 0;
        for (int i = 0; i < b.size(); i++) {
            res += __builtin_popcount(b[i]);
        }
        return res;
    }

    bool operator<(const Bitset& t) const { return b < t.b; }

    bool operator()(int i) const { return b[i >> 5] >> (i & 31) & 1; }

    int find_first() const {
        for (int i = 0; i < (int)b.size(); i++) {
            if (b[i] != 0u) {
                return __builtin_ctz(b[i]) + (i << 5);
            }
        }
        return -1;
    }
    int find_next(int i) const {
        int x = i >> 5, y = i & 31;
        if (i >= 0 && x < (int)b.size()) {
            unsigned v = (b[x] >> y) & (~1u);
            if (v) {
                return __builtin_ctz(v) + i;
            }
            for (int i = x + 1; i < (int)b.size(); i++) {
                if (b[i] != 0u) {
                    return __builtin_ctz(b[i]) + (i << 5);
                }
            }
        }
        return -1;
    }

    struct iterator {
        const Bitset* const bs;
        int pos;
        iterator(const Bitset* bs, int pos) : bs(bs), pos(pos) {}
        bool operator!=(const iterator& it) const {
            return pair{bs, pos} != pair{it.bs, it.pos};
        }
        int operator*() const { return pos; }
        iterator& operator++() {
            pos = bs->find_next(pos);
            return *this;
        }
    };
    iterator begin() const { return iterator(this, find_first()); }
    iterator end() const { return iterator(this, -1); }
    iterator find(int i) const {
        if (this->operator()(i)) {
            return iterator(this, i);
        }
        return end();
    }
};

using Set = Bitset;

constexpr int inf = 1E9;
constexpr double eps = 0.1;
constexpr int alpha = 5;

using TabuTab = vector<vector<vector<int>>>;

struct Conflict {
    int edge = inf, color = inf;
    Conflict(int edge = inf, int color = inf) : edge(edge), color(color) {}

    bool operator<(const Conflict& t) const {
        return edge == t.edge ? color < t.color : edge < t.edge;
    }

    bool operator==(const Conflict& t) const {
        return edge == t.edge && color == t.color;
    }

    friend Conflict operator+(Conflict a, Conflict b) {
        return Conflict(a.edge + b.edge, a.color + b.color);
    }

    friend ostream& operator<<(ostream& os, const Conflict& t) {
        os << "(" << t.edge << ", " << t.color << ")";
        return os;
    }
    bool operator!() const { return !edge && !color; }
};

const Conflict success = Conflict(0, 0);

struct Reduce {
    Conflict r;
    int row, i, j;
    Reduce(Conflict r = Conflict(), int row = 0, int i = 0, int j = 0)
        : r(r), row(row), i(i), j(j) {}
};

struct Sol {
    static vector<vector<vector<bool>>> fea;
    static vector<vector<vector<int>>> D;
    static int n;
    static vector<vector<int>> fixed, flexiblePos, flexibleVal;

    vector<vector<int>> lsc;

    Conflict conflict;
    vector<vector<int>> conflictC, conflictR;

    vector<Set> bi;
    Set rows;

#define id(i, j) ((i) * n + (j))

    Sol() : lsc(fixed) {
        for (int i = 0; i < n; i++) {
            Shuffle(i);
        }
    }

    Sol(vector<vector<int>> lsc) : lsc(lsc) {}

    void init() {
        conflict = Conflict(0, 0);
        bi.assign(n, {});
        conflictC.assign(n, vector<int>(n));
        conflictR.assign(n, vector<int>(n));

        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                if (conflictC[j][lsc[i][j]]) {
                    if (fixed[i][j] == -1) {
                        bi[i].insert(j);
                    }
                    if (conflictC[j][lsc[i][j]] == 1) {
                        if (fixed[conflictR[j][lsc[i][j]]][j] == -1) {
                            bi[conflictR[j][lsc[i][j]]].insert(j);
                        }
                    }
                }
                if (fixed[i][j] == -1 && fea[i][j][lsc[i][j]] == false) {
                    conflict.color++;
                    bi[i].insert(j);
                }
                conflict.edge += (conflictC[j][lsc[i][j]]++);
                conflictR[j][lsc[i][j]] ^= i;
            }
        }

        for (int i = 0; i < n; i++) {
            if (!bi[i].empty()) {
                rows.insert(i);
            }
        }
    }

    void Set(int i, int j, int c) {
        int src = lsc[i][j];
        conflictC[j][src]--;
        conflictR[j][src] ^= i;

        if (conflictC[j][src] == 1) {
            int row = conflictR[j][src];
            if (fixed[row][j] == -1 && fea[row][j][lsc[row][j]]) {
                bi[row].erase(j);
                if (bi[row].empty()) {
                    rows.erase(row);
                }
            }
        }

        lsc[i][j] = c;

        if (conflictC[j][c] == 1) {
            if (fixed[conflictR[j][c]][j] == -1) {
                bi[conflictR[j][c]].insert(j);
                rows.insert(conflictR[j][c]);
            }
        }
        conflictC[j][c]++;
        conflictR[j][c] ^= i;

        if (conflictC[j][c] > 1 || !fea[i][j][c]) {
            bi[i].insert(j);
            rows.insert(i);
        } else {
            bi[i].erase(j);
            if (bi[i].empty()) {
                rows.erase(i);
            }
        }
    }

    tuple<Reduce, Reduce> getReduce(TabuTab& tabu, int iter) {
        Reduce tb, ntb, rd, nrd;
        RandSelect stb, sntb, srd(1), snrd(1);

        auto Tabu = [&](int row, int i, int j) -> bool {
            return tabu[row][i][lsc[row][j]] > iter ||
                   tabu[row][j][lsc[row][i]] > iter;
        };

        for (auto row : rows) {
            for (auto i : bi[row]) {
                Reduce best, nbest;
                RandSelect sbest, snbest;
                for (auto j : flexiblePos[row]) {
                    if (i == j) {
                        continue;
                    }
                    int re = conflictC[i][lsc[row][i]] +
                             conflictC[j][lsc[row][j]] - 2 -
                             conflictC[i][lsc[row][j]] -
                             conflictC[j][lsc[row][i]];
                    int rc =
                        -fea[row][i][lsc[row][i]] - fea[row][j][lsc[row][j]] +
                        fea[row][i][lsc[row][j]] + fea[row][j][lsc[row][i]];
                    Conflict r = Conflict(-re, -rc);
                    auto& R = Tabu(row, i, j) ? tb : ntb;
                    auto& s = Tabu(row, i, j) ? stb : sntb;

                    auto& bR = Tabu(row, i, j) ? best : nbest;
                    auto& bs = Tabu(row, i, j) ? sbest : snbest;
                    if (r < bR.r || (r == bR.r && bs.isSelect())) {
                        if (r < bR.r) {
                            bs.reset();
                        }
                        bR = Reduce(r, row, i, j);
                    }
                    if (r < R.r || (r == R.r && s.isSelect())) {
                        if (r < R.r) {
                            s.reset();
                        }
                        R = Reduce(r, row, i, j);
                    }
                }
                if (best.r.edge < inf && srd.isSelect()) {
                    rd = best;
                }
                if (nbest.r.edge < inf && snrd.isSelect()) {
                    nrd = nbest;
                }
            }
        }
        if (tb.r.edge > 0) {
            tb = rd;
        }
        if (ntb.r.edge > 0) {
            ntb = nrd;
        }
        return tuple<Reduce, Reduce>{tb, ntb};
    }

    friend Sol operator+(const Sol& father, const Sol& mother) {
        vector<vector<int>> count(n, vector<int>(n));
        vector<vector<int>> lsc(n);
        for (int i = 0; i < n; i++) {
            bool flag = true;
            for (int j = 0; j < n; j++) {
                if (count[j][father.lsc[i][j]]) {
                    flag = false;
                }
            }
            if (flag) {
                lsc[i] = father.lsc[i];
                for (int j = 0; j < n; j++) {
                    count[j][father.lsc[i][j]]++;
                }
            } else {
                lsc[i] = mother.lsc[i];
            }
        }
        return Sol(lsc);
    }

    void Shuffle(int i) {
        shuffle(flexibleVal[i].begin(), flexibleVal[i].end(), rnd);
        int k = 0;
        for (auto j : flexiblePos[i]) {
            lsc[i][j] = flexibleVal[i][k++];
        }
    }

    Sol& operator=(const Sol& t) {
        conflict = t.conflict;
        lsc = t.lsc;
        return *this;
    }

    bool operator<(const Sol& t) const { return conflict < t.conflict; }
    bool operator==(const Sol& t) const { return conflict == t.conflict; }
};

int Sol::n = 0;
vector<vector<int>> Sol::fixed = {};
vector<vector<int>> Sol::flexiblePos = {};
vector<vector<int>> Sol::flexibleVal = {};
vector<vector<vector<bool>>> Sol::fea = {};
vector<vector<vector<int>>> Sol::D = {};

bool Reduction(vector<vector<bool>*> p) {
    bool res = false;
    int n = Sol::n;
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            if ((*p[x])[y]) {
                int s = 2 * n, t = s + 1;
                MaxFlow::FlowGraph<Constants::i64> g(2 * n + 2, s, t);

                for (int i = 0; i < n; i++) {
                    g.addedge(s, i, 1);
                    g.addedge(i + n, t, 1);
                    if (i == x) {
                        continue;
                    }
                    for (int j = 0; j < n; j++) {
                        if (j == y) {
                            continue;
                        }
                        if ((*p[i])[j]) {
                            g.addedge(i, j + n, 1);
                        }
                    }
                }
                if (g.hlpp() != n - 1) {
                    (*p[x])[y] = false;
                    res = true;
                }
            }
        }
    }
    return res;
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(0);

    const double T = stod(argv[1]);
    const u32 seed = stoi(argv[2]);
    const clock_t start = clock();

    auto checkTime = [&]() -> bool {
        return double(clock() - start) / CLOCKS_PER_SEC < T - eps;
    };

    rnd = mt19937(seed);

    int n;
    cin >> n;

    Sol::n = n;
    Sol::fixed.assign(n, vector<int>(n, -1));
    Sol::fea.assign(n, vector<vector<bool>>(n, vector<bool>(n, true)));
    Sol::D.assign(n, vector<vector<int>>(n));
    Sol::flexiblePos.assign(n, vector<int>());
    Sol::flexibleVal.assign(n, vector<int>());

    int u, v, w;
    while (cin >> u >> v >> w) {
        Sol::fixed[u][v] = w;
        for (int i = 0; i < n; i++) {
            if (i != u) {
                Sol::fea[i][v][w] = false;
            }
            if (i != v) {
                Sol::fea[u][i][w] = false;
            }
        }
        Sol::fea[u][v].assign(n, false);
        Sol::fea[u][v][w] = true;
    }

    for (int i = 0; i < n; i++) {
        vector<bool> vis(n);
        for (int j = 0; j < n; j++) {
            if (Sol::fixed[i][j] != -1) {
                vis[Sol::fixed[i][j]] = true;
            } else {
                Sol::flexiblePos[i].push_back(j);
            }
        }
        for (int w = 0; w < n; w++) {
            if (!vis[w]) {
                Sol::flexibleVal[i].push_back(w);
            }
        }
    }

    constexpr int preT = 100000;
    constexpr int genT = 100000;

    auto Tabu = [&](auto& sol, int T) {
        sol.init();
        // cerr << "Tabu start " << sol.conflict << '\n';
        auto best = sol;
        const int P = rnd() % 4 + 1;
        const int base = rnd() % 3;

        TabuTab tabu(n, vector<vector<int>>(n, vector<int>(n)));

        int t = 0;
        for (; t < T && checkTime(); t++) {
            auto [tb, ntb] = sol.getReduce(tabu, t);

            auto maxR = (tb.r < ntb.r && sol.conflict + tb.r < best.conflict)
                            ? tb
                            : ntb;
            if (maxR.r.edge >= inf) {
                continue;
            }

            auto Set = [&](int i, int j, int c) {
                int src = sol.lsc[i][j];
                sol.Set(i, j, c);

                tabu[i][j][src] = t + sol.conflict.edge + sol.conflict.color +
                                  base + rnd() % P;
            };

            int di = sol.lsc[maxR.row][maxR.j], dj = sol.lsc[maxR.row][maxR.i];

            Set(maxR.row, maxR.i, di);
            Set(maxR.row, maxR.j, dj);

            sol.conflict = sol.conflict + maxR.r;

            // auto temp = sol.conflict;
            // sol.init();
            // assert(temp == sol.conflict);

            if (sol < best) {
                best = sol;
            }

            if (!best.conflict) {
                break;
            }

            // if (t % 100000 == 0) {
            //     cerr << tot << ' ' << best.conflict << '\n';
            // }
        }
        sol = best;
        // cerr << "Tabu end " << sol.conflict << '\n';
    };

    Sol ans;
    Tabu(ans, preT);

    constexpr int M = 8;
    const int K = M / 4;
    vector<Sol> groups(M);
    if (ans.conflict.edge) {
        for (auto& sol : groups) {
            Tabu(sol, preT);
            if (sol < ans) {
                ans = sol;
            }
            if (!ans.conflict) {
                break;
            }
        }
    }

    if (ans.conflict.edge) {
        auto start_reduction = clock();
        bool loop = true;
        vector<bool> R(n, true), C(n, true);
        while (loop) {
            loop = false;
            for (int i = 0; i < n; i++) {
                if (!R[i]) {
                    continue;
                }
                R[i] = false;
                vector<vector<bool>> temp = Sol::fea[i];
                vector<vector<bool>*> p;
                for (int j = 0; j < n; j++) {
                    p.push_back(&Sol::fea[i][j]);
                }
                auto now = Reduction(p);
                loop |= now;
                if (now) {
                    for (int j = 0; j < n; j++) {
                        if (temp[j] != Sol::fea[i][j]) {
                            C[j] = true;
                        }
                    }
                }
            }
            for (int j = 0; j < n; j++) {
                if (!C[j]) {
                    continue;
                }
                C[j] = false;
                vector<vector<bool>> temp(n);
                for (int i = 0; i < n; i++) {
                    temp[i] = Sol::fea[i][j];
                }
                vector<vector<bool>*> p;
                for (int i = 0; i < n; i++) {
                    p.push_back(&Sol::fea[i][j]);
                }
                auto now = Reduction(p);
                loop |= now;
                if (now) {
                    for (int i = 0; i < n; i++) {
                        if (temp[i] != Sol::fea[i][j]) {
                            R[i] = true;
                        }
                    }
                }
            }
        }

        int fixed = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (Sol::fixed[i][j] != -1) {
                    Sol::D[i][j] = {Sol::fixed[i][j]};
                    continue;
                }
                int tot = 0;
                for (int w = 0; w < n; w++) {
                    if (Sol::fea[i][j][w]) {
                        Sol::D[i][j].push_back(w);
                        tot++;
                    }
                }
                if (tot == 1) {
                    fixed++;
                    Sol::fixed[i][j] = Sol::D[i][j].back();
                    Sol::flexiblePos[i].erase(find(Sol::flexiblePos[i].begin(),
                                                   Sol::flexiblePos[i].end(),
                                                   j));
                    Sol::flexibleVal[i].erase(find(Sol::flexibleVal[i].begin(),
                                                   Sol::flexibleVal[i].end(),
                                                   Sol::D[i][j].back()));
                }
            }
        }
        // cerr << "Reduction Success... " << fixed << " cell fixed\n";

        // cerr << double(clock() - start_reduction) / CLOCKS_PER_SEC << '\n';
        // cerr << double(clock()) / CLOCKS_PER_SEC << '\n';

        // for (int i = 0; i < n; i++) {
        //     for (int j = 0; j < n; j++) {
        //         cerr << Sol::D[i][j].size() << " \n"[j == n - 1];
        //     }
        // }

        ans.init();
        for (auto& sol : groups) {
            sol = Sol();
            Tabu(sol, preT);
            if (sol < ans) {
                ans = sol;
            }
            if (!ans.conflict) {
                break;
            }
        }
        sort(groups.begin(), groups.end());

        vector<int> best;
        for (int i = 0; i < M; i++) {
            if (groups[i] == groups[0]) {
                best.push_back(i);
            }
        }

        vector<Sol> elites;

        for (; checkTime() && ans.conflict.edge;) {
            int x = best[rnd() % best.size()];
            int y = rnd() % groups.size();
            while (x == y) {
                x = best[rnd() % best.size()];
                y = rnd() % groups.size();
            }

            Sol p = groups[x] + groups[y];
            Tabu(p, genT);

            if (p < ans) {
                ans = p;
            }
            if (!ans.conflict) {
                break;
            }

            auto maxG = min_element(groups.begin(), groups.end());
            if (p < *maxG) {
                *maxG = p;
            } else {
                groups.push_back(p);
            }
            if (groups.size() == 2 * M) {
                sort(groups.begin(), groups.end());
                groups.resize(M);

                if (groups[0] == groups.back()) {
                    elites.push_back(groups[rnd() % M]);
                    for (auto& g : groups) {
                        for (int x = 0; x < n; x++) {
                            if (rnd() % 2) {
                                g.Shuffle(x);
                            }
                        }
                        Tabu(g, genT);
                    }
                } else {
                    gp_hash_table<int, null_type> vis;
                    for (int k = 0; k < K; k++) {
                        int i = rnd() % M;
                        while (vis.find(i) != vis.end()) {
                            i = rnd() % M;
                        }
                        vis.insert(i);

                        for (int x = 0; x < n; x++) {
                            if (rnd() % 2) {
                                groups[i].Shuffle(x);
                            }
                        }
                        Tabu(groups[i], genT);
                    }
                }

                if (elites.size() == M) {
                    groups = elites;
                    best.resize(M);
                    iota(best.begin(), best.end(), 0);
                    elites.clear();
                } else {
                    auto minG = min_element(groups.begin(), groups.end());

                    best.clear();
                    for (int i = 0; i < M; i++) {
                        if (groups[i] == *minG) {
                            best.push_back(i);
                        }
                    }
                }
            }
            if (!ans.conflict) {
                break;
            }

            // cerr << elites.size() << ' ' << ans.conflict << '\n';
        }
    }

    for (auto& row : ans.lsc) {
        for (int j = 0; j < n; j++) {
            cout << row[j] << " \n"[j == n - 1];
        }
    }

    auto check = [&](vector<vector<int>>& lsc) -> bool {
        auto is_permutation = [&](vector<int> p) {
            sort(p.begin(), p.end());
            for (int i = 0; i < n; i++) {
                if (p[i] != i) {
                    return false;
                }
            }
            return true;
        };
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (Sol::fixed[i][j] != -1 && lsc[i][j] != Sol::fixed[i][j]) {
                    return false;
                }
            }
            // if (!is_permutation(lsc[i])) {
            //     return false;
            // }
            // vector<int> p(n);
            // for (int j = 0; j < n; j++) {
            //     p[j] = lsc[j][i];
            // }
            // if (!is_permutation(p)) {
            //     return false;
            // }
        }
        return true;
    };

#if defined(MESSAGE)
    ans.init();
    // cerr << ans.conflict << '\n';
    // assert(check(ans.lsc));

    ofstream csvFile(logFile, ios::app);
    if (!check(ans.lsc)) {
        csvFile << "[LogicError] ";
    }
    csvFile << argv[3] << ", "
            << "HEA-ELITE, " << seed << ", "
            << double(clock() - start) / CLOCKS_PER_SEC << ", "
            << ans.conflict.edge << "\n";
    cerr << double(clock() - start) / CLOCKS_PER_SEC << '\n';
#endif

    return 0;
}