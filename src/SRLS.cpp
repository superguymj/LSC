#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#define rep(i,x,y) for (int i = x; i <= y; i++) 
#define repd(i,x,y) for (int i = x; i >= y; i--)

#define MESSAGE

using namespace std;
using namespace __gnu_pbds;

using u32 = unsigned int;
using i64 = long long;
using u64 = unsigned long long;

const string logFile = "results.csv";

mt19937 rnd;

struct RandSelect {
    int count = 2;

    RandSelect() {}
    RandSelect(int c) : count(c) {}

    void reset() {
        count = 2;
    }
    bool isSelect() {
        return rnd() % (count++) == 0;
    }
};

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
    bool operator!() const {
        return !edge && !color;
    }
};

const Conflict success = Conflict(0, 0);

struct Reduce {
    Conflict r;
    int row, i, j;
    Reduce(Conflict r = Conflict(), int row = 0, int i = 0, int j = 0) : r(r), row(row), i(i), j(j) {}
};

struct Sol {
    static vector<vector<vector<bool>>> fea;
    static vector<vector<vector<int>>> D;
    static int n;
    static vector<vector<int>> fixed, flexiblePos, flexibleVal;


    vector<vector<int>> lsc;

    Conflict conflict;
    vector<vector<int>> conflictC, conflictR;

    vector<gp_hash_table<int, null_type>> bi;
    gp_hash_table<int, null_type> rows;

    #define id(i, j) ((i) * n + (j))

    Sol() : lsc(fixed) {
        for (int i = 0; i < n; i++) {
            Shuffle(i);
        }
    }

    Sol(vector<vector<int>> lsc) : lsc(lsc) {
    }

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
            if (bi[i].size()) {
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
                if (!bi[row].size()) {
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
            if (!bi[i].size()) {
                rows.erase(i);
            }
        }

    }

    tuple<Reduce, Reduce, Reduce> getReduce(TabuTab& tabu, int iter) {
        Reduce tb, ntb, rd;
        RandSelect stb, sntb, srd(1);
        for (auto row : rows) {
            for (auto i : bi[row]) {
                Reduce best;
                RandSelect sbest;
                for (auto j : flexiblePos[row]) {
                    if (i == j) {
                        continue;
                    }
                    int re = conflictC[i][lsc[row][i]] + conflictC[j][lsc[row][j]] - 2 - conflictC[i][lsc[row][j]] - conflictC[j][lsc[row][i]];
                    int rc = - fea[row][i][lsc[row][i]] - fea[row][j][lsc[row][j]] + fea[row][i][lsc[row][j]] + fea[row][j][lsc[row][i]];
                    Conflict r = Conflict(-re, -rc);
                    auto& R = (tabu[row][i][lsc[row][j]] > iter || tabu[row][j][lsc[row][i]] > iter) ? tb : ntb;
                    auto& s = (tabu[row][i][lsc[row][j]] > iter || tabu[row][j][lsc[row][i]] > iter) ? stb : sntb;
                    if (r < best.r || (r == best.r && sbest.isSelect())) {
                        if (r < best.r) {
                            sbest.reset();
                        }
                        best = Reduce(r, row, i, j);
                    }
                    if (r < R.r || (r == R.r && s.isSelect())) {
                        if (r < R.r) {
                            s.reset();
                        }
                        R = Reduce(r, row, i, j);
                    }
                }
                if (srd.isSelect()) {
                    rd = best;
                }
            }
        }
        return tuple<Reduce, Reduce, Reduce>{tb, ntb, rd};
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

    bool operator<(const Sol& t) const {
        return conflict < t.conflict;
    }
    bool operator==(const Sol& t) const {
        return conflict == t.conflict;
    }
};

int Sol::n = 0;
vector<vector<int>> Sol::fixed = {};
vector<vector<int>> Sol::flexiblePos = {};
vector<vector<int>> Sol::flexibleVal = {};
vector<vector<vector<bool>>> Sol::fea = {};
vector<vector<vector<int>>> Sol::D = {};

struct Bitset {
    vector<u32> b;
    Bitset(int n) : b(((n - 1) >> 5) + 1) {}
    Bitset(vector<bool> p) : b(((p.size() - 1) >> 5) + 1) {
        for (int i = 0; i < p.size(); i++) {
            if (p[i]) {
                set(i);
            }
        }
    }

    void set(int x) {
        b[x >> 5] |= 1 << (x & 31);
    }
    void reset(int x) {
        b[x >> 5] &= ~(1 << (x & 31));
    }
    void clear() {
        b.assign(b.size(), 0);
    }

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

    bool operator<(const Bitset& t) const {
        return b < t.b;
    }

    bool operator()(int i) {
        return b[i >> 5] >> (i & 31) & 1;
    }
};

bool Reduction(vector<vector<bool>*> p) {
    if (p.size() <= 1) {
        return false;
    }
    bool res = false;
    for (int i = 0; i < p.size(); i++) {
        Bitset v(*p[i]);
        if (v.count() == 1) {
            p.erase(p.begin() + i);
            for (int w = 0; w < Sol::n; w++) {
                if (v(w)) {
                    for (auto it : p) {
                        if ((*it)[w]) {
                            (*it)[w] = false;
                            res = true;
                        }
                    }
                }
            }
            return res | Reduction(p);
        }
    }
    
    vector<map<Bitset, vector<int>>> f(p.size() + 1);
    for (int i = 0; i < p.size(); i++) {
        Bitset v(*p[i]);
        for (int sz = min(i, (int)p.size() - 2); sz; sz--) {
            for (auto& [b, _] : f[sz]) {
                auto u = b | v;
                if (u.count() == sz + 1) {
                    vector<vector<bool>*> rp{p[i]};
                    p.erase(p.begin() + i);
                    reverse(_.begin(), _.end());
                    for (auto x : _) {
                        rp.push_back(p[x]);
                        p.erase(p.begin() + x);
                    }
                    for (int w = 0; w < Sol::n; w++) {
                        if (u(w)) {
                            for (auto it : p) {
                                if ((*it)[w]) {
                                    (*it)[w] = false;
                                    res = true;
                                }
                            }
                        }
                    }
                    // cerr << "Part\n";
                    return res | Reduction(p) | Reduction(rp);
                }
                if (!f[sz + 1].count(u)) {
                    auto& r = f[sz + 1][u];
                    r = _, r.push_back(i);
                }
            }
        }
        f[1][v] = {i};
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
                Sol::flexiblePos[i].erase(find(Sol::flexiblePos[i].begin(), Sol::flexiblePos[i].end(), j));
                Sol::flexibleVal[i].erase(find(Sol::flexibleVal[i].begin(), Sol::flexibleVal[i].end(), Sol::D[i][j].back()));
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

    auto Tabu = [&](auto& sol) {
        sol.init();
        // cerr << "Tabu start " << sol.conflict << '\n';
        auto best = sol;
        const int P = 10;
        const int base = 1;

        constexpr int rt0 = 10;
        constexpr int rt_ub = 15;
        constexpr int accu_ub = 1000;

        TabuTab tabu(n, vector<vector<int>>(n, vector<int>(n)));
        int t = 0, accu = 0, rt = rt0;

        auto AdaptiveRestart = [&]() -> bool {
            if (sol.conflict.edge - best.conflict.edge > rt) {
                if (rt < rt_ub) {
                    accu++;
                    if (accu == accu_ub) {
                        accu = 0;
                        rt++;
                    }
                }
                return true;
            }
            return false;
        };

        for (; checkTime(); t++) {
            auto [tb, ntb, rd] = sol.getReduce(tabu, t);

            auto maxR = (tb.r < ntb.r && sol.conflict + tb.r < best.conflict) ? tb : ntb;
            if (maxR.r.edge > 0) {
                maxR = rd;       
            }

            auto Set = [&](int i, int j, int c) {
                int src = sol.lsc[i][j];
                sol.Set(i, j, c);

                tabu[i][j][src] = t + sol.conflict.edge * 2 / 5 + base + rnd() % P;
            };

            int di = sol.lsc[maxR.row][maxR.j], dj = sol.lsc[maxR.row][maxR.i];

            Set(maxR.row, maxR.i, di);
            Set(maxR.row, maxR.j, dj);
            
            sol.conflict = sol.conflict + maxR.r;

            // auto temp = sol.conflict;
            // sol.init();
            // assert(temp == sol.conflict);

            if (sol < best || sol == best) {
                best = sol;
            } else {
                if (AdaptiveRestart()) {
                    t = 0;
                    tabu.assign(n, vector<vector<int>>(n, vector<int>(n)));
                    sol = best;
                    sol.init();
                }
            }

            if (!best.conflict) {
                break;
            } 
            
            if (t % 100000 == 0) {
                cerr << best.conflict << '\n';
            }
        }
        sol = best;
        // cerr << "Tabu end " << sol.conflict << '\n';
    };

    Sol ans;
    Tabu(ans);

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
            << "SRLS, "
            << seed << ", "
            << double(clock() - start) / CLOCKS_PER_SEC << ", "
            << ans.conflict.edge << "\n";
    cerr << double(clock() - start) / CLOCKS_PER_SEC << '\n';
#endif

    return 0;
}