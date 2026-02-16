#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <algorithm>
#include <ext/pb_ds/assoc_container.hpp>
#include <fstream>
#include <iostream>
#include <vector>

#include "MaxFlow.h"
#include "reduction.h"
#include "utils.h"
#include "Array.h"
#include "Bitset.h"
#include "MinCostFlow.h"

#define MESSAGE

using namespace std;
using namespace __gnu_pbds;

using u32 = unsigned int;
using i64 = long long;
using u64 = unsigned long long;

constexpr int fix = 1E5;

const string logFile = "results.csv";

mt19937 rnd;

constexpr double eps = 0.1;

using TabuTab = Array3D<int>;

struct Conflict {
    int edge = 1E9, color = 1E9;
    Conflict() {}
    Conflict(int edge, int color) : edge(edge), color(color) {}

    bool operator<(const Conflict &t) const {
        return edge == t.edge ? color < t.color : edge < t.edge;
    }

    bool operator==(const Conflict &t) const { return edge == t.edge && color == t.color; }

    friend Conflict operator+(Conflict a, Conflict b) {
        return Conflict(a.edge + b.edge, a.color + b.color);
    }

    friend ostream &operator<<(ostream &os, const Conflict &t) {
        os << "(" << t.edge << ", " << t.color << ")";
        return os;
    }
    bool operator!() const { return !edge; }
};

struct Reduce {
    Conflict r;
    short row, i, j;
    Reduce(Conflict r = Conflict(), short row = 0, short i = 0, short j = 0)
        : r(r), row(row), i(i), j(j) {}
};

struct Sol {
    static vector<vector<vector<bool>>> fea;
    static int n;
    static vector<vector<short>> fixed, flexiblePos, flexibleVal;

    static vector<vector<vector<int>>> den;

    vector<vector<short>> lsc;

    Conflict conflict;
    vector<vector<short>> conflictC, conflictR;

    vector<Bitset> bi;
    Bitset rows;

    Sol() : lsc(fixed), bi(n, Bitset(n)), rows(n) {
		for (int i = 0; i < n; i++) {
			Shuffle(i);
		}


		// vector<vector<short>> val(n, vector<short>(n));
		// vector<int> Row(n);
		// iota(Row.begin(), Row.end(), 0);
		// shuffle(Row.begin(), Row.end(), rnd);
        
		// vector<int> ord(Row);

		// for (auto row : Row) {
		// 	int s = 2 * n, t = s + 1;
		// 	MinCostFlow<int> mcf(t + 1);
		// 	for (int i = 0; i < n; i++) {
		// 		mcf.addEdge(s, i, 1, 0);
		// 		mcf.addEdge(i + n, t, 1, 0);
		// 		shuffle(ord.begin(), ord.end(), rnd);
		// 		for (auto j : ord) {
		// 			if (fea[row][i][j]) {
		// 				mcf.addEdge(i, j + n, 1, val[i][j]);
		// 			}
		// 		}
		// 	}
		// 	mcf.flow(s, t);
		// 	for (int i = 0; i < n; i++) {
		// 		for (int j : mcf.g[i]) {
		// 			auto [v, cap, cost] = mcf.e[j];
		// 			if (n <= v && v < 2 * n && !cap) {
		// 				lsc[row][i] = v - n;
		// 				val[i][v - n]++;
		// 				break;
		// 			}
		// 		}
		// 	}
		// }
    }

    void init() {
        conflict = Conflict(0, 0);
        bi.assign(n, Bitset(n));
		rows.clear();
        conflictC.assign(n, vector<short>(n));
        conflictR.assign(n, vector<short>(n));

        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                auto x = lsc[i][j];
                if (conflictC[j][x]) {
                    if (fixed[i][j] == -1) {
                        bi[i].set(j);
                    }
                    if (conflictC[j][x] == 1) {
                        if (fixed[conflictR[j][x]][j] == -1) {
                            bi[conflictR[j][x]].set(j);
                        }
                    }
                }
                if (fixed[i][j] == -1 && fea[i][j][x] == false) {
                    bi[i].set(j);
                }
                conflict.edge += (conflictC[j][x]++);
                conflictR[j][x] ^= i;

                conflict.color += den[i][j][x];
            }
        }

        for (int i = 0; i < n; i++) {
            if (bi[i].any()) {
                rows.set(i);
            }
        }
    }

    void Set(short i, short j, short c) {
        auto src = lsc[i][j];
        conflictC[j][src]--;
        conflictR[j][src] ^= i;

        if (conflictC[j][src] == 1) {
            auto row = conflictR[j][src];
            if (fixed[row][j] == -1 && fea[row][j][lsc[row][j]]) {
                bi[row].reset(j);
                if (!bi[row].any()) {
                    rows.reset(row);
                }
            }
        }

        lsc[i][j] = c;

        if (conflictC[j][c] == 1) {
            if (fixed[conflictR[j][c]][j] == -1) {
                bi[conflictR[j][c]].set(j);
                rows.set(conflictR[j][c]);
            }
        }
        conflictC[j][c]++;
        conflictR[j][c] ^= i;

        if (conflictC[j][c] > 1 || !fea[i][j][c]) {
            bi[i].set(j);
            rows.set(i);
        } else {
            bi[i].reset(j);
            if (!bi[i].any()) {
                rows.reset(i);
            }
        }
    }

    tuple<Reduce, Reduce, Reduce> getReduce(const TabuTab &tabu, int iter) {
        Reduce tb, ntb, rd;
        RandSelect stb(1), sntb(1), srd(1);
        for (int row = rows.begin(); row < n; row = rows.next(row)) {
            for (int i = bi[row].begin(); i < n; i = bi[row].next(i)) {
                Reduce best;
                RandSelect sbest(1);
                for (auto j : flexiblePos[row]) {
                    if (i == j) {
                        continue;
                    }
                    auto x = lsc[row][i], y = lsc[row][j];
                    int re = conflictC[i][x] + conflictC[j][y] - 2 - conflictC[i][y] - conflictC[j][x];
                    int rc = den[row][i][x] + den[row][j][y] - den[row][i][y] - den[row][j][x];

                    Conflict r = Conflict(-re, -rc);
					bool TabuFlag = tabu(row, i, y) > iter || tabu(row, j, x) > iter;
                    auto &R = TabuFlag ? tb : ntb;
                    auto &s = TabuFlag ? stb : sntb;
                    if (r < best.r || (r == best.r && sbest.isSelect(rnd))) {
                        if (r < best.r) {
                            sbest.reset();
                        }
                        best = Reduce(r, row, i, j);
                    }
                    if (r < R.r || (r == R.r && s.isSelect(rnd))) {
                        if (r < R.r) {
                            s.reset();
                        }
                        R = Reduce(r, row, i, j);
                    }
                }
                if (srd.isSelect(rnd)) {
                    rd = best;
                }
            }
        }
        return {tb, ntb, rd};
    }

    void Shuffle(int row) {
        // shuffle(flexibleVal[row].begin(), flexibleVal[row].end(), rnd);
        // int k = 0;
        // for (auto j : flexiblePos[row]) {
        //     lsc[row][j] = flexibleVal[row][k++];
        // }

        vector<int> ord(n);
        iota(ord.begin(), ord.end(), 0);
        auto &p = fea[row];
        int s = 2 * n + 1, t = s + 1;
        MaxFlow<int> f(t + 1);
        for (int i = 0; i < n; i++) {
        	shuffle(ord.begin(), ord.end(), rnd);
            f.addEdge(s, i, 1);
            f.addEdge(i + n, t, 1);
            for (auto j : ord) {
                if (p[i][j]) {
                    f.addEdge(i, j + n, 1);
                }
            }
        }
        f.flow(s, t);
        for (int i = 0; i < n; i++) {
            for (auto j : f.g[i]) {
                auto [v, c] = f.e[j];
                if (n <= v && v < 2 * n && !c) {
                    lsc[row][i] = v - n;
                    break;
                }
            }
        }
    }

    Sol &operator=(const Sol &t) {
        conflict = t.conflict;
        lsc = t.lsc;
        return *this;
    }

    bool operator<(const Sol &t) const { return conflict < t.conflict; }
    bool operator==(const Sol &t) const { return conflict == t.conflict; }
};

int Sol::n = 0;
vector<vector<short>> Sol::fixed = {};
vector<vector<short>> Sol::flexiblePos = {};
vector<vector<short>> Sol::flexibleVal = {};
vector<vector<vector<bool>>> Sol::fea = {};
vector<vector<vector<int>>> Sol::den = {};

int main(int argc, char *argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(0);

    const double T = stod(argv[1]);
    const u32 seed = stoi(argv[2]);
    const clock_t start = clock();

    auto checkTime = [&]() -> bool { return double(clock() - start) / CLOCKS_PER_SEC < T - eps; };

    rnd = mt19937(seed);

    int n;
    cin >> n;

    Sol::n = n;
    Sol::fixed.assign(n, vector<short>(n, -1));
    Sol::fea.assign(n, vector<vector<bool>>(n, vector<bool>(n, true)));
    Sol::flexiblePos.assign(n, vector<short>());
    Sol::flexibleVal.assign(n, vector<short>());
    Sol::den.assign(n, vector<vector<int>>(n, vector<int>(n)));

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

	RotateReduction(Sol::fea);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int tot = accumulate(Sol::fea[i][j].begin(), Sol::fea[i][j].end(), 0);
            if (tot == 1) {
                Sol::fixed[i][j] = find(Sol::fea[i][j].begin(), Sol::fea[i][j].end(), true) -
                                   Sol::fea[i][j].begin();
            }
        }
    }
    for (int i = 0; i < n; i++) {
        vector<bool> vis(n);
        for (int j = 0; j < n; j++) {
            if (Sol::fixed[i][j] != -1) {
                vis[Sol::fixed[i][j]] = true;
            } else {
                Sol::flexiblePos[i].push_back(j);
            }
            for (int w = 0; w < n; w++) {
                if (Sol::fea[i][j][w]) {
                    Sol::den[i][j][w] = 0;
                } else {
                    Sol::den[i][j][w] = fix;
                }
            }
        }
        for (int w = 0; w < n; w++) {
            if (!vis[w]) {
                Sol::flexibleVal[i].push_back(w);
            }
        }
    }

    auto den = Sol::den;

    constexpr int genT = 1500;

    set<vector<vector<short>>> st;

    auto Tabu = [&](auto &sol, int T) {
        sol.init();
        // cerr << "Tabu start " << sol.conflict << '\n';
        auto best = sol;
        const int P = rnd() % 4 + 1;
        constexpr int base = 10;

        TabuTab tabu(n);

        int iter = 0, count = 1;
        for (; count < T && checkTime(); count++, iter++) {
            auto [tb, ntb, rd] = sol.getReduce(tabu, iter);

            auto maxR = (tb.r < ntb.r && sol.conflict + tb.r < best.conflict) ? tb : ntb;
            if (maxR.r.edge > 0) {
                maxR = rd;
            }

            auto Set = [&](short i, short j, short c) {
                auto src = sol.lsc[i][j];
                sol.Set(i, j, c);

                tabu(i, j, src) = iter + sol.conflict.edge + base + rnd() % P;
            };

            auto di = sol.lsc[maxR.row][maxR.j], dj = sol.lsc[maxR.row][maxR.i];

            Set(maxR.row, maxR.i, di);
            Set(maxR.row, maxR.j, dj);

            sol.conflict = sol.conflict + maxR.r;

            if (sol < best) {
				if (count > 100) {
					// cerr << "Tabu: " << count << '\n';
				}
                count = 1;
                best = sol;
            }

            if (!best.conflict) {
                break;
            }
        }
        sol = best;
        // cerr << "Tabu end " << sol.conflict << '\n';

        {
            if (sol.conflict.edge < 6 && !st.count(sol.lsc)) {
                st.insert(sol.lsc);
                // cerr << sol.conflict << '\n';
                int c = 1;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (sol.conflictC[j][sol.lsc[i][j]] == 1) {
                            Sol::den[i][j][sol.lsc[i][j]] -= c;
                        } else {
                            Sol::den[i][j][sol.lsc[i][j]] += c;
                        }
                    }
                }
            }
        }
    };

    Sol ans;
    ans.init();
    int iter = 0;
    for (; checkTime() && ans.conflict.edge; iter++) {
        Sol sol;
        Tabu(sol, genT);
        if (sol.conflict < ans.conflict) {
            ans = sol;
        }
        if (iter == 400) {
            iter = 0;
            Sol::den = den;
            ans.init();
        }
    }
    cerr << "Iter = " << iter << '\n';

    for (auto &row : ans.lsc) {
        for (int j = 0; j < n; j++) {
            cout << row[j] << " \n"[j == n - 1];
        }
    }

    auto check = [&](vector<vector<short>> &lsc) -> bool {
        auto is_permutation = [&](vector<short> p) {
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
            if (!is_permutation(lsc[i])) {
                return false;
            }
            vector<short> p(n);
            for (int j = 0; j < n; j++) {
                p[j] = lsc[j][i];
            }
            if (!is_permutation(p)) {
                return false;
            }
        }
        return true;
    };

    // {
    //     ofstream fout("data.txt");
    //     for (int i = 0; i < n; i++) {
    //         for (int j = 0; j < n; j++) {
    //             for (int k = 0; k < n; k++) {
    //                 fout << Sol::den[i][j][k] << " \n"[k == n - 1];
    //             }
    //         }
    //     }
    // }

#if defined(MESSAGE)
    ans.init();
    cerr << ans.conflict << '\n';
    // assert(check(ans.lsc));

    ofstream csvFile(logFile, ios::app);
    // if (!check(ans.lsc)) {
    //     csvFile << "[LogicError] ";
    // }
    csvFile << now_str() << ", " << argv[3] << ", "
            << "TABU(flow-400), " << seed << ", " << double(clock() - start) / CLOCKS_PER_SEC
            << ", " << ans.conflict.edge << "\n";
    cerr << double(clock() - start) / CLOCKS_PER_SEC << '\n';
#endif

    return 0;
}