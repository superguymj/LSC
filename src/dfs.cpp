#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

#include <algorithm>
#include <ext/pb_ds/assoc_container.hpp>
#include <fstream>
#include <iostream>
#include <vector>

#include "Array.h"
#include "Bitset.h"
#include "MaxFlow.h"
#include "MinCostFlow.h"
#include "reduction.h"
#include "utils.h"

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
    int n;
    vector<vector<vector<bool>>> fea;
    vector<vector<short>> fixed, flexiblePos, flexibleVal;

    vector<vector<vector<int>>> den;

    vector<vector<short>> lsc;

    Conflict conflict;
    vector<vector<short>> conflictC, conflictR;

    vector<Bitset> bi;
    Bitset rows;

    Sol(vector<vector<vector<bool>>> _fea) : n(_fea.size()), fea(_fea), rows(n) {
        fixed.assign(n, vector<short>(n, -1));
        flexiblePos.assign(n, {});
        flexibleVal.assign(n, {});
        den.assign(n, vector<vector<int>>(n, vector<int>(n)));

        for (int i = 0; i < n; i++) {
            vector<bool> vis(n);
            for (int j = 0; j < n; j++) {
                int tot = accumulate(fea[i][j].begin(), fea[i][j].end(), 0);
                if (tot == 1) {
                    int w = find(fea[i][j].begin(), fea[i][j].end(), true) - fea[i][j].begin();
                    fixed[i][j] = w;
                    vis[w] = true;
                } else {
                    flexiblePos[i].push_back(j);
                }
                for (int w = 0; w < n; w++) {
                    if (!fea[i][j][w]) {
                        den[i][j][w] = fix;
                    }
                }
            }
            for (int w = 0; w < n; w++) {
                if (!vis[w]) {
                    flexibleVal[i].push_back(w);
                }
            }
        }
        lsc = fixed;

        for (int i = 0; i < n; i++) {
            Shuffle(i);
        }

        init();
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
                    int re =
                        conflictC[i][x] + conflictC[j][y] - 2 - conflictC[i][y] - conflictC[j][x];
                    int rc = den[row][i][x] + den[row][j][y] - den[row][i][y] - den[row][j][x];

                    Conflict r = Conflict(-re, -rc);
                    bool TabuFlag = tabu(row, i, y) > iter && tabu(row, j, x) > iter;
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

    vector<vector<vector<bool>>> fea(n, vector<vector<bool>>(n, vector<bool>(n, true)));
	vector<vector<short>> fixed(n, vector<short>(n, -1));
    int u, v, w;
    while (cin >> u >> v >> w) {
		fixed[u][v] = w;
        for (int i = 0; i < n; i++) {
            if (i != u) {
                fea[i][v][w] = false;
            }
            if (i != v) {
                fea[u][i][w] = false;
            }
        }
        fea[u][v].assign(n, false);
        fea[u][v][w] = true;
    }

    RotateReduction(fea);

    constexpr int genT = 1000;

    set<vector<vector<short>>> st;

    auto Tabu = [&](auto &sol, int T) {
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
            if (sol.conflict.edge < 10) {
                // cerr << sol.conflict << '\n';
                int c = 1;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (sol.conflictC[j][sol.lsc[i][j]] == 1) {
                            sol.den[i][j][sol.lsc[i][j]] -= c;
                        } else {
                            sol.den[i][j][sol.lsc[i][j]] += c;
                        }
                    }
                }
            }
        }
    };

    Sol ans(fea);
    queue<vector<vector<vector<bool>>>> q;
    q.push(fea);

    while (q.size() && ans.conflict.edge) {
        auto fea = q.front();
        q.pop();

		Sol sol(fea);
        for (int iter = 0; checkTime() && iter < 1; iter++) {
            Tabu(sol, genT);

			// cerr << "sol = " << sol.conflict << '\n';
            if (sol < ans) {
                ans = sol;
                break;
            }

			for (int i = 0; i < n; i++) {
				sol.Shuffle(i);
			}
			sol.init();
        }

		cerr << "ans = " << ans.conflict << '\n';
        if (!ans.conflict) {
            break;
        }

		RandSelect s(1);
		int u = -1, v = -1, w = n + 1;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int tot = accumulate(fea[i][j].begin(), fea[i][j].end(), 0);
				if (tot > 1) {
					if (tot < w || (tot == w && s.isSelect(rnd))) {
						if (tot < w) {
							s.reset();
						}
						u = i, v = j, w = tot;
					}
				}
			}
		}


		auto temp = fea;
		vector<int> p;
		for (int w = 0; w < n; w++) {
			if (fea[u][v][w]) {
				p.push_back(w);
				fea[u][v][w] = temp[u][v][w] = false;
			}
		}
		shuffle(p.begin(), p.end(), rnd);

		for (auto w : p) {
			fea[u][v][w] = true;
			if (!RotateReduction(fea)) {
				q.push(fea);
			}
			fea = temp;
		}
		cerr << u << ' ' << v << ' ' << w << '\n';
    }

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
                if (fixed[i][j] != -1 && lsc[i][j] != fixed[i][j]) {
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

#if defined(MESSAGE)
    ans.init();
    cerr << ans.conflict << '\n';
    // assert(check(ans.lsc));

    ofstream csvFile(logFile, ios::app);
    // if (!check(ans.lsc)) {
    //     csvFile << "[LogicError] ";
    // }
    csvFile << now_str() << ", " << argv[3] << ", "
            << "dfs(flow-400), " << seed << ", " << double(clock() - start) / CLOCKS_PER_SEC << ", "
            << ans.conflict.edge << "\n";
    cerr << double(clock() - start) / CLOCKS_PER_SEC << '\n';
#endif

    return 0;
}