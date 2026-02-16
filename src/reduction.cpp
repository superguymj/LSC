#include "reduction.h"
#include "Bitset.h"
#include "MaxFlow.h"
#include "tarjan.h"

#include <algorithm>
#include <map>
#include <iostream>

bool dfs(std::vector<std::vector<bool> *> p) {
    if (p.size() <= 1) {
        return false;
    }
	const int n = p[0]->size();
    bool res = false;
    for (int i = 0; i < (int)p.size(); i++) {
        Bitset v(*p[i]);
        if (v.count() == 1) {
            p.erase(p.begin() + i);
            for (int w = 0; w < n; w++) {
                if (v(w)) {
                    for (auto it : p) {
                        if ((*it)[w]) {
                            (*it)[w] = false;
                            res = true;
                        }
                    }
                }
            }
            return res || dfs(p);
        }
    }

    std::vector<std::map<Bitset, std::vector<int>>> f(p.size() + 1);
    for (int i = 0; i < (int)p.size(); i++) {
        Bitset v(*p[i]);
        for (int sz = std::min(i, (int)p.size() - 2); sz; sz--) {
            for (auto &[b, _] : f[sz]) {
                auto u = b | v;
                if (u.count() == sz + 1) {
                    std::vector<std::vector<bool> *> rp{p[i]};
                    p.erase(p.begin() + i);
                    std::reverse(_.begin(), _.end());
                    for (auto x : _) {
                        rp.push_back(p[x]);
                        p.erase(p.begin() + x);
                    }
                    for (int w = 0; w < n; w++) {
                        if (u(w)) {
                            for (auto it : p) {
                                if ((*it)[w]) {
                                    (*it)[w] = false;
                                    res = true;
                                }
                            }
                        }
                    }
                    return res || dfs(p) || dfs(rp);
                }
                if (!f[sz + 1].count(u)) {
                    auto &r = f[sz + 1][u];
                    r = _, r.push_back(i);
                }
            }
        }
        f[1][v] = {i};
    }
    return res;
}

bool flow(std::vector<std::vector<bool> *> p) {
	const int n = p.size();
	int s = 2 * n + 1, t = s + 1;
	MaxFlow<int> f(t + 1);
	for (int i = 0; i < n; i++) {
		f.addEdge(s, i, 1);
		f.addEdge(i + n, t, 1);
		for (int j = 0; j < n; j++) {
			if ((*p[i])[j]) {
				f.addEdge(i, j + n, 1);
			}
		}
	}
	f.flow(s, t);

	std::vector<int> to(n);
	for (int i = 0; i < n; i++) {
		for (auto j : f.g[i]) {
			auto [v, c] = f.e[j];
			if (n <= v && v < 2 * n && !c) {
				to[v - n] = i;
				break;
			}
		}
	}

	SCC scc(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if ((*p[i])[j]) {
				scc.addEdge(i, to[j]);
			}
		}
	}
	scc.work();
	
	bool flag = false;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if ((*p[i])[j] && scc.bel[i] != scc.bel[to[j]]) {
				flag = true;
				(*p[i])[j] = false;
			}
		}
	}

	return flag;
}

bool Reduction(std::vector<std::vector<std::vector<bool>>> &fea) {
	const int n = fea.size();
    bool loop = true, modify = false;
    std::vector<bool> R(n, true), C(n, true);
    while (loop) {
        loop = false;
        for (int i = 0; i < n; i++) {
            if (!R[i]) {
                continue;
            }
            R[i] = false;
            std::vector<std::vector<bool>> temp = fea[i];
            std::vector<std::vector<bool> *> p;
            for (int j = 0; j < n; j++) {
                p.push_back(&fea[i][j]);
            }
            auto now = flow(p);
            loop |= now;
			modify |= now;
            if (now) {
                for (int j = 0; j < n; j++) {
                    if (temp[j] != fea[i][j]) {
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
            std::vector<std::vector<bool>> temp(n);
            for (int i = 0; i < n; i++) {
                temp[i] = fea[i][j];
            }
            std::vector<std::vector<bool> *> p;
            for (int i = 0; i < n; i++) {
                p.push_back(&fea[i][j]);
            }
            auto now = flow(p);
            loop |= now;
			modify |= now;
            if (now) {
                for (int i = 0; i < n; i++) {
                    if (temp[i] != fea[i][j]) {
                        R[i] = true;
                    }
                }
            }
        }
    }
	// std::cerr << modify << '\n';
	return modify;
}

std::vector<std::vector<std::vector<bool>>> Rotate(const std::vector<std::vector<std::vector<bool>>> fea) {
	const int n = fea.size();
	std::vector<std::vector<std::vector<bool>>> rfea(fea);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				rfea[j][k][i] = fea[i][j][k];
			}
		}
	}
	return rfea;
}

void RotateReduction(std::vector<std::vector<std::vector<bool>>>& fea) {
	int iter = 0;
	while (Reduction(fea) || iter) {
		iter = (iter + 1) % 3;
		Rotate(fea).swap(fea);
	}
}