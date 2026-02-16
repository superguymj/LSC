#include "reduction.h"
#include "Bitset.h"
#include "MaxFlow.h"
#include "tarjan.h"

#include <algorithm>
#include <map>
#include <iostream>
#include <numeric>

int flow(std::vector<std::vector<bool> *> p) {
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
	if (f.flow(s, t) < n) {
		return -1;
	}

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

int Reduction(std::vector<std::vector<std::vector<bool>>> &fea) {
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
			if (now == -1) {
				return -1;
			}
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
			if (now == -1) {
				return -1;
			}
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

int RotateReduction(std::vector<std::vector<std::vector<bool>>>& fea) {
	int iter = 0, temp;
	while ((temp = Reduction(fea)) || iter) {
		if (temp == -1) {
			return -1;
		}
		iter = (iter + 1) % 3;
		Rotate(fea).swap(fea);
	}
	return 0;
}