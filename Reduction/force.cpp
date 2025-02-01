#include <bits/stdc++.h>

using namespace std;
int n;

using u32 = unsigned int;

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

    friend Bitset operator|(Bitset lhs, Bitset rhs) {
        for (int i = 0; i < lhs.b.size(); i++) {
            lhs.b[i] |= rhs.b[i];
        }
        return lhs;
    }
    void operator|=(Bitset rhs) {
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

int main() {
    cin >> n;

    vector<vector<Bitset>> a(n, vector<Bitset>(n, Bitset(n)));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int x;
            cin >> x;
            if (x != -1) {
                a[i][j].set(x);
            } else {
                for (int k = 0; k < n; k++) {
                    a[i][j].set(k);
                }
            }
        }
    }

    vector<Bitset*> p(n);
    auto dfs = [&](auto&& self, int S) -> void {
        for (int i = (S - 1) & S; i; i = (i - 1) & S) {
            Bitset u(n);
            for (int j = 0; j < n; j++) {
                if (i >> j & 1) {
                    u |= *(p[j]);
                }
            }
            if (u.count() == __builtin_popcount(i)) {
                for (int j = 0; j < n; j++) {
                    if (~i >> j & 1) {
                        for (int k = 0; k < n; k++) {
                            if (u(k)) {
                                p[j]->reset(k);
                            }
                        }
                    }
                }
                self(self, i);
                self(self, S ^ i);
                return;
            }
        }
    };

    for (int t = 0; t < n * n; t++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                p[j] = &a[i][j];
            }
            dfs(dfs, (1 << n) - 1);
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                p[j] = &a[j][i];
            }
            dfs(dfs, (1 << n) - 1);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int w = -1;
            int tot = 0;
            for (int k = 0; k < n; k++) {
                if (a[i][j](k)) {
                    tot++;
                    w = k;
                }
            }
            if (tot > 1) {
                w = -1;
            }
            cout << w << " \n"[j == n - 1];
        }
    }

    return 0;
}