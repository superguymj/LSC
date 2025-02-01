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

bool Reduction(vector<vector<bool>*> p) {
    if (p.size() <= 1) {
        return false;
    }
    bool res = false;
    for (int i = 0; i < p.size(); i++) {
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

int main() {
    cin >> n;

    vector<vector<vector<bool>>> vis(n, vector<vector<bool>>(n, vector<bool>(n, true)));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int x;
            cin >> x;
            if (x != -1) {
                vis[i][j].assign(n, false);
                vis[i][j][x] = true;
            }
        }
    }

    bool loop = true;
    vector<bool> R(n, true), C(n, true);
    while (loop) {
        loop = false;
        for (int i = 0; i < n; i++) {
            if (!R[i]) {
                continue;
            }
            R[i] = false;
            vector<vector<bool>> temp = vis[i];
            vector<vector<bool>*> p;
            for (int j = 0; j < n; j++) {
                p.push_back(&vis[i][j]);
            }
            auto now = Reduction(p);
            loop |= now;
            if (now) {
                for (int j = 0; j < n; j++) {
                    if (temp[j] != vis[i][j]) {
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
                temp[i] = vis[i][j];
            }
            vector<vector<bool>*> p;
            for (int i = 0; i < n; i++) {
                p.push_back(&vis[i][j]);
            }
            auto now = Reduction(p);
            loop |= now;
            if (now) {
                for (int i = 0; i < n; i++) {
                    if (temp[i] != vis[i][j]) {
                        R[i] = true;
                    }
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int w = -1;
            int tot = 0;
            for (int k = 0; k < n; k++) {
                if (vis[i][j][k]) {
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