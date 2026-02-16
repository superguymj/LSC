#include <algorithm>
#include <vector>

struct Bitset {
    using u64 = unsigned long long;

    static constexpr int W = 64;
    static constexpr int LG = 6;
    static constexpr int MASK = 63;

    std::vector<u64> b;

    Bitset(int n) : b(((n - 1) >> LG) + 1) {}

    Bitset(std::vector<bool> p) : b(((static_cast<int>(p.size()) - 1) >> LG) + 1) {
        for (int i = 0; i < (int)p.size(); i++) {
            if (p[i]) {
                set(i);
            }
        }
    }

    void set(int x) { b[x >> LG] |= (1ULL << (x & MASK)); }
    void reset(int x) { b[x >> LG] &= ~(1ULL << (x & MASK)); }
    void clear() { std::fill(b.begin(), b.end(), 0ULL); }

    friend Bitset operator|(Bitset lhs, const Bitset &rhs) {
        for (int i = 0; i < (int)lhs.b.size(); i++) {
            lhs.b[i] |= rhs.b[i];
        }
        return lhs;
    }

    void operator|=(const Bitset &rhs) {
        for (int i = 0; i < (int)b.size(); i++) {
            b[i] |= rhs.b[i];
        }
    }

    int count() const {
        int res = 0;
        for (int i = 0; i < (int)b.size(); i++) {
            res += __builtin_popcountll(b[i]);
        }
        return res;
    }

    bool any() const {
        for (int i = 0; i < (int)b.size(); i++) {
            if (b[i]) {
                return true;
            }
        }
        return false;
    }

    bool operator<(const Bitset &t) const { return b < t.b; }

    bool operator()(int i) const { return (b[i >> LG] >> (i & MASK)) & 1ULL; }

    int begin() const {
        for (int i = 0; i < (int)b.size(); i++) {
            if (b[i]) {
                return i * W + __builtin_ctzll(b[i]);
            }
        }
        return (int)b.size() * W;
    }

    int next(int i) const {
        int x = i >> LG, y = i & MASK;

        if (y + 1 < W) {
            u64 shifted = b[x] >> (y + 1);
            if (shifted) {
                return i + 1 + __builtin_ctzll(shifted);
            }
        }

        for (x++; x < (int)b.size(); x++) {
            if (b[x]) {
                return x * W + __builtin_ctzll(b[x]);
            }
        }
        return (int)b.size() * W;
    }
};