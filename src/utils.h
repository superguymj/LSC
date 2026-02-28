#include <random>
#include <iostream>

std::string now_str();

struct RandSelect {
    int count = 2;

    RandSelect();
    explicit RandSelect(int c);

    void reset();
    bool isSelect(uint32_t (*rng)(uint32_t mod));
};

static uint64_t rng_s;
static inline uint32_t rng32(uint32_t mod) {
    rng_s = rng_s * 6364136223846793005ULL + 1442695040888963407ULL;
	if (mod <= 1) {
		return 0;
	}
	return (rng_s >> 32) % mod;
}
static void rng_seed(uint32_t seed) {
    uint64_t z = (uint64_t)seed;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    rng_s = z ^ (z >> 31);
    for (int i = 0; i < 16; i++)
        rng32(0);
}

template<class InputIt, class = std::_RequireInputIter<InputIt>>
void shuffle(InputIt first, InputIt last, uint32_t (*rng)(uint32_t mod)) {
	int n = last - first;
	for (int i = 0; i < n; i++) {
		int x = rng(i + 1);
		std::iter_swap(first + i, first + x);
	}
}