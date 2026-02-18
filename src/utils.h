#include <random>

std::string now_str();

struct RandSelect {
    int count = 2;

    RandSelect();
    explicit RandSelect(int c);

    void reset();
    bool isSelect(uint32_t (*rng)());
};

static uint64_t rng_s;
static inline uint32_t rng32() {
    rng_s = rng_s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (uint32_t)(rng_s >> 32);
}
static void rng_seed(uint32_t seed) {
    uint64_t z = (uint64_t)seed;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    rng_s = z ^ (z >> 31);
    for (int i = 0; i < 16; i++)
        rng32();
}