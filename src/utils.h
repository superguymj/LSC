#include <random>

std::string now_str();

struct RandSelect {
    int count = 2;

    RandSelect();
    explicit RandSelect(int c);

    void reset();
    bool isSelect(std::mt19937 &rnd);
};