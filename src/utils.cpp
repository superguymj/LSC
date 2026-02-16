#include "utils.h"
#include <random>
#include <string>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>

std::string now_str() {
    using namespace std::chrono;
    auto now = system_clock::now();
    std::time_t t = system_clock::to_time_t(now);

    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm, "%F %T"); // %F=YYYY-MM-DD, %T=HH:MM:SS
    return oss.str();
}

RandSelect::RandSelect() = default;
RandSelect::RandSelect(int c) : count(c) {}

void RandSelect::reset() { count = 2; }

bool RandSelect::isSelect(std::mt19937 &rnd) { 
    return rnd() % (count++) == 0; 
}