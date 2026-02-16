# 编译器
CXX = g++

# 编译选项
CXXFLAGS = -std=c++23 -O3 -march=native -funroll-loops -ftree-vectorize

# 源文件目录
SRC_DIR = src

# 目标文件
TARGET = lsc

# 源文件
SRC = $(SRC_DIR)/Tabu.cpp $(SRC_DIR)/utils.cpp $(SRC_DIR)/reduction.cpp

# 规则部分
all: $(TARGET)

# 生成可执行文件
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

# 清理目标文件
clean:
	rm -f $(TARGET)

debug: $(SRC)
	$(CXX) -o lsc $^ -g -Wall -Wextra

# PHONY 目标，防止和文件名冲突
.PHONY: all clean