# 编译器
CXX = g++

# 编译选项
<<<<<<< HEAD
CXXFLAGS = -std=c++23 -O3 -march=native -funroll-loops
=======
CXXFLAGS = -std=c++23 -O3 -funroll-loops -ftree-vectorize
>>>>>>> d4e97354f8c35fd757f1f7dc28825f4963012b02

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