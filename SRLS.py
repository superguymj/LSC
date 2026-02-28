#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import datetime as dt
import math
import re
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import List


def usage() -> None:
    print("Usage: python SRLS.py <seed> <instances_name> <instances_path> <output_path>")


def read_square_from_n2_lines(file_path: Path) -> List[List[int]]:
    """
    读取 SRLS 输出：
    - 格式为 n*n 行
    - 每行一个整数（颜色）
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Result file not found: {file_path}")

    lines = file_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    vals = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        # 每行一个颜色，兼容行里有其他字符时提取整数
        m = re.search(r"-?\d+", line)
        if m:
            vals.append(int(m.group()))

    total = len(vals)
    if total == 0:
        raise ValueError("No values found in output file.")

    n = int(math.isqrt(total))
    if n * n != total:
        raise ValueError(f"Output count is {total}, not a perfect square (n*n).")

    square = [vals[i * n:(i + 1) * n] for i in range(n)]
    return square


def count_conflict_pairs(square: List[List[int]]) -> int:
    """
    冲突格子对数 = 行冲突对 + 列冲突对
    若某值在一行/列中出现 c 次，则贡献 C(c,2)
    """
    n = len(square)
    if any(len(row) != n for row in square):
        raise ValueError("Square is not n x n.")

    def pairs(values: List[int]) -> int:
        cnt = Counter(values)
        return sum(v * (v - 1) // 2 for v in cnt.values() if v > 1)

    ans = 0
    for row in square:
        ans += pairs(row)

    for j in range(n):
        col = [square[i][j] for i in range(n)]
        ans += pairs(col)

    return ans


def append_results_csv(csv_path: Path, instance_name: str, seed: str, duration: float, optima: int):
    header = "Date, Instance, Algorithm, RandSeed, Duration, Optima"
    line = ", ".join([
        dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        instance_name,
        "SRLS",
        str(seed),
        f"{duration:.3f}",
        str(optima),
    ])

    need_header = (not csv_path.exists()) or csv_path.stat().st_size == 0

    with csv_path.open("a", encoding="utf-8", newline="") as f:
        if need_header:
            f.write(header + "\n")
        f.write(line + "\n")


def main():
    if len(sys.argv) != 5:
        usage()
        sys.exit(1)

    seed = sys.argv[1]
    instances_name = sys.argv[2]
    instances_path = Path(sys.argv[3])
    output_path = Path(sys.argv[4])

    script_dir = Path(__file__).resolve().parent
    exe_path = script_dir / "SRLS.exe"
    results_csv = script_dir / "results.csv"

    if not exe_path.exists():
        print(f"ERROR: SRLS.exe not found: {exe_path}")
        sys.exit(2)

    if not instances_path.exists():
        print(f"ERROR: Instance file not found: {instances_path}")
        sys.exit(3)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [str(exe_path), str(seed), str(instances_path), str(output_path)]

    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=str(script_dir), capture_output=True, text=True)
    duration = time.perf_counter() - t0

    if proc.returncode != 0:
        print(f"ERROR: SRLS.exe failed (code={proc.returncode})")
        if proc.stdout:
            print("STDOUT:\n", proc.stdout)
        if proc.stderr:
            print("STDERR:\n", proc.stderr)
        sys.exit(proc.returncode)

    # 读取输出并计算冲突
    square = read_square_from_n2_lines(output_path)
    optima = count_conflict_pairs(square)

    append_results_csv(results_csv, instances_name, seed, duration, optima)

    print(f"Done. Duration={duration:.3f}s, Optima={optima}")
    print(f"results.csv -> {results_csv}")


if __name__ == "__main__":
    main()