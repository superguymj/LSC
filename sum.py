#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd


def load_csv(path: str) -> pd.DataFrame:
    """
    兼容以下情况：
    - 标准 CSV：a,b,c
    - 你这种：a, b, c（逗号后有空格）
    - Date 列可解析为 datetime
    """
    df = pd.read_csv(
        path,
        sep=r"\s*,\s*",     # 逗号两侧允许空白
        engine="python"
    )

    # 基本列检查
    required = {"Date", "Algorithm", "Duration"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"缺少必要列: {missing}. 当前列: {list(df.columns)}")

    # 类型处理
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    df["Duration"] = pd.to_numeric(df["Duration"], errors="coerce")

    # 丢掉 Duration 为空的行
    df = df.dropna(subset=["Duration"])

    return df


def summarize_by_algorithm(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby("Algorithm")["Duration"]

    # 你要的“平均值、方差等等”
    summary = g.agg(
        n="count",
        mean="mean",
        var="var",       # 样本方差（ddof=1）
        std="std",       # 样本标准差（ddof=1）
        min="min",
        q10=lambda s: s.quantile(0.10),
        q25=lambda s: s.quantile(0.25),
        median="median",
        q75=lambda s: s.quantile(0.75),
        q90=lambda s: s.quantile(0.90),
        max="max",
    )

    # 偏度/峰度（pandas 默认是样本偏度/超额峰度）
    summary["skew"] = g.skew()
    summary["kurtosis"] = g.apply(pd.Series.kurt)

    # 变异系数（std/mean），均值为0时避免除零
    summary["cv"] = summary["std"] / summary["mean"].replace(0, pd.NA)

    # 让输出更好看
    summary = summary.sort_values("mean", ascending=True)

    return summary.reset_index()


def main():
    parser = argparse.ArgumentParser(description="按 Algorithm 统计 Duration 运行时间参数")
    parser.add_argument("input_csv", help="输入CSV路径")
    parser.add_argument("-o", "--output", default="", help="输出统计结果CSV路径（可选）")
    args = parser.parse_args()

    df = load_csv(args.input_csv)
    summary = summarize_by_algorithm(df)

    # 控制台打印
    pd.set_option("display.max_rows", 200)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 50)
    print(summary.to_string(index=False))

    # 可选导出
    if args.output:
        summary.to_csv(args.output, index=False)
        print(f"\n已写出: {args.output}")


if __name__ == "__main__":
    main()