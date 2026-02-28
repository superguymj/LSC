#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd


def load_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s*,\s*",
        engine="python"
    )

    required = {"Date", "Algorithm", "Instance", "Duration"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"缺少必要列: {missing}. 当前列: {list(df.columns)}")

    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    df["Duration"] = pd.to_numeric(df["Duration"], errors="coerce")
    df = df.dropna(subset=["Duration"])
    return df


def summarize_by_algorithm_instance(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby(["Algorithm", "Instance"])["Duration"]

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

    summary["skew"] = g.skew()
    summary["kurtosis"] = g.apply(pd.Series.kurt)

    summary["cv"] = summary["std"] / summary["mean"].replace(0, pd.NA)

    # 排序：先按 Algorithm/Instance，再按 mean（你也可以按需要改）
    summary = summary.sort_values(["Algorithm", "Instance", "mean"], ascending=[True, True, True])

    return summary.reset_index()


def main():
    parser = argparse.ArgumentParser(description="按 Algorithm+Instance 统计 Duration 运行时间参数")
    parser.add_argument("input_csv", help="输入CSV路径")
    parser.add_argument("-o", "--output", default="", help="输出统计结果CSV路径（可选）")
    args = parser.parse_args()

    df = load_csv(args.input_csv)
    summary = summarize_by_algorithm_instance(df)

    pd.set_option("display.max_rows", 500)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 80)
    print(summary.to_string(index=False))

    if args.output:
        summary.to_csv(args.output, index=False)
        print(f"\n已写出: {args.output}")


if __name__ == "__main__":
    main()