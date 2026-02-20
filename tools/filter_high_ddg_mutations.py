#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import glob
import sys
import numpy as np

def find_file_by_pattern(directory, keyword):
    pattern = os.path.join(directory, f"*{keyword}*.csv")
    files = glob.glob(pattern)
    if not files:
        print(f"未找到包含 '{keyword}' 的 CSV 文件，请检查路径。", file=sys.stderr)
        sys.exit(1)
    if len(files) > 1:
        print(f"发现多个文件包含 '{keyword}'，将使用第一个：{files[0]}")
    return files[0]

def main():
    parser = argparse.ArgumentParser(description="筛选 binding_ddG > 阈值 且 stability_ddG < 阈值 的突变并计算加权score")
    parser.add_argument("--dir", required=True, help="包含 ddG CSV 文件的目录")
    parser.add_argument("--chain", required=True, help="突变所在的链名（例如 G）")
    parser.add_argument("--outname", default="filtered_ddg_mutations.csv", help="输出 CSV 文件名")
    parser.add_argument("--bind_threshold", type=float, default=0, help="binding_ddG 阈值（默认 0.5）")
    parser.add_argument("--stab_threshold", type=float, default=1.5, help="stability_ddG 阈值（默认 0.5）")
    parser.add_argument("--w_binding", type=float, default=0.5, help="binding_ddG 权重（默认 0.5）")
    parser.add_argument("--w_stability", type=float, default=0.5, help="stability_ddG 权重（默认 0）")
    args = parser.parse_args()

    binding_file = find_file_by_pattern(args.dir, "binding_ddg")
    stability_file = find_file_by_pattern(args.dir, "stability_ddg")

    print(f"使用文件：\n  binding: {binding_file}\n  stability: {stability_file}")

    bind_df = pd.read_csv(binding_file)
    stab_df = pd.read_csv(stability_file)

    if "Position" not in bind_df.columns or "Position" not in stab_df.columns:
        raise ValueError("输入文件中未找到 'Position' 列，请检查文件格式。")

    bind_melt = bind_df.melt(id_vars=["Position"], var_name="mut_aa", value_name="binding_ddg")
    stab_melt = stab_df.melt(id_vars=["Position"], var_name="mut_aa", value_name="stability_ddg")
    merged = pd.merge(bind_melt, stab_melt, on=["Position", "mut_aa"], how="inner")

    merged["binding_ddg"] = pd.to_numeric(merged["binding_ddg"], errors="coerce")
    merged["stability_ddg"] = pd.to_numeric(merged["stability_ddg"], errors="coerce")


    filtered = merged[
        (merged["binding_ddg"] > args.bind_threshold) &
        (merged["stability_ddg"] < args.stab_threshold)
    ].copy()

    if filtered.empty:
        print("无符合条件的突变，程序结束。")
        sys.exit(0)

    def zscore(series):
        return (series - series.mean()) / series.std(ddof=0)

    filtered["z_binding"] = zscore(filtered["binding_ddg"])
    filtered["z_stability"] = zscore(filtered["stability_ddg"])

    filtered["score"] = (
        args.w_binding * filtered["z_binding"]
        - args.w_stability * filtered["z_stability"]
    )

    out_df = filtered.assign(chain=args.chain)
    out_df.rename(columns={"Position": "resi"}, inplace=True)
    out_df["resi"] = pd.to_numeric(out_df["resi"], errors="coerce")

    out_df = out_df[["chain", "resi", "mut_aa", "binding_ddg", "stability_ddg", "score"]]
    numeric_cols = ["binding_ddg", "stability_ddg", "score"]
    out_df[numeric_cols] = out_df[numeric_cols].round(2)
    out_df = out_df.sort_values(by="resi").reset_index(drop=True)
    
    out_path = os.path.join(args.dir, args.outname)
    out_df.to_csv(out_path, index=False)

    print(f"已提取 {len(out_df)} 个突变满足条件：binding_ddG > {args.bind_threshold} 且 stability_ddG < {args.stab_threshold}")
    print(f"已执行 z-score 归一化并生成加权综合得分 (binding {args.w_binding} / stability {args.w_stability})")
    print(f"输出文件：{out_path}")

if __name__ == "__main__":
    main()
