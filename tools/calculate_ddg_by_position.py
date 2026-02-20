#!/usr/bin/env python3

import argparse
import glob
import csv
import re
import os
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

AA3_ORDER = [
    "ALA","ARG","ASN","ASP","CYS",
    "GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO",
    "SER","THR","TRP","TYR","VAL"
]

def read_foldx_energies(filename):
    with open(filename) as f:
        for line in f:
            if line.startswith("./"):
                cols = line.split()
                try:
                    interaction = float(cols[5])
                    stability2 = float(cols[7])
                    return interaction, stability2
                except (IndexError, ValueError):
                    raise ValueError(f"{filename} 格式错误：无法读取能量列")
    raise ValueError(f"{filename} 未找到能量行")

def main():
    parser = argparse.ArgumentParser(description="计算蛋白复合物结合能与稳定性变化")
    parser.add_argument("--wt", required=True, help="野生型 Summary 文件路径")
    parser.add_argument("--indir", default=".", help="包含突变体 Summary 文件的目录")
    parser.add_argument("--pattern", default="Summary_*.fxout", help="突变体文件匹配模式")
    parser.add_argument("--outdir", help="输出路径")
    parser.add_argument("--name", help="输出文件名")
    args = parser.parse_args()

    # 读取 WT
    wt_inter, wt_stab2 = read_foldx_energies(args.wt)

    # 获取突变体文件
    all_files = glob.glob(os.path.join(args.indir, args.pattern))
    mut_files = [f for f in all_files if os.path.abspath(f) != os.path.abspath(args.wt)]

    # 数据存储
    binding_ddg = defaultdict(dict)
    stability_ddg = defaultdict(dict)

    for f in sorted(mut_files):
        fname = os.path.basename(f)
        # 提取突变信息
        m = re.search(r"Summary_([A-Z]{3})(\d+)_", fname)
        if not m:
            continue
        aa3, pos = m.groups()
        if aa3 not in AA3_ORDER:
            continue

        inter, stab2 = read_foldx_energies(f)

        # ΔΔG计算
        binding_ddg[pos][aa3] = round(inter - wt_inter, 2)
        stability_ddg[pos][aa3] = round(stab2 - wt_stab2, 2)

    # 写入CSV
    def write_ddg_csv(out_file, data_dict):
        header = ["Position"] + AA3_ORDER
        with open(out_file, "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            for pos in sorted(data_dict.keys(), key=lambda x: int(x)):
                row = [pos]
                for aa in AA3_ORDER:
                    row.append(f"{data_dict[pos][aa]:.2f}" if aa in data_dict[pos] else "")
                writer.writerow(row)
        print(f"已输出 {out_file}")

    write_ddg_csv(f"{args.outdir}/{args.name}_binding_ddg.csv", binding_ddg)
    write_ddg_csv(f"{args.outdir}/{args.name}_stability_ddg.csv", stability_ddg)

if __name__ == "__main__":
    main()
