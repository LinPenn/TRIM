import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import argparse
import sys
import os
import glob

# 阈值设置
COLOR_LIMIT = 5.0          # 颜色映射范围 (-5 到 5)
SIZE_MIN = 10              # 最小气泡大小
SIZE_MAX = 400             # 最大气泡大小

def parse_args():
    parser = argparse.ArgumentParser(description="绘制气泡热图 (Bubble Heatmap) - 自动扫描目录")
    parser.add_argument('-d', '--directory', type=str, required=True,
                        help='包含 binding 和 stability CSV 文件的目录路径')
    parser.add_argument('-o', '--output', type=str, default='Bubble_Heatmap.pdf',
                        help='输出图片文件名')
    
    return parser.parse_args()

def find_input_files(directory):
    """在指定目录中自动查找包含 binding 和 stability 的 csv 文件"""
    if not os.path.isdir(directory):
        print(f"错误: 目录不存在 -> {directory}")
        sys.exit(1)
        
    # 获取目录下所有csv文件
    csv_files = glob.glob(os.path.join(directory, "*.csv"))
    
    binding_file = None
    stability_file = None
    
    for f in csv_files:
        filename = os.path.basename(f).lower()
        if "binding" in filename:
            binding_file = f
        elif "stability" in filename:
            stability_file = f
            
    # 检查是否找到
    if not binding_file:
        print(f"错误: 在 '{directory}' 下未找到文件名包含 'binding' 的 CSV 文件。")
        sys.exit(1)
    if not stability_file:
        print(f"错误: 在 '{directory}' 下未找到文件名包含 'stability' 的 CSV 文件。")
        sys.exit(1)
        
    print(f"已自动识别 Binding 文件: {os.path.basename(binding_file)}")
    print(f"已自动识别 Stability 文件: {os.path.basename(stability_file)}")
    
    return binding_file, stability_file

def main():
    args = parse_args()
    
    binding_file, stability_file = find_input_files(args.directory)
    
    output_path = args.output
    if os.path.dirname(output_path) == '':
        output_path = os.path.join(args.directory, output_path)

    try:
        df_binding = pd.read_csv(binding_file)
        df_stability = pd.read_csv(stability_file)
    except Exception as e:
        print(f"读取 CSV 文件时出错: {e}")
        return

    # 确保 Position 是整数
    df_binding['Position'] = df_binding['Position'].astype(int)
    df_stability['Position'] = df_stability['Position'].astype(int)

    # 数据处理
    def melt_df(df, val_name):
        return df.melt(id_vars='Position', var_name='Amino_Acid', value_name=val_name)

    df_bind_long = melt_df(df_binding, 'Binding_DDG')
    df_stab_long = melt_df(df_stability, 'Stability_DDG')

    # 合并数据
    df_merged = pd.merge(df_bind_long, df_stab_long, on=['Position', 'Amino_Acid'])
    df_merged.dropna(subset=['Binding_DDG', 'Stability_DDG'], inplace=True) 

    aa_order = sorted(df_merged['Amino_Acid'].unique())
    pos_order = sorted(df_merged['Position'].unique(), reverse=True)

    aa_map = {aa: i for i, aa in enumerate(aa_order)}
    pos_map = {p: i for i, p in enumerate(pos_order)}

    df_merged['X_Index'] = df_merged['Amino_Acid'].map(aa_map)
    df_merged['Y_Index'] = df_merged['Position'].map(pos_map)

    # 计算气泡大小
    stab_min = -5
    stab_max = 5
    
    df_merged['Bubble_Size'] = (
        (df_merged['Stability_DDG'] - stab_min) / (stab_max - stab_min)
    ) * (SIZE_MAX - SIZE_MIN) + SIZE_MIN

    # 绘图
    plt.figure(figsize=(14, 12))
    ax = plt.gca()
    cmap = mcolors.LinearSegmentedColormap.from_list(
    "binding_ddg",
    ["#3B4CC0", "white", "#B40426"]
    )
    # 绘制主气泡
    scatter = ax.scatter(
        x=df_merged['X_Index'],
        y=df_merged['Y_Index'],
        s=df_merged['Bubble_Size'],
        c=df_merged['Binding_DDG'], 
        cmap=cmap,          
        edgecolors='gray',        
        linewidth=0.5,
        alpha=0.9,
        vmin=-COLOR_LIMIT, 
        vmax=COLOR_LIMIT,
        zorder=2                  
    )

    ax.set_xticks(range(len(aa_order)))
    ax.set_xticklabels(aa_order, fontsize=12, rotation=0)
    
    ax.set_yticks(range(len(pos_order)))
    ax.set_yticklabels(pos_order, fontsize=14)

    ax.set_xticks(np.arange(len(aa_order) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(pos_order) + 1) - 0.5, minor=True)

    ax.grid(which='minor', color='black', linestyle='-', linewidth=1, alpha=0.7)
    ax.grid(which='major', visible=False) 

    ax.set_xlim(-0.5, len(aa_order) - 0.5)
    ax.set_ylim(-0.5, len(pos_order) - 0.5)

    cbar = plt.colorbar(scatter, ax=ax, fraction=0.025, pad=0.03)
    cbar.set_label(r'$\Delta\Delta G_{binding}$', fontweight='bold', fontsize=14)

    # 选取代表性的数值用于展示图例
    legend_vals = np.array([-5, 0, 2.5, 5])
    legend_sizes = ((legend_vals - stab_min) / (stab_max - stab_min)) * (SIZE_MAX - SIZE_MIN) + SIZE_MIN
    
    legend_elements = []
    for val, size in zip(legend_vals, legend_sizes):
        legend_elements.append(
            plt.scatter([], [], s=size, c='white', edgecolors='gray', label=f'{val}')
        )
    
    ax.legend(handles=legend_elements, title= r'$\Delta\Delta G_{Stability}$', 
              loc='upper right', bbox_to_anchor=(1.12, 1), 
              frameon=False, labelspacing=1.5, borderpad=1, title_fontsize=14)

    plt.title('Bubble Heatmap of Mutational Effects', fontsize=20, fontweight='bold', pad=30)
    plt.xlabel('Mutation', fontsize=14, fontweight='bold', labelpad=10)
    plt.ylabel('Position', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    print(f"绘图完成！图片已保存为: {output_path}")

if __name__ == "__main__":
    main()