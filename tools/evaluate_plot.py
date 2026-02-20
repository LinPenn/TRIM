import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


ENERGY_COL = "delta_dG_separated"

# Energy bar colors
COLOR_DESTABILIZING = "#E64B35"  # delta > 0
COLOR_STABILIZING = "#4DBBD5"    # delta <= 0

INTERACTIONS = [
    ("Hydrogen Bonds", "hydrogen_bonds"),
    ("Hydrophobic", "hydrophobic_interactions"),
    ("Salt Bridges", "salt_bridges"),
    ("Pi-Stacks", "pi_stacks"),
    ("Pi-Cation", "pi_cation_interactions"),
    ("Halogen Bonds", "halogen_bonds"),
]
# Line styles per interaction
STYLE_MAP = {
    "Hydrogen Bonds": dict(color="#1F77B4", marker="o"),
    "Hydrophobic": dict(color="#7F7F7F", marker="s"),
    "Salt Bridges": dict(color="#D62728", marker="^"),
    "Pi-Stacks": dict(color="#9467BD", marker="D"),
    "Pi-Cation": dict(color="#17BECF", marker="v"),
    "Halogen Bonds": dict(color="#8C564B", marker="P"),
}

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate PDF figures (energy barplot + interaction delta lineplot).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-n", "--name", required=True, help="WT label in 'mutation' column (e.g., 6M0J)")
    p.add_argument("-i", "--input", required=True, help="Input CSV path")
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("--font", default="Arial", help="Font family (fallbacks included)")
    return p.parse_args()

def set_style(font_family: str):
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [font_family, "Helvetica", "DejaVu Sans"]
    plt.rcParams["font.size"] = 8
    plt.rcParams["axes.linewidth"] = 1
    plt.rcParams["xtick.major.width"] = 1
    plt.rcParams["ytick.major.width"] = 1
    plt.rcParams["axes.labelsize"] = 10
    plt.rcParams["xtick.labelsize"] = 8
    plt.rcParams["ytick.labelsize"] = 8
    sns.set_context("paper", rc={"lines.linewidth": 1.5})

def plot_energy(df: pd.DataFrame, outdir: str, wt_name: str):
    if ENERGY_COL not in df.columns:
        raise ValueError(f"Missing column: {ENERGY_COL}")

    df_mut = df[df["mutation"] != wt_name].copy()
    if df_mut.empty:
        print("Warning: no mutants found (only WT?). Skip energy plot.")
        return

    df_mut = df_mut.sort_values(ENERGY_COL, ascending=False)

    colors = [
        COLOR_DESTABILIZING if v > 0 else COLOR_STABILIZING
        for v in df_mut[ENERGY_COL].astype(float).values
    ]

    plt.figure(figsize=(10, 4))
    ax = sns.barplot(x="mutation", y=ENERGY_COL, data=df_mut, palette=colors)

    ax.set_xlabel("Mutants", fontweight="bold")
    ax.set_ylabel(r"$\Delta\Delta G_{separated}$", fontweight="bold")

    plt.xticks(rotation=90, fontsize=6)
    plt.axhline(0, color="black", linewidth=0.8, linestyle="--")
    sns.despine()

    out = os.path.join(outdir, "Mutation_Effect.pdf")
    plt.tight_layout()
    plt.savefig(out, format="pdf")
    plt.close()
    print(f"Saved: {out}")

def plot_interaction_delta(df: pd.DataFrame, outdir: str, wt_name: str):

    wt_df = df[df["mutation"] == wt_name]
    if wt_df.empty:
        raise ValueError(f"WT row not found in CSV: mutation == {wt_name}")
    wt_row = wt_df.iloc[0]

    to_plot = []
    for display, base in INTERACTIONS:
        if base in df.columns:
            to_plot.append((display, base))
    if not to_plot:
        raise ValueError("No interaction count columns found in CSV (e.g., hydrogen_bonds, salt_bridges...).")

    if ENERGY_COL not in df.columns:
        raise ValueError(f"Missing column: {ENERGY_COL}")

    df_mut = df[df["mutation"] != wt_name].copy()
    if df_mut.empty:
        print("Warning: no mutants found (only WT?). Skip interaction plot.")
        return

    df_mut = df_mut.sort_values(ENERGY_COL, ascending=False)
    mutants = df_mut["mutation"].astype(str).tolist()

    plt.figure(figsize=(12, 5))

    y_max = 0.0
    for display, col in to_plot:
        style = STYLE_MAP.get(display, dict(color="black", marker="o"))

        y = df_mut[col].astype(float).values
        wt_val = float(wt_row[col])

        y_max = max(y_max, float(np.nanmax(y)), wt_val)

        plt.plot(
            mutants, y,
            label=display,
            color=style["color"],
            marker=style["marker"],
            markersize=3,
            linewidth=1
        )

        plt.axhline(
            wt_val,
            color=style["color"],
            linestyle="--",
            linewidth=1,
            alpha=0.5
        )

        plt.text(len(mutants) - 1, wt_val, " WT", color=style["color"], fontsize=7, va="center")

    plt.xlabel("Mutants", fontweight="bold")
    plt.ylabel("Number of Interactions", fontweight="bold")
    plt.xticks(rotation=90, fontsize=6)

    plt.ylim(0, y_max + 2)

    plt.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=3)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    out = os.path.join(outdir, "Interaction_Analysis.pdf")
    plt.tight_layout()
    plt.savefig(out, format="pdf")
    plt.close()
    print(f"Saved: {out}")

def main():
    args = parse_args()

    if not os.path.exists(args.input):
        print(f"Error: input CSV not found: {args.input}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    set_style(args.font)

    df = pd.read_csv(args.input)

    if "mutation" not in df.columns:
        print("Error: CSV must have a 'mutation' column.")
        sys.exit(1)

    plot_energy(df, args.outdir, args.name)
    plot_interaction_delta(df, args.outdir, args.name)

if __name__ == "__main__":
    main()
