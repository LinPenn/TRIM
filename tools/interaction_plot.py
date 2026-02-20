import argparse
import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import matplotlib.path as mpath


INTERACTION_TYPES = {
    "Hydrophobic":    {"suffix": "_hydrophobic_interactions.csv", "color": "#7F7F7F", "alpha": 0.5, "zorder": 1},
    "Hydrogen Bonds": {"suffix": "_hydrogen_bonds.csv",           "color": "#1F77B4", "alpha": 0.85, "zorder": 2},
    "Pi-Stacks":      {"suffix": "_pi_stacks.csv",                "color": "#9467BD", "alpha": 0.85, "zorder": 3},
    "Pi-Cation":      {"suffix": "_pi_cation_interactions.csv",   "color": "#2CA02C", "alpha": 0.85, "zorder": 4},
    "Halogen Bonds":  {"suffix": "_halogen_bonds.csv",            "color": "#BCBD22", "alpha": 0.9, "zorder": 5},
    "Salt Bridges":   {"suffix": "_salt_bridges.csv",             "color": "#D62728", "alpha": 0.9, "zorder": 6},
}

FIG_WIDTH = 7
FIG_HEIGHT = 10
X_REC = 0.2
X_LIG = 0.8
BAR_WIDTH = 0.015

AA3_TO_AA1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H",
    "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
    "TYR":"Y","VAL":"V",
    "HIE":"H","HID":"H","HIP":"H","CYX":"C","MSE":"M",
}

def get_res_num(label):
    """Extract residue number for sorting (e.g., QA42 -> 42)."""
    m = re.search(r"\d+", str(label))
    return int(m.group()) if m else 0

def draw_sigmoid(ax, y1, y2, x1, x2, color, alpha, zorder, lw=1.5):
    verts = [
        (x1, y1),
        (x1 + (x2 - x1) / 2, y1),
        (x1 + (x2 - x1) / 2, y2),
        (x2, y2),
    ]
    codes = [mpath.Path.MOVETO, mpath.Path.CURVE4, mpath.Path.CURVE4, mpath.Path.CURVE4]
    path = mpath.Path(verts, codes)
    patch = patches.PathPatch(
        path,
        facecolor="none",
        edgecolor=color,
        lw=lw,
        alpha=alpha,
        capstyle="round",
        zorder=zorder,
    )
    ax.add_patch(patch)

def build_y_map(items):
    n = len(items)
    if n <= 1:
        return {items[0]: 0.5} if n == 1 else {}
    return {item: 0.05 + 0.90 * (i / (n - 1)) for i, item in enumerate(items)}

def parse_residue_label(label):

    s = str(label).strip()

    m = re.match(r"^([A-Za-z])([A-Za-z])(\d+)$", s)
    if m:
        aa1 = m.group(1).upper()
        chain = m.group(2).upper()
        resnum = int(m.group(3))
        return chain, resnum, aa1

    m = re.match(r"^([A-Za-z]{3})([A-Za-z])(\d+)$", s)
    if m:
        aa3 = m.group(1).upper()
        chain = m.group(2).upper()
        resnum = int(m.group(3))
        aa1 = AA3_TO_AA1.get(aa3, "X")
        return chain, resnum, aa1

    m_chain_num = re.search(r"([A-Za-z])\s*[:_\-]?\s*(\d+)", s)
    if m_chain_num:
        chain = m_chain_num.group(1).upper()
        resnum = int(m_chain_num.group(2))

        m_aa1 = re.search(r"^([A-Za-z])\s*[:_\-]?\s*[A-Za-z]\s*[:_\-]?\s*\d+", s)
        if m_aa1:
            aa1 = m_aa1.group(1).upper()
            return chain, resnum, aa1

        m_aa3 = re.search(r"([A-Za-z]{3})", s)
        if m_aa3:
            aa3 = m_aa3.group(1).upper()
            aa1 = AA3_TO_AA1.get(aa3, "X")
            return chain, resnum, aa1

        return chain, resnum, "X"

    return None, None, None

def export_chain_residues(all_data, chain_id, out_path):

    chain_id = chain_id.upper()
    found = set()

    def consider(label):
        ch, resnum, aa1 = parse_residue_label(label)
        if ch == chain_id and resnum is not None and aa1 is not None:
            found.add((resnum, aa1))

    for item in all_data:
        consider(item["rec"])
        consider(item["lig"])

    if not found:
        print(f"[export] No residues found for chain {chain_id}. Nothing written.")
        return

    ordered = sorted(found, key=lambda x: x[0])

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        for resnum, aa1 in ordered:
            f.write(f"{resnum},{aa1}\n")

    print(f"[export] Wrote {len(ordered)} residues for chain {chain_id} -> {out_path}")

def parse_args():
    p = argparse.ArgumentParser(
        description="Draw multi-interaction map (receptor vs ligand) from PLIP CSVs and optionally export residues for a chain."
    )
    p.add_argument("-p", "--prefix", required=True, help="File prefix")
    p.add_argument("-o", "--out", default="Multi_Interaction.pdf", help="Output figure (.pdf recommended; .png also ok).")
    p.add_argument("--font", default="Arial", help="Font family.")
    p.add_argument("--dpi", type=int, default=600, help="DPI for PNG output (ignored for PDF).")
    p.add_argument("-e", "--export_chain", default=None, help="Chain ID to export residues for (e.g., A).")
    p.add_argument("-eo", "--export_out", default=None, help="Output file for exported residues (default: <prefix>_chain<id>_residues.txt).")
    return p.parse_args()

def main():
    args = parse_args()
    plt.rcParams["font.family"] = args.font

    all_data = []
    all_residues_rec = set()
    all_residues_lig = set()

    print("Reading CSV files...")
    for itype, spec in INTERACTION_TYPES.items():
        fname = args.prefix + spec["suffix"]
        if not os.path.exists(fname):
            print(f"  - {itype}: missing {fname}, skipped.")
            continue

        try:
            df = pd.read_csv(fname)
        except Exception as e:
            print(f"  - {itype}: failed to read {fname}: {e}")
            continue

        if "rec" not in df.columns or "lig" not in df.columns:
            print(f"  - {itype}: skipped {fname} (needs columns: rec, lig).")
            continue

        df = df[["rec", "lig"]].drop_duplicates()
        print(f"  - {itype}: {len(df)} interactions")

        for _, row in df.iterrows():
            r, l = row["rec"], row["lig"]
            all_data.append({"rec": r, "lig": l, "type": itype})
            all_residues_rec.add(r)
            all_residues_lig.add(l)

    if not all_data:
        print("Error: no valid data found.")
        sys.exit(1)

    if args.export_chain:
        out_path = args.export_out
        if not out_path:
            out_path = f"{args.prefix}_chain{args.export_chain.upper()}_residues.txt"
        export_chain_residues(all_data, args.export_chain, out_path)

    rec_sorted = sorted(all_residues_rec, key=get_res_num, reverse=True)
    lig_sorted = sorted(all_residues_lig, key=get_res_num, reverse=True)

    rec_y = build_y_map(rec_sorted)
    lig_y = build_y_map(lig_sorted)

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    all_data.sort(key=lambda x: INTERACTION_TYPES[x["type"]]["zorder"])
    for item in all_data:
        r = item["rec"]
        l = item["lig"]
        itype = item["type"]
        if r not in rec_y or l not in lig_y:
            continue
        spec = INTERACTION_TYPES[itype]
        draw_sigmoid(
            ax,
            rec_y[r],
            lig_y[l],
            X_REC,
            X_LIG,
            color=spec["color"],
            alpha=spec["alpha"],
            zorder=spec["zorder"],
            lw=1.5,
        )

    for name, y in rec_y.items():
        ax.add_patch(
            patches.Rectangle(
                (X_REC - BAR_WIDTH, y - 0.006),
                BAR_WIDTH,
                0.012,
                facecolor="#555555",
                edgecolor="none",
                zorder=10,
            )
        )
        ax.text(X_REC - BAR_WIDTH - 0.02, y, str(name), ha="right", va="center", fontsize=10)

    for name, y in lig_y.items():
        ax.add_patch(
            patches.Rectangle(
                (X_LIG, y - 0.006),
                BAR_WIDTH,
                0.012,
                facecolor="#555555",
                edgecolor="none",
                zorder=10,
            )
        )
        ax.text(X_LIG + BAR_WIDTH + 0.02, y, str(name), ha="left", va="center", fontsize=10)

    ax.text(X_REC, 0.98, "Receptor", ha="center", va="bottom", fontsize=14, fontweight="bold")
    ax.text(X_LIG, 0.98, "Ligand", ha="center", va="bottom", fontsize=14, fontweight="bold")

    legend_handles = []
    appeared = {d["type"] for d in all_data}
    for itype, spec in sorted(INTERACTION_TYPES.items(), key=lambda kv: kv[1]["zorder"]):
        if itype not in appeared:
            continue
        legend_handles.append(mlines.Line2D([], [], color=spec["color"], lw=3, label=itype))

    ax.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=3,
        frameon=False,
        fontsize=11,
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    plt.tight_layout()

    out = args.out
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    ext = os.path.splitext(out)[1].lower()
    if ext == ".png":
        plt.savefig(out, dpi=args.dpi, bbox_inches="tight", transparent=False)
        print(f"Done! Saved PNG: {out} (dpi={args.dpi})")
    else:
        plt.savefig(out, format="pdf", bbox_inches="tight")
        print(f"Done! Saved PDF: {out}")

    plt.show()

if __name__ == "__main__":
    main()
