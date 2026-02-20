#!/usr/bin/env python3
import os
import re
import csv
import argparse

parser = argparse.ArgumentParser(
    description="Collect PLIP interaction counts and Rosetta dG_separated, then compute MT-WT deltas (no pandas). "
                "Also convert mutation name from <Chain><Pos><MutAA> (e.g., A24R) to <WTAA><Pos><MutAA> (e.g., Q24R) "
                "using a chain residue map file like: 24,Q"
)

parser.add_argument("-p", "--plip_dir", required=True, help="Directory containing PLIP csv files")
parser.add_argument("-s", "--score_file", required=True, help="Rosetta score.sc file")
parser.add_argument("-n", "--wt_name", default="6M0J", help="WT structure name (e.g., 6M0J)")
parser.add_argument("-m", "--wt_residue_map", required=True,
                    help="Residue map file: position,WT_AA1 per line, e.g. PLIP_6M0J_chainA_residues.txt")
parser.add_argument("-o", "--out", default="interaction_dG_summary.csv", help="Output csv filename")

args = parser.parse_args()

PLIP_DIR = args.plip_dir
SCORE_FILE = args.score_file
WT_NAME = args.wt_name
WT_MAP_FILE = args.wt_residue_map
OUT_FILE = args.out

INTERACTION_TYPES = [
    "hydrogen_bonds",
    "hydrophobic_interactions",
    "salt_bridges",
    "pi_stacks",
    "pi_cation_interactions",
    "halogen_bonds",
]

def load_wt_residue_map(fname):
    wt_map = {}
    with open(fname, newline="") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [x.strip() for x in line.split(",")]
            if len(parts) != 2:
                continue
            try:
                pos = int(parts[0])
            except ValueError:
                continue
            aa1 = parts[1].upper()
            if aa1:
                wt_map[pos] = aa1
    if not wt_map:
        raise ValueError(f"WT residue map '{fname}' is empty or invalid.")
    return wt_map

WT_MAP = load_wt_residue_map(WT_MAP_FILE)

_mut_pat = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")

def convert_mutation_name(raw_mut):
    """
    Convert <Chain><Pos><MutAA> to <WTAA><Pos><MutAA>.
    Example: A24R + WT_MAP[24]=Q => Q24R.
    If pattern doesn't match or pos not found, return raw_mut unchanged.
    """
    m = _mut_pat.match(raw_mut)
    if not m:
        return raw_mut
    _chain, pos_str, mut_aa = m.groups()
    pos = int(pos_str)
    wt_aa = WT_MAP.get(pos)
    if not wt_aa:
        return raw_mut
    return f"{wt_aa}{pos}{mut_aa.upper()}"

interaction_counts = {}

for fname in os.listdir(PLIP_DIR):
    if not fname.endswith(".csv"):
        continue

    for itype in INTERACTION_TYPES:
        if fname.endswith(f"{itype}.csv"):
            core = fname.replace("PLIP_", "").replace(f"_{itype}.csv", "")
            parts = core.split("_")

            if len(parts) == 1:
                mutation = parts[0]
            else:
                raw_mut = parts[1] 
                mutation = convert_mutation_name(raw_mut)

            path = os.path.join(PLIP_DIR, fname)

            with open(path, newline="") as f:
                reader = csv.reader(f)
                next(reader, None) 
                n = sum(1 for _ in reader)

            interaction_counts.setdefault(
                mutation, {k: 0 for k in INTERACTION_TYPES}
            )[itype] = n

            break

dg_dict = {}

with open(SCORE_FILE) as f:
    for line in f:
        if not line.startswith("SCORE:"):
            continue

        parts = line.split()
        if len(parts) < 2:
            continue
        if parts[1] == "total_score":
            continue

        try:
            dG_sep = float(parts[5])
        except Exception:
            continue

        description = parts[-1]
        desc_parts = description.split("_")

        if len(desc_parts) == 2:
            mutation = desc_parts[0]  
        else:
            raw_mut = desc_parts[1]    
            mutation = convert_mutation_name(raw_mut)

        dg_dict[mutation] = dG_sep

rows = []

def make_row(mutation):
    row = {"mutation": mutation}
    for k in INTERACTION_TYPES:
        row[k] = interaction_counts.get(mutation, {}).get(k, 0)
    row["dG_separated"] = dg_dict.get(mutation)
    return row

if WT_NAME in interaction_counts or WT_NAME in dg_dict:
    rows.append(make_row(WT_NAME))
else:
    rows.append({"mutation": WT_NAME, **{k: 0 for k in INTERACTION_TYPES}, "dG_separated": dg_dict.get(WT_NAME)})

for mut in sorted(interaction_counts):
    if mut == WT_NAME:
        continue
    rows.append(make_row(mut))

for mut in sorted(dg_dict):
    if mut == WT_NAME:
        continue
    if mut not in interaction_counts:
        rows.append(make_row(mut))

wt_row = rows[0]

for row in rows:
    for k in INTERACTION_TYPES:
        row[f"delta_{k}"] = row.get(k, 0) - wt_row.get(k, 0)

    if row.get("dG_separated") is not None and wt_row.get("dG_separated") is not None:
        row["delta_dG_separated"] = round(
            row["dG_separated"] - wt_row["dG_separated"], 3
        )
    else:
        row["delta_dG_separated"] = None

def sort_key(row):
    if row["mutation"] == WT_NAME:
        return (-1, 0, 0)

    m = re.search(r"\d+", row["mutation"])
    site = int(m.group()) if m else float("inf")

    return (
        0,
        -(row["delta_dG_separated"] or 0),
        site
    )

rows = sorted(rows, key=sort_key)

fieldnames = (
    ["mutation"]
    + INTERACTION_TYPES
    + ["dG_separated"]
    + [f"delta_{k}" for k in INTERACTION_TYPES]
    + ["delta_dG_separated"]
)

with open(OUT_FILE, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print("Done!")
