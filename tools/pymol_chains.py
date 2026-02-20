#!/usr/bin/env python3
from pymol import cmd
import argparse
import os
import sys


def normalize_chains(chains: str) -> str:

    chains = chains.strip()

    if "+" not in chains and "," not in chains:
        # ABC -> A+B+C
        return "+".join(list(chains))
    else:
        # A,B,C -> A+B+C
        return chains.replace(",", "+")

def extract_chains(input_pdb, chains, output_pdb):

    cmd.reinitialize()

    # 支持 A,B 或 A+B
    chains = normalize_chains(chains)
    obj_name = os.path.basename(input_pdb).replace(".pdb", "")

    cmd.load(input_pdb, obj_name)
    cmd.select("sel_chain", f"{obj_name} and chain {chains}")
    cmd.save(output_pdb, "sel_chain")

    print(f"[OK] Extract chains {chains} -> {output_pdb}")


def merge_pdbs(pdb_list, output_pdb):

    cmd.reinitialize()

    for i, pdb in enumerate(pdb_list):
        obj_name = f"obj_{i}"
        cmd.load(pdb, obj_name)

    cmd.sort()
    cmd.save(output_pdb, "all")

    print(f"[OK] Merged {len(pdb_list)} PDBs -> {output_pdb}")


def main():
    parser = argparse.ArgumentParser(description="PyMOL utility: extract chains or merge PDBs")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ---------- extract ----------
    p_extract = subparsers.add_parser("extract", help="Extract specified chains from a complex")
    p_extract.add_argument("-i", "--input", required=True, help="Input complex PDB")
    p_extract.add_argument("-c", "--chains", required=True, help="Chains to extract, e.g. A or A,B or A+B")
    p_extract.add_argument("-o", "--output", required=True, help="Output PDB")

    # ---------- merge ----------
    p_merge = subparsers.add_parser("merge", help="Merge multiple PDBs into one")
    p_merge.add_argument("-i", "--inputs", nargs="+", required=True, help="Input PDB files (space separated)")
    p_merge.add_argument("-o", "--output", required=True, help="Output merged PDB")

    args = parser.parse_args()

    if args.command == "extract":
        extract_chains(args.input, args.chains, args.output)
    elif args.command == "merge":
        merge_pdbs(args.inputs, args.output)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
