#!/usr/bin/env python3
"""
create_matrice_featurecounts.py
================================
Assemble per-sample featureCounts ENSG-remapped count files into a single
gene × sample matrix TSV, ready for OUTRIDER.

The input files are produced by map_refseq_to_ensembl and look like:

    # Program:featureCounts ...
    Geneid  Chr  Start  End  Strand  Length  <bam_path>
    ENSG00000001  chr1  ...  45
    ENSG00000002  chr2  ...  12
    ...

The output is a TSV with genes as rows and sample IDs as columns:

    ENSG00000001  25D1001  25D1002  ...
    45            12       ...

Usage:
    python create_matrice_featurecounts.py \\
        --input_dir  /prod/RUN48/pipeline_v0/featureCounts_gencode/ \\
        --output     /prod/RUN48/pipeline_v0/featureCounts_gencode/matrice_fc.txt \\
        --blacklist  /prod/blacklist_outrider.txt \\
        --runs_dir   /datawork2/genetique/RNASeq/diag/prod \\
        [--pattern   "*_gene_cds_counts_ensembl.txt"]
"""

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Assemble featureCounts ENSG matrices for OUTRIDER")
    p.add_argument("--runs_dir",   required=True,
                   help="Root production directory containing all run folders")
    p.add_argument("--output",     required=True,
                   help="Output matrix TSV path")
    p.add_argument("--blacklist",  default="",
                   help="Blacklist file (one sample ID per line, optional)")
    p.add_argument("--pattern",    default="*_gene_cds_counts_ensembl.txt",
                   help="Glob pattern for count files under each sample folder")
    p.add_argument("--run_filter", default="",
                   help="Restrict to a single run folder name (optional)")
    p.add_argument("--pipeline_dir", default="pipeline_v0",
                   help="Pipeline output subdirectory name (default: pipeline_v0)")
    p.add_argument("--keywords",   default="MOINS,PUROMINS",
                   help="Comma-separated keywords — only samples whose folder name "
                        "contains at least one keyword are included (case-insensitive). "
                        "Set to empty string to include all samples.")
    return p.parse_args()


def load_blacklist(path):
    blacklist = set()
    if path and os.path.exists(path):
        with open(path) as f:
            for line in f:
                sample = line.strip().split()[0]
                if sample:
                    blacklist.add(sample)
        print(f"[blacklist] {len(blacklist)} samples excluded", file=sys.stderr)
    return blacklist


def extract_run_number(run_folder):
    m = re.search(r"RUN(\d+)", run_folder, re.IGNORECASE)
    return int(m.group(1)) if m else 0


def shorten_sample_id(sample_id):
    """25D1001-MOINS → 25D1001  (keep first 7 chars / up to first '-')"""
    return sample_id.split("-")[0][:7]


def load_count_file(path):
    """
    Parse a featureCounts output file (possibly ENSG-remapped).
    Returns a Series with ENSG gene IDs as index and counts as values.
    Skips comment lines (starting with #) and the Geneid header.
    """
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            gene_id = parts[0].split(".")[0]   # strip version suffix
            try:
                count = int(parts[-1])
            except ValueError:
                continue
            rows[gene_id] = count
    return pd.Series(rows, name=path)


def main():
    args = parse_args()
    blacklist = load_blacklist(args.blacklist)
    keywords  = [kw.upper() for kw in args.keywords.split(",") if kw.strip()]
    if keywords:
        print(f"[filter] Keeping only samples matching keywords: {keywords}",
              file=sys.stderr)
    else:
        print("[filter] No keyword filter — including all samples", file=sys.stderr)

    # ------------------------------------------------------------------ #
    # Discover all count files across all runs
    # ------------------------------------------------------------------ #
    pattern = os.path.join(
        args.runs_dir,
        args.run_filter if args.run_filter else "20*_RUN*",
        args.pipeline_dir,
        "featureCounts_gencode",
        "*",                  # sample subfolder
        args.pattern
    )

    all_files = sorted(glob.glob(pattern))
    if not all_files:
        print(f"[ERROR] No files found matching: {pattern}", file=sys.stderr)
        sys.exit(1)

    print(f"[info] Found {len(all_files)} count files", file=sys.stderr)

    # ------------------------------------------------------------------ #
    # Group by short sample ID, keeping most recent run
    # ------------------------------------------------------------------ #
    # key: short_sample_id → (run_number, file_path)
    best = {}
    for fpath in all_files:
        parts = Path(fpath).parts
        # extract run folder from path
        run_folder = next(
            (p for p in parts if re.match(r"20\d{6}_RUN\d+", p)), None)
        if run_folder is None:
            continue
        # sample folder is the direct parent of the file
        sample_folder = Path(fpath).parent.name
        short_id = shorten_sample_id(sample_folder)

        if short_id in blacklist:
            print(f"[skip] {short_id} is blacklisted", file=sys.stderr)
            continue

        # The sample subfolder only has the short ID (e.g. 26D0198).
        # Look for keywords in the full path which includes the original
        # sample name from the featureCounts input BAM.
        if keywords:
            # Try to find the original BAM path from the count file header
            # (featureCounts embeds it as a comment on line 1)
            full_name_found = False
            try:
                with open(fpath) as fh:
                    header = fh.readline()   # e.g. # Program:featureCounts ... BAM path
                    if any(kw in header.upper() for kw in keywords):
                        full_name_found = True
            except Exception:
                pass

            if not full_name_found:
                print(f"[skip] {short_id} ({sample_folder}) — no keyword match in BAM header",
                      file=sys.stderr)
                continue

        run_num = extract_run_number(run_folder)
        if short_id not in best or run_num > best[short_id][0]:
            best[short_id] = (run_num, fpath)

    if not best:
        print("[ERROR] No samples remaining after blacklist filtering",
              file=sys.stderr)
        sys.exit(1)

    print(f"[info] Assembling matrix for {len(best)} samples", file=sys.stderr)

    # ------------------------------------------------------------------ #
    # Load all count files and assemble matrix
    # ------------------------------------------------------------------ #
    series_list = []
    for short_id, (run_num, fpath) in sorted(best.items()):
        s = load_count_file(fpath)
        s.name = short_id
        series_list.append(s)
        print(f"  {short_id}: {len(s)} genes  ({fpath})", file=sys.stderr)

    matrix = pd.concat(series_list, axis=1)
    matrix = matrix.fillna(0).astype(int)
    matrix.index.name = "geneID"

    # ------------------------------------------------------------------ #
    # Write output
    # ------------------------------------------------------------------ #
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(out, sep="\t")
    print(f"[done] Matrix written: {matrix.shape[0]} genes × "
          f"{matrix.shape[1]} samples → {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
