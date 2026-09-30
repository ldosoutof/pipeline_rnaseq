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
        --runs_dir   /path/to/runs_dir \\
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
    p.add_argument("--pattern",    default="*_gene_cds_counts*.txt",
                   help="Glob pattern for count files under each sample folder")
    p.add_argument("--run_filter", default="",
                   help="Restrict to a single run folder name (optional)")
    p.add_argument("--pipeline_dir", default="pipeline_v0",
                   help="Pipeline output subdirectory name (default: pipeline_v0)")
    p.add_argument("--keywords",   default="MOINS,PUROMOINS",
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
    """
    Extrait l'ID court depuis n'importe quel format :
      26D0643                              -> 26D0643
      26D0643-MOINS                        -> 26D0643
      25D2693-STEMC-PUROMOINS-AVITI        -> 25D2693
      25D2693.STEMC.PUROMOINS.AVITI        -> 25D2693  (format R)
      X26D0643                             -> 26D0643  (ancien R)
    """
    import re as _re
    s = str(sample_id).lstrip('X')
    return _re.split(r'[.\-]', s)[0]


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
    # Discover all count files across ALL runs (courant + historiques).
    # CORRECTIF : le scan utilise TOUJOURS le motif "20*_RUN*" pour agréger
    # toute la cohorte, comme la règle matrix (HTSeq). run_filter ne doit servir
    # QU'À la validation en aval (vérifier que le run courant est présent), PAS à
    # restreindre le scan — sinon la matrice ne contient que le run courant
    # (16 échantillons au lieu de la cohorte complète).
    # ------------------------------------------------------------------ #
    pattern = os.path.join(
        args.runs_dir,
        "20*_RUN*",
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
            # featureCounts embeds the BAM path in its first comment line — search there.
            # Also check the sample folder name itself (covers cases where the subfolder
            # keeps the full name, e.g. 25D2693-STEMC-PUROMOINS-AVITI).
            full_name_found = any(kw in sample_folder.upper() for kw in keywords)
            if not full_name_found:
                try:
                    with open(fpath) as fh:
                        for _ in range(3):          # check first 3 lines (Program, Command, header)
                            line = fh.readline()
                            if any(kw in line.upper() for kw in keywords):
                                full_name_found = True
                                break
                except Exception:
                    pass

            if not full_name_found:
                print(f"[skip] {short_id} ({sample_folder}) — no keyword match",
                      file=sys.stderr)
                continue

        run_num = extract_run_number(run_folder)
        if short_id not in best or run_num > best[short_id][0]:
            best[short_id] = (run_num, fpath, sample_folder)

    if not best:
        print("[ERROR] No samples remaining after blacklist filtering",
              file=sys.stderr)
        sys.exit(1)

    print(f"[info] Assembling matrix for {len(best)} samples", file=sys.stderr)

    # ------------------------------------------------------------------ #
    # Validate current-run samples are all present in the matrix.
    # Historical samples missing on disk (archived runs) are warned and
    # skipped.  But samples from the CURRENT run must never be absent —
    # that would mean featurecounts_gene failed silently.
    # ------------------------------------------------------------------ #
    current_run = args.run_filter   # non-empty only when called for current run
    missing_historical = []
    missing_current    = []

    for sid, (rnum, fp, sample_folder) in list(best.items()):
        if not os.path.exists(fp):
            run_folder = next(
                (p for p in Path(fp).parts if re.match(r"20\d{6}_RUN\d+", p)), ""
            )
            if current_run and run_folder == current_run:
                missing_current.append((sid, fp))
            else:
                missing_historical.append((sid, fp))
            del best[sid]

    # Historical missing → warn only (run may have been archived)
    if missing_historical:
        print(
            f"[WARN] {len(missing_historical)} historical sample(s) have no featureCounts "
            "file on disk (run may be archived) — excluded from matrix:",
            file=sys.stderr
        )
        for sid, fp in missing_historical:
            print(f"  [ARCHIVED] {sid}: {fp}", file=sys.stderr)

    # Current-run missing → hard error: featurecounts_gene must have failed silently
    if missing_current:
        print(
            f"[ERROR] {len(missing_current)} sample(s) from the CURRENT run ({current_run}) "
            "are missing their featureCounts file.\n"
            "  This means featurecounts_gene or map_refseq_to_ensembl did not run for them.\n"
            "  Check Snakemake logs under log/featureCounts/<sample>/",
            file=sys.stderr
        )
        for sid, fp in missing_current:
            print(f"  [MISSING] {sid}: {fp}", file=sys.stderr)
        sys.exit(2)

    if not best:
        print("[ERROR] No samples remaining after removing missing files", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # Load all count files and assemble matrix
    # ------------------------------------------------------------------ #
    series_list = []
    for short_id, (run_num, fpath, sample_folder) in sorted(best.items()):
        s = load_count_file(fpath)
        # Use full sample name from filename instead of short ID
        full_name = os.path.basename(fpath).split("_gene_cds_counts")[0]
        s.name = full_name
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
