#!/usr/bin/env python3
"""
fraser_config_create.py
=======================
Generate a FRASER2 sample config TSV from a directory of STAR BAMs.

Output format (TSV):
    sampleID  bamFile  group  pairedEnd

Columns:
  group     — always 1 (whole cohort shares one group).
              In FRASER2 the autoencoder covariation correction runs over ALL
              samples in the FraserDataSet; the group column is a legacy field
              from FRASER1 and is unused by FRASER ≥ 2.x.

  pairedEnd — always TRUE (pipeline produces paired-end STAR BAMs).

Usage:
    python fraser_config_create.py \\
        --output      fraser_config.txt \\
        --pattern     MOINS \\
        --root_dir    /datawork2/.../runs_dir \\
        --fraser_dir  /datawork2/.../pipeline_v0/fraser \\
        [--blacklist  unified_blacklist.tsv] \\
        [--excluded_runs RUN17,RUN23]
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate FRASER2 sample config from STAR BAM directory")
    p.add_argument("--output",         required=True,
                   help="Output config TSV path")
    p.add_argument("--pattern",        required=True,
                   help="Filename pattern to select BAMs (e.g. MOINS)")
    p.add_argument("--root_dir",       required=True,
                   help="Root directory to scan for BAM files")
    p.add_argument("--fraser_dir",     required=True,
                   help="FRASER working directory (passed to backup logic)")
    p.add_argument("--blacklist",      default="",
                   help="Unified blacklist TSV (sample_id<TAB>tool<TAB>reason). "
                        "Samples with tool='fraser' or tool='all' are excluded.")
    p.add_argument("--excluded_runs",  default="",
                   help="Comma-separated list of run tags to exclude "
                        "(e.g. RUN17,RUN23). Replaces the hardcoded RUN17 filter.")
    # Legacy positional-argument compatibility (old CLI: script out pattern root fraser [bl])
    p.add_argument("positional", nargs="*",
                   help=argparse.SUPPRESS)
    return p.parse_args()


def extract_run_number(bam_path):
    """Extract numeric run index from RUNxx in BAM path."""
    basename = os.path.basename(bam_path)
    for part in basename.split("_"):
        if part.startswith("RUN"):
            try:
                return int(part[3:])
            except ValueError:
                continue
    return 0


def load_blacklist(blacklist_file, tool="fraser"):
    """
    Load the unified blacklist TSV and return the set of sample_ids to exclude
    for the given tool. Entries with tool='all' are always excluded.
    Ignores comment lines and the header row.
    """
    excluded = set()
    if not blacklist_file or not os.path.exists(blacklist_file):
        return excluded

    with open(blacklist_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[0] == "sample_id":   # header row
                continue
            sample_id   = parts[0]
            entry_tool  = parts[1].lower() if len(parts) > 1 else "all"
            if entry_tool in (tool, "all"):
                excluded.add(sample_id)

    print(f"[blacklist] fraser: {len(excluded)} sample(s) excluded",
          file=sys.stderr)
    return excluded


def generate_config_fraser(output_filename, pattern, root_dir,
                           fraser_output_dir, blacklist_file=None,
                           excluded_runs=None):
    excluded_runs = set(excluded_runs or [])

    blacklist = load_blacklist(blacklist_file or "", tool="fraser")

    # Collect BAMs — keep most recent run per sample_id
    sample_to_bams = defaultdict(list)
    for dirpath, _dirnames, filenames in os.walk(root_dir):
        if "RNASEQ" not in dirpath or "star" not in dirpath:
            continue
        # Skip explicitly excluded runs
        if any(tag in dirpath for tag in excluded_runs):
            continue
        for filename in filenames:
            if not filename.endswith("Aligned.sortedByCoord.out.bam"):
                continue
            if pattern not in filename:
                continue
            # Derive short sample_id from BAM filename
            if "-" in filename:
                sample_id = filename.split("-")[0]
            elif "_" in filename:
                sample_id = filename.split("_")[0]
            else:
                sample_id = filename
            if sample_id in blacklist:
                continue
            bam_path = os.path.join(dirpath, filename)
            run_num  = extract_run_number(bam_path)
            sample_to_bams[sample_id].append((run_num, bam_path))

    # Keep most recent BAM per sample
    final_samples = {
        sid: sorted(bams, key=lambda x: x[0], reverse=True)[0][1]
        for sid, bams in sample_to_bams.items()
    }

    # Write config — group=1 for all (FRASER2 ignores this column; one cohort)
    os.makedirs(os.path.dirname(os.path.abspath(output_filename)), exist_ok=True)
    with open(output_filename, "w") as out:
        out.write("sampleID\tbamFile\tgroup\tpairedEnd\n")
        for sample_id, bam_path in sorted(final_samples.items()):
            out.write(f"{sample_id}\t{bam_path}\t1\tTRUE\n")

    print(f"[fraser_config] Written to {output_filename} "
          f"— {len(final_samples)} samples included.",
          file=sys.stderr)
    if excluded_runs:
        print(f"[fraser_config] Excluded run tags: {sorted(excluded_runs)}",
              file=sys.stderr)


def main():
    args = parse_args()

    # Legacy positional CLI support: out pattern root fraser [blacklist]
    if args.positional:
        pos = args.positional
        generate_config_fraser(
            output_filename  = pos[0] if len(pos) > 0 else args.output,
            pattern          = pos[1] if len(pos) > 1 else args.pattern,
            root_dir         = pos[2] if len(pos) > 2 else args.root_dir,
            fraser_output_dir= pos[3] if len(pos) > 3 else args.fraser_dir,
            blacklist_file   = pos[4] if len(pos) > 4 else args.blacklist,
            excluded_runs    = [r.strip() for r in args.excluded_runs.split(",") if r.strip()],
        )
        return

    generate_config_fraser(
        output_filename   = args.output,
        pattern           = args.pattern,
        root_dir          = args.root_dir,
        fraser_output_dir = args.fraser_dir,
        blacklist_file    = args.blacklist,
        excluded_runs     = [r.strip() for r in args.excluded_runs.split(",") if r.strip()],
    )


if __name__ == "__main__":
    main()


