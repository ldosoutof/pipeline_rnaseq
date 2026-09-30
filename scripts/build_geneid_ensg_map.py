#!/usr/bin/env python3
"""
build_geneid_ensg_map.py
========================
Build a GeneID → ENSG mapping TSV from an Ensembl GTF file.

The Ensembl GTF contains lines like:
  gene_id "ENSG00000223972"; ... db_xref "GeneID:100287102"; ...

featureCounts run against a RefSeq GTF uses NCBI GeneIDs as its row index.
This script extracts every (GeneID, ENSG) pair so those counts can be
remapped to Ensembl IDs before being fed into OUTRIDER.

Output format (TSV, no header):
  100287102    ENSG00000223972

Usage:
  python build_geneid_ensg_map.py --gtf /path/to/Homo_sapiens.GRCh38.106.gtf \
                                  --out /datawork/.../geneid_to_ensg.tsv

Or via the Snakemake rule build_geneid_map in 03bis_featureCount.smk.
"""

import argparse
import re
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Build GeneID → ENSG mapping from Ensembl GTF")
    p.add_argument("--gtf", required=True, help="Ensembl GTF file (can be .gtf or .gtf.gz)")
    p.add_argument("--out", required=True, help="Output TSV path (GeneID<TAB>ENSG)")
    p.add_argument("--feature", default="gene",
                   help="GTF feature type to scan (default: gene)")
    return p.parse_args()


def _open_gtf(path):
    """Open plain or gzip-compressed GTF transparently."""
    import gzip
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, encoding="utf-8", errors="replace")


def extract_attr(attr_string, key):
    """
    Extract the value of `key` from a GTF attribute string.
    Handles both quoted  (key "value")  and unquoted  (key value)  forms.
    Returns None if the key is absent.
    """
    m = re.search(r'\b' + re.escape(key) + r'\s+"?([^";]+)"?', attr_string)
    return m.group(1).strip() if m else None


def build_map(gtf_path, feature="gene"):
    """
    Parse the GTF and return a dict  {gene_id_int_str: ensg_id}.

    Strategy:
      • For every `gene` line, extract gene_id (ENSG) and all db_xref "GeneID:XXXXXX".
      • One ENSG can map to multiple GeneIDs (rare) — we keep all pairs.
      • If a GeneID maps to more than one ENSG (extremely rare version conflicts),
        we keep the first one encountered and log a warning.
    """
    mapping = {}          # GeneID_str -> ENSG
    conflicts = {}        # GeneID_str -> [ENSG, ...] for warning
    n_lines = 0
    n_genes = 0
    n_pairs = 0

    geneid_re = re.compile(r'db_xref\s+"GeneID:(\d+)"')

    with _open_gtf(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != feature:
                continue

            attrs = fields[8]
            ensg = extract_attr(attrs, "gene_id")
            if not ensg:
                continue
            # Strip version suffix if present (ENSG00000001.5 → ENSG00000001)
            ensg_clean = ensg.split(".")[0]

            n_genes += 1
            for m in geneid_re.finditer(attrs):
                gene_id = m.group(1)
                if gene_id in mapping:
                    if mapping[gene_id] != ensg_clean:
                        conflicts.setdefault(gene_id, [mapping[gene_id]]).append(ensg_clean)
                else:
                    mapping[gene_id] = ensg_clean
                    n_pairs += 1

    print(f"[build_geneid_map] Scanned {n_lines:,} GTF lines, "
          f"{n_genes:,} gene records, {n_pairs:,} GeneID→ENSG pairs extracted.",
          file=sys.stderr)

    if conflicts:
        print(f"[WARN] {len(conflicts)} GeneID(s) map to multiple ENSGs "
              f"(keeping first encountered):", file=sys.stderr)
        for gid, ensgs in list(conflicts.items())[:10]:
            print(f"  GeneID {gid}: {mapping[gid]} vs {ensgs}", file=sys.stderr)
        if len(conflicts) > 10:
            print(f"  ... and {len(conflicts) - 10} more", file=sys.stderr)

    return mapping


def write_map(mapping, out_path):
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        for gene_id, ensg in sorted(mapping.items(), key=lambda x: int(x[0])):
            fh.write(f"{gene_id}\t{ensg}\n")
    print(f"[build_geneid_map] Wrote {len(mapping):,} pairs to {out_path}",
          file=sys.stderr)


def main():
    args = parse_args()

    if not Path(args.gtf).exists():
        print(f"[ERROR] GTF not found: {args.gtf}", file=sys.stderr)
        sys.exit(1)

    print(f"[build_geneid_map] Parsing {args.gtf} ...", file=sys.stderr)
    mapping = build_map(args.gtf, feature=args.feature)

    if not mapping:
        print("[ERROR] No GeneID→ENSG pairs found. "
              "Check that --gtf is an Ensembl GTF containing db_xref \"GeneID:...\" attributes.",
              file=sys.stderr)
        sys.exit(1)

    write_map(mapping, args.out)


if __name__ == "__main__":
    main()
