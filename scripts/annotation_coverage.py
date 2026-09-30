#!/usr/bin/env python3
"""
annotation_coverage.py
======================
QC : vérifie que les gènes du panel DI (PanelApp DI-Green) sont couverts
par les deux GTF utilisés dans le pipeline RNA-seq.

Sorties (TSV) :
  --out_summary  : tableau récapitulatif par GTF (n_total, n_found, n_missing, pct_coverage)
  --out_missing  : liste des gènes DI absents de chaque GTF (pour inspection)

Usage :
    python annotation_coverage.py \\
        --panelapp  /dataref/.../DI_green.tsv \\
        --gtf       /dataref/.../Homo_sapiens.GRCh38.106.gtf \\
        --gtf_refseq /dataref/.../gencode.v49.basic.annotation.nochr.gtf \\
        --geneid_map /datawork/.../geneid_to_ensg.tsv \\
        --out_summary  annotation_coverage_summary.tsv \\
        --out_missing  annotation_coverage_missing.tsv

Pour lancer depuis le pipeline Snakemake :
    python {PIPELINE_DIR}/scripts/annotation_coverage.py \\
        --panelapp   {PANELAPP} \\
        --gtf        {GTF} \\
        --gtf_refseq {GTF_REFSEQ} \\
        --geneid_map {GENEID_MAP} \\
        --out_summary {output.summary} \\
        --out_missing {output.missing}
"""

import argparse
import gzip
import re
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="QC coverage of DI panel genes in pipeline annotation files")
    p.add_argument("--panelapp",    required=True,
                   help="PanelApp DI-Green TSV (gene_symbol or ENSG column)")
    p.add_argument("--gtf",         required=True,
                   help="Ensembl GTF (pipeline normal — HTSeq/Kallisto/OUTRIDER)")
    p.add_argument("--gtf_refseq",  required=True,
                   help="GENCODE or RefSeq GTF (pipeline hyper — featureCounts CDS)")
    p.add_argument("--geneid_map",  default="",
                   help="GeneID→ENSG TSV from build_geneid_ensg_map.py (optional)")
    p.add_argument("--out_summary", default="annotation_coverage_summary.tsv")
    p.add_argument("--out_missing", default="annotation_coverage_missing.tsv")
    return p.parse_args()


def _open(path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, encoding="utf-8", errors="replace")


def load_panelapp_genes(path):
    """
    Returns a set of ENSG IDs from a PanelApp TSV.
    Accepts files where the ENSG column is named 'EnsemblId(GRch38)',
    'gene_id', 'ensembl_id', or similar; falls back to gene_symbol if absent.
    """
    ensg_set = set()
    sym_set  = set()
    ensg_col = None
    sym_col  = None

    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for i, h in enumerate(header):
            hl = h.lower()
            if "ensembl" in hl and ensg_col is None:
                ensg_col = i
            if ("symbol" in hl or hl == "gene") and sym_col is None:
                sym_col = i

        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if ensg_col is not None and ensg_col < len(cols):
                ensg = cols[ensg_col].split(".")[0].strip()
                if ensg.startswith("ENSG"):
                    ensg_set.add(ensg)
            if sym_col is not None and sym_col < len(cols):
                sym = cols[sym_col].strip()
                if sym:
                    sym_set.add(sym)

    print(f"[panelapp] {len(ensg_set)} ENSG IDs, {len(sym_set)} gene symbols loaded",
          file=sys.stderr)
    return ensg_set, sym_set


def extract_gene_ids_from_gtf(path):
    """
    Returns a set of clean ENSG IDs from a GTF file.
    Handles:
      - Ensembl GTF : gene_id "ENSG..."
      - GENCODE GTF : gene_id "ENSG....version"
      - RefSeq GTF  : gene_id numeric → ignored (no ENSG to extract)
    Also returns a set of gene_names (gene_name attribute) for symbol fallback.
    """
    ensg_set = set()
    name_set = set()
    gene_id_re   = re.compile(r'gene_id\s+"([^"]+)"')
    gene_name_re = re.compile(r'gene_name\s+"([^"]+)"')

    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            m = gene_id_re.search(attrs)
            if m:
                gid = m.group(1).split(".")[0]
                if gid.startswith("ENSG"):
                    ensg_set.add(gid)
            n = gene_name_re.search(attrs)
            if n:
                name_set.add(n.group(1))

    return ensg_set, name_set


def load_geneid_map(path):
    """Returns set of ENSG values covered by the GeneID→ENSG map."""
    covered = set()
    if not path or not Path(path).exists():
        return covered
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                covered.add(parts[1].split(".")[0])
    return covered


def coverage_report(label, di_ensg, gtf_ensg, gtf_names, di_symbols):
    """Returns (summary_dict, missing_ensg_set)."""
    found_by_ensg   = di_ensg & gtf_ensg
    missing_by_ensg = di_ensg - gtf_ensg

    # Secondary check via gene_name for genes without ENSG in PanelApp
    found_by_sym    = {s for s in di_symbols if s in gtf_names}
    missing_by_sym  = di_symbols - gtf_names if di_symbols else set()

    n_total   = len(di_ensg) if di_ensg else len(di_symbols)
    n_found   = len(found_by_ensg) if di_ensg else len(found_by_sym)
    n_missing = len(missing_by_ensg) if di_ensg else len(missing_by_sym)
    pct       = 100 * n_found / n_total if n_total else 0.0

    summary = {
        "source": label,
        "n_DI_genes": n_total,
        "n_found": n_found,
        "n_missing": n_missing,
        "pct_coverage": round(pct, 2),
    }
    missing = missing_by_ensg if di_ensg else set()
    return summary, missing


def main():
    args = parse_args()

    print("[annotation_coverage] Loading PanelApp DI-Green ...", file=sys.stderr)
    di_ensg, di_symbols = load_panelapp_genes(args.panelapp)

    print(f"[annotation_coverage] Parsing GTF (normal): {args.gtf}", file=sys.stderr)
    gtf_ensg, gtf_names = extract_gene_ids_from_gtf(args.gtf)
    print(f"  → {len(gtf_ensg):,} ENSG IDs, {len(gtf_names):,} gene names", file=sys.stderr)

    print(f"[annotation_coverage] Parsing GTF (hyper):  {args.gtf_refseq}", file=sys.stderr)
    gtf_ref_ensg, gtf_ref_names = extract_gene_ids_from_gtf(args.gtf_refseq)
    print(f"  → {len(gtf_ref_ensg):,} ENSG IDs, {len(gtf_ref_names):,} gene names",
          file=sys.stderr)

    gid_map_ensg = load_geneid_map(args.geneid_map)
    if gid_map_ensg:
        print(f"[annotation_coverage] GeneID map covers {len(gid_map_ensg):,} ENSG IDs",
              file=sys.stderr)

    # ── Compute coverage ─────────────────────────────────────────────────────
    rows = []
    missing_rows = []

    for label, ensg_set, name_set in [
        ("GTF_normal (Ensembl)",  gtf_ensg,     gtf_names),
        ("GTF_hyper (GENCODE/RefSeq)", gtf_ref_ensg, gtf_ref_names),
        ("GeneID_map",            gid_map_ensg, set()),
    ]:
        summ, missing = coverage_report(label, di_ensg, ensg_set, name_set, di_symbols)
        rows.append(summ)
        for ensg in sorted(missing):
            missing_rows.append({"source": label, "missing_ENSG": ensg})

    # ── Write summary ─────────────────────────────────────────────────────────
    Path(args.out_summary).parent.mkdir(parents=True, exist_ok=True)
    cols = ["source", "n_DI_genes", "n_found", "n_missing", "pct_coverage"]
    with open(args.out_summary, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    # ── Write missing genes ───────────────────────────────────────────────────
    with open(args.out_missing, "w") as fh:
        fh.write("source\tmissing_ENSG\n")
        for row in missing_rows:
            fh.write(f"{row['source']}\t{row['missing_ENSG']}\n")

    print(f"\n[annotation_coverage] Summary written to {args.out_summary}",
          file=sys.stderr)
    print(f"[annotation_coverage] Missing genes written to {args.out_missing}",
          file=sys.stderr)

    # ── Warn if coverage drops below threshold ────────────────────────────────
    MIN_COVERAGE = 95.0
    for row in rows:
        if row["pct_coverage"] < MIN_COVERAGE:
            print(f"[WARN] {row['source']}: coverage {row['pct_coverage']}% "
                  f"< {MIN_COVERAGE}% threshold ({row['n_missing']} genes missing)",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
