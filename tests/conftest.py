"""
conftest.py — shared pytest fixtures for the RNA-seq pipeline test suite.

All fixtures are file-system free by default (in-memory DataFrames / tmp_path).
Heavy I/O fixtures use tmp_path so they are cleaned up automatically.
"""

import json
import textwrap
from pathlib import Path

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# OUTRIDER fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def outrider_df():
    """Minimal OUTRIDER result table with two samples and three genes."""
    return pd.DataFrame({
        "geneID":          ["ENSG00000001", "ENSG00000001", "ENSG00000002",
                            "ENSG00000002", "ENSG00000003"],
        "sampleID":        ["25D1001-HOL-Hay", "25D1002-HOL-Hay",
                            "25D1001-HOL-Hay", "25D1002-HOL-Hay",
                            "25D1001-HOL-Hay"],
        "pValue":          [0.001, 0.5,   0.002, 0.8,   0.03],
        "padjust":         [0.01,  0.9,   0.02,  0.95,  0.1],
        "zScore":          [-3.2,  0.1,   4.5,   -0.2,  -2.1],
        "l2fc":            [-1.5,  0.05,  2.1,   -0.1,  -0.9],
        "rawcounts":       [5,     120,   300,   110,   8],
        "normcounts":      [4.8,   118.0, 295.0, 108.0, 7.5],
        "meanCorrected":   [100.0, 100.0, 100.0, 100.0, 100.0],
        "aberrant":        [True,  False, True,  False, True],
        "AberrantBySample":[1,     0,     1,     0,     1],
        "AberrantByGene":  [2,     2,     2,     2,     1],
    })


@pytest.fixture
def outrider_file(tmp_path, outrider_df):
    """OUTRIDER TSV file on disk."""
    p = tmp_path / "outrider_htseq.tab"
    outrider_df.to_csv(p, sep="\t", index=False)
    return p


# ---------------------------------------------------------------------------
# FRASER fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def fraser_df():
    """Minimal FRASER result table with two samples."""
    return pd.DataFrame({
        "seqnames":     ["chr1",  "chr2",  "chrX"],
        "start":        [1000,    2000,    3000],
        "end":          [1100,    2100,    3100],
        "width":        [100,     100,     100],
        "strand":       ["+",     "-",     "+"],
        "sampleID":     ["25D1001-HOL", "25D1001-HOL", "25D1002-HOL"],
        "type":         ["DS",    "SS",    "DS"],
        "pValue":       [0.0001,  0.0002,  0.5],
        "padjust":      [0.001,   0.002,   0.8],
        "psiValue":     [0.3,     0.7,     0.5],
        "deltaPsi":     [-0.4,    0.3,     0.01],
        "hgncSymbol":   ["BRCA2", "TP53",  "ACTB"],
        "counts":       [10,      20,      50],
        "totalCounts":  [100,     200,     500],
    })


@pytest.fixture
def fraser_file(tmp_path, fraser_df):
    """FRASER TSV file on disk."""
    p = tmp_path / "fraser.tab"
    fraser_df.to_csv(p, sep="\t", index=False)
    return p


# ---------------------------------------------------------------------------
# GTF fixture
# ---------------------------------------------------------------------------

GTF_CONTENT = textwrap.dedent("""\
    #!genome-build GRCh38
    1\thavana\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000001"; gene_name "BRCA2";
    2\thavana\tgene\t2000\t3000\t.\t-\t.\tgene_id "ENSG00000002"; gene_name "TP53";
    X\thavana\tgene\t3000\t4000\t.\t+\t.\tgene_id "ENSG00000003"; gene_name "ACTB";
""")


@pytest.fixture
def gtf_file(tmp_path):
    """Tiny GTF file with 3 genes."""
    p = tmp_path / "genes.gtf"
    p.write_text(GTF_CONTENT)
    return p


# ---------------------------------------------------------------------------
# Blacklist fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def blacklist_file(tmp_path):
    """Blacklist with one sample."""
    p = tmp_path / "blacklist.txt"
    p.write_text("25D1002 excluded_reason\n")
    return p


# ---------------------------------------------------------------------------
# Mendeliome JSON fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def mendeliome_file(tmp_path):
    """Minimal Mendeliome JSON (Australia format)."""
    payload = {
        "version": "2024.01",
        "genes": [
            {
                "gene_data": {"gene_symbol": "BRCA2"},
                "confidence_level": "3",
                "mode_of_inheritance": "MONOALLELIC",
                "phenotypes": ["Breast cancer", "Ovarian cancer"],
            },
            {
                "gene_data": {"gene_symbol": "TP53"},
                "confidence_level": "3",
                "mode_of_inheritance": "MONOALLELIC",
                "phenotypes": ["Li-Fraumeni syndrome"],
            },
        ],
    }
    p = tmp_path / "mendeliome.json"
    p.write_text(json.dumps(payload))
    return p


# ---------------------------------------------------------------------------
# gnomAD fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def gnomad_df():
    """Minimal gnomAD v2-style constraint table."""
    return pd.DataFrame({
        "gene":  ["BRCA2", "TP53", "ACTB"],
        "pLI":   [0.9,     0.99,   0.01],
        "oe_lof":[0.12,    0.08,   0.85],
        "lof_z": [3.1,     4.2,    -0.5],
        "mis_z": [2.0,     3.0,     0.1],
        "syn_z": [0.5,     1.0,    -0.2],
        "oe_mis":[0.7,     0.6,     1.1],
        "oe_syn":[1.0,     0.9,     1.0],
    })


@pytest.fixture
def gnomad_file(tmp_path, gnomad_df):
    """gnomAD TSV on disk."""
    p = tmp_path / "gnomad.tsv"
    gnomad_df.to_csv(p, sep="\t", index=False)
    return p


# ---------------------------------------------------------------------------
# Fraser count directory (backup_fraser_counts)
# ---------------------------------------------------------------------------

@pytest.fixture
def fraser_count_dir(tmp_path):
    """Minimal FRASER count directory structure."""
    base = tmp_path / "fraser_counts"
    for sub in ["nonSplitCounts", "splitCounts"]:
        d = base / "savedObjects" / "Data_Analysis" / sub
        d.mkdir(parents=True)
        (d / "dummy.rds").write_text("dummy")
    return base


# ---------------------------------------------------------------------------
# HTSeq count files (create_matrice_by_run)
# ---------------------------------------------------------------------------

@pytest.fixture
def htseq_run_dir(tmp_path):
    """
    Fake run directory tree expected by create_matrice_by_run.py.

        tmp_path/
          20240101_RUN25_NextSeq_High_8RNASEQ/
            pipeline_v0/
              htseq/
                25D1001-MOINS/
                  25D1001-MOINS_gene_counts.txt
                25D1002-MOINS/
                  25D1002-MOINS_gene_counts.txt
    """
    genes = ["ENSG00000001", "ENSG00000002", "ENSG00000003"]

    for sample, counts in [("25D1001-MOINS", [10, 20, 30]),
                            ("25D1002-MOINS", [5,  15, 25])]:
        d = (tmp_path / "20240101_RUN25_NextSeq_High_8RNASEQ"
             / "pipeline_v0" / "htseq" / sample)
        d.mkdir(parents=True)
        rows = "\n".join(f"{g}\t{c}" for g, c in zip(genes, counts))
        (d / f"{sample}_gene_counts.txt").write_text(rows + "\n")

    return tmp_path
