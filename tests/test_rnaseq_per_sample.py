"""
test_rnaseq_per_sample.py
Tests for scripts/rnaseq_analysis_per_sample.py

Covers:
  - _gtf_to_dict        : GTF DataFrame → picklable dict
  - _gnomad_to_dict     : gnomAD DataFrame → picklable dict
  - _mendeliome_to_dict : Mendeliome DataFrame → picklable dict
  - _process_and_save_sample : annotation + TSV write per sample
  - RNASeqProcessorPerSample.load_fraser    : column aliasing
  - RNASeqProcessorPerSample.load_outrider  : Unnamed column removal, sampleID normalisation
  - RNASeqProcessorPerSample.load_gtf       : GTF parsing
  - RNASeqProcessorPerSample.load_mendeliome: JSON parsing
  - RNASeqProcessorPerSample._filter_data   : mode / pvalue filtering
  - RNASeqProcessorPerSample._get_matched_samples : exact + partial match
  - RNASeqProcessorPerSample.create_zip_archive
"""

import json
import sys
import textwrap
from pathlib import Path

import pandas as pd
import pytest

# Make scripts/ importable
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from rnaseq_analysis_per_sample import (
    RNASeqProcessorPerSample,
    _gnomad_to_dict,
    _gtf_to_dict,
    _mendeliome_to_dict,
    _process_and_save_sample,
)


# =============================================================================
# _gtf_to_dict
# =============================================================================

class TestGtfToDict:
    def test_returns_dict_with_expected_keys(self, gtf_file, tmp_path):
        proc = _make_processor(tmp_path, gtf_file)
        proc.load_gtf()
        result = _gtf_to_dict(proc.gtf_data)
        assert "by_gene" in result
        assert "by_gene_name" in result

    def test_by_gene_indexed_by_clean_ensg(self, gtf_file, tmp_path):
        proc = _make_processor(tmp_path, gtf_file)
        proc.load_gtf()
        d = _gtf_to_dict(proc.gtf_data)
        assert "ENSG00000001" in d["by_gene"]

    def test_by_gene_name_indexed_by_symbol(self, gtf_file, tmp_path):
        proc = _make_processor(tmp_path, gtf_file)
        proc.load_gtf()
        d = _gtf_to_dict(proc.gtf_data)
        assert "BRCA2" in d["by_gene_name"]

    def test_gene_record_contains_coordinates(self, gtf_file, tmp_path):
        proc = _make_processor(tmp_path, gtf_file)
        proc.load_gtf()
        d = _gtf_to_dict(proc.gtf_data)
        rec = d["by_gene"]["ENSG00000001"]
        assert rec["gene_name"] == "BRCA2"
        assert rec["chrom"] == "1"

    def test_none_input_returns_empty_dict(self):
        assert _gtf_to_dict(None) == {}


# =============================================================================
# _gnomad_to_dict
# =============================================================================

class TestGnomadToDict:
    def test_returns_dict_keyed_by_gene(self, gnomad_df):
        d = _gnomad_to_dict(gnomad_df)
        assert "BRCA2" in d
        assert "TP53" in d

    def test_contains_pli(self, gnomad_df):
        d = _gnomad_to_dict(gnomad_df)
        assert d["BRCA2"]["pLI"] == pytest.approx(0.9)

    def test_deduplicates_by_highest_pli(self):
        df = pd.DataFrame({
            "gene": ["BRCA2", "BRCA2"],
            "pLI":  [0.5,     0.9],
            "oe_lof": [0.2, 0.1],
        })
        d = _gnomad_to_dict(df)
        assert len([k for k in d if k == "BRCA2"]) == 1
        assert d["BRCA2"]["pLI"] == pytest.approx(0.9)

    def test_none_returns_empty_dict(self):
        assert _gnomad_to_dict(None) == {}

    def test_missing_optional_columns_still_works(self):
        df = pd.DataFrame({"gene": ["ACTB"], "pLI": [0.01]})
        d = _gnomad_to_dict(df)
        assert "ACTB" in d


# =============================================================================
# _mendeliome_to_dict
# =============================================================================

class TestMendeliomeToDict:
    def test_returns_dict_keyed_by_gene(self):
        df = pd.DataFrame({
            "gene_symbol":         ["BRCA2", "TP53"],
            "confidence_level":    ["3",     "3"],
            "Mode_Of_Inheritance": ["MONO",  "MONO"],
            "Phenotypes":          ["Cancer","LFS"],
        })
        d = _mendeliome_to_dict(df)
        assert "BRCA2" in d

    def test_phenotype_value_preserved(self):
        df = pd.DataFrame({
            "gene_symbol":         ["BRCA2"],
            "confidence_level":    ["3"],
            "Mode_Of_Inheritance": ["MONOALLELIC"],
            "Phenotypes":          ["Breast cancer | Ovarian cancer"],
        })
        d = _mendeliome_to_dict(df)
        assert "Breast cancer" in d["BRCA2"]["Phenotypes"]

    def test_none_returns_empty_dict(self):
        assert _mendeliome_to_dict(None) == {}


# =============================================================================
# _process_and_save_sample
# =============================================================================

class TestProcessAndSaveSample:
    def _gtf_dict_fixture(self):
        gtf_df = pd.DataFrame({
            "gene_id":   ["ENSG00000001", "ENSG00000002"],
            "gene_name": ["BRCA2",        "TP53"],
            "chrom":     ["1",            "2"],
            "start":     [1000,           2000],
            "end":       [2000,           3000],
            "strand":    ["+",            "-"],
        })
        return _gtf_to_dict(gtf_df)

    def test_outrider_file_created(self, tmp_path):
        sample_data = {
            "geneID":   ["ENSG00000001"],
            "sampleID": ["25D1001-HOL"],
            "pValue":   [0.001],
            "padjust":  [0.01],
            "zScore":   [-3.2],
            "l2fc":     [-1.5],
            "rawcounts":[5],
            "normcounts":[4.8],
            "meanCorrected":[100.0],
            "aberrant":  [True],
            "AberrantBySample": [1],
            "AberrantByGene": [1],
        }
        args = (
            "25D1001-HOL",
            sample_data,
            "outrider",
            str(tmp_path),
            self._gtf_dict_fixture(),
            {},   # no gnomad
            {},   # no mendeliome
            "geneID",
        )
        filepath, n, short = _process_and_save_sample(args)
        assert Path(filepath).exists()
        assert n == 1
        assert short == "25D1001"

    def test_fraser_file_created(self, tmp_path):
        sample_data = {
            "seqnames":   ["chr1"],
            "start":      [1000],
            "end":        [1100],
            "width":      [100],
            "strand":     ["+"],
            "sampleID":   ["25D1001-HOL"],
            "type":       ["DS"],
            "pValue":     [0.0001],
            "padjust":    [0.001],
            "psiValue":   [0.3],
            "deltaPsi":   [-0.4],
            "hgncSymbol": ["BRCA2"],
            "counts":     [10],
            "totalCounts":[100],
        }
        args = (
            "25D1001-HOL",
            sample_data,
            "fraser",
            str(tmp_path),
            self._gtf_dict_fixture(),
            {},
            {},
            "hgncSymbol",
        )
        filepath, n, short = _process_and_save_sample(args)
        assert Path(filepath).exists()
        df = pd.read_csv(filepath, sep="\t")
        assert "seqnames" in df.columns
        assert len(df) == 1

    def test_sample_id_correctly_shortened(self, tmp_path):
        sample_data = {"geneID": ["ENSG00000001"], "sampleID": ["25D1001-HOL-Hay"],
                       "pValue": [0.01], "padjust": [0.05], "zScore": [1.0],
                       "l2fc": [0.5], "rawcounts": [50], "normcounts": [48.0],
                       "meanCorrected": [100.0], "aberrant": [False],
                       "AberrantBySample": [0], "AberrantByGene": [0]}
        args = ("25D1001-HOL-Hay", sample_data, "outrider",
                str(tmp_path), {}, {}, {}, "geneID")
        _, _, short = _process_and_save_sample(args)
        assert short == "25D1001"

    def test_gnomad_columns_added(self, tmp_path, gnomad_df):
        gtf_df = pd.DataFrame({
            "gene_id":   ["ENSG00000001"],
            "gene_name": ["BRCA2"],
            "chrom":     ["1"],
            "start":     [1000],
            "end":       [2000],
            "strand":    ["+"],
        })
        sample_data = {
            "geneID":   ["ENSG00000001"],
            "sampleID": ["25D1001"],
            "pValue":   [0.001],
            "padjust":  [0.01],
            "zScore":   [-3.0],
            "l2fc":     [-1.0],
            "rawcounts":[5],
            "normcounts":[4.0],
            "meanCorrected":[100.0],
            "aberrant":  [True],
            "AberrantBySample": [1],
            "AberrantByGene": [1],
        }
        args = (
            "25D1001",
            sample_data,
            "outrider",
            str(tmp_path),
            _gtf_to_dict(gtf_df),
            _gnomad_to_dict(gnomad_df),
            {},
            "geneID",
        )
        filepath, _, _ = _process_and_save_sample(args)
        df = pd.read_csv(filepath, sep="\t")
        assert "pLI" in df.columns
        assert df["pLI"].iloc[0] == pytest.approx(0.9)


# =============================================================================
# RNASeqProcessorPerSample — load_fraser
# =============================================================================

class TestLoadFraser:
    def test_loads_basic_file(self, fraser_file, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, fraser_file=fraser_file)
        df = proc.load_fraser()
        assert df is not None
        assert len(df) == 3

    def test_renames_padjvalue_to_padjust(self, tmp_path, gtf_file):
        df = pd.DataFrame({
            "sampleID":   ["25D1001"],
            "pValue":     [0.001],
            "padjValue":  [0.01],    # old column name
            "hgncSymbol": ["BRCA2"],
        })
        p = tmp_path / "fraser_alias.tab"
        df.to_csv(p, sep="\t", index=False)
        proc = _make_processor(tmp_path, gtf_file, fraser_file=p)
        result = proc.load_fraser()
        assert "padjust" in result.columns
        assert "padjValue" not in result.columns

    def test_none_when_no_file(self, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file)
        assert proc.load_fraser() is None


# =============================================================================
# RNASeqProcessorPerSample — load_outrider
# =============================================================================

class TestLoadOutrider:
    def test_loads_basic_file(self, outrider_file, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, outrider_file=outrider_file)
        df = proc.load_outrider()
        assert df is not None
        assert "sampleID" in df.columns

    def test_removes_unnamed_r_index_column(self, tmp_path, gtf_file):
        df = pd.DataFrame({
            "Unnamed: 0": [0, 1],
            "sampleID":   ["25D1001", "25D1002"],
            "geneID":     ["ENSG00000001", "ENSG00000001"],
            "pValue":     [0.001, 0.5],
        })
        p = tmp_path / "outrider_r.tab"
        df.to_csv(p, sep="\t", index=False)
        proc = _make_processor(tmp_path, gtf_file, outrider_file=p)
        result = proc.load_outrider()
        assert not any(c.startswith("Unnamed") for c in result.columns)

    def test_raises_if_no_sampleid(self, tmp_path, gtf_file):
        df = pd.DataFrame({"geneID": ["ENSG00000001"], "pValue": [0.001]})
        p = tmp_path / "bad.tab"
        df.to_csv(p, sep="\t", index=False)
        proc = _make_processor(tmp_path, gtf_file, outrider_file=p)
        with pytest.raises(ValueError, match="sampleID"):
            proc.load_outrider()

    def test_normalises_sampleid_case(self, tmp_path, gtf_file):
        df = pd.DataFrame({
            "sampleid": ["25D1001"],   # lowercase
            "geneID":   ["ENSG00000001"],
            "pValue":   [0.001],
        })
        p = tmp_path / "outrider_lc.tab"
        df.to_csv(p, sep="\t", index=False)
        proc = _make_processor(tmp_path, gtf_file, outrider_file=p)
        result = proc.load_outrider()
        assert "sampleID" in result.columns


# =============================================================================
# RNASeqProcessorPerSample — load_gtf
# =============================================================================

class TestLoadGtf:
    def test_parses_correct_number_of_genes(self, gtf_file, tmp_path):
        proc = _make_processor(tmp_path, gtf_file)
        df = proc.load_gtf()
        assert len(df) == 3

    def test_skips_comment_lines(self, gtf_file, tmp_path):
        proc = _make_processor(tmp_path, gtf_file)
        df = proc.load_gtf()
        assert not any(str(v).startswith("#") for v in df["gene_id"])

    def test_contains_expected_columns(self, gtf_file, tmp_path):
        proc = _make_processor(tmp_path, gtf_file)
        df = proc.load_gtf()
        for col in ("chrom", "start", "end", "strand", "gene_id", "gene_name"):
            assert col in df.columns, f"Missing column: {col}"

    def test_gene_name_parsed_correctly(self, gtf_file, tmp_path):
        proc = _make_processor(tmp_path, gtf_file)
        df = proc.load_gtf()
        names = set(df["gene_name"])
        assert "BRCA2" in names
        assert "TP53" in names


# =============================================================================
# RNASeqProcessorPerSample — load_mendeliome
# =============================================================================

class TestLoadMendeliome:
    def test_parses_genes_from_json(self, mendeliome_file, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, mendeliome_file=mendeliome_file)
        df = proc.load_mendeliome()
        assert df is not None
        assert "BRCA2" in df["gene_symbol"].values

    def test_deduplicates_gene_symbols(self, tmp_path, gtf_file):
        payload = {
            "version": "1",
            "genes": [
                {"gene_data": {"gene_symbol": "BRCA2"}, "confidence_level": "3",
                 "mode_of_inheritance": "MONO", "phenotypes": ["A"]},
                {"gene_data": {"gene_symbol": "BRCA2"}, "confidence_level": "2",
                 "mode_of_inheritance": "BI", "phenotypes": ["B"]},
            ],
        }
        p = tmp_path / "m.json"
        p.write_text(json.dumps(payload))
        proc = _make_processor(tmp_path, gtf_file, mendeliome_file=p)
        df = proc.load_mendeliome()
        assert df["gene_symbol"].duplicated().sum() == 0

    def test_returns_none_for_missing_file(self, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file,
                               mendeliome_file=tmp_path / "nonexistent.json")
        result = proc.load_mendeliome()
        assert result is None


# =============================================================================
# RNASeqProcessorPerSample — _filter_data
# =============================================================================

class TestFilterData:
    def test_mode_all_returns_all_samples(self, outrider_df, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, mode="all")
        proc.outrider_data = outrider_df
        result = proc._filter_data(outrider_df.copy(), "OUTRIDER")
        assert len(result) == len(outrider_df)

    def test_pvalue_filter_removes_non_significant(self, outrider_df, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, mode="all", pvalue_filter=0.05)
        result = proc._filter_data(outrider_df.copy(), "OUTRIDER")
        assert all(result["padjust"] < 0.05)

    def test_pvalue_filter_keeps_nothing_when_all_fail(self, outrider_df, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, mode="all", pvalue_filter=0.0001)
        result = proc._filter_data(outrider_df.copy(), "OUTRIDER")
        assert len(result) == 0

    def test_sample_mode_filters_to_listed_samples(self, outrider_df, tmp_path, gtf_file):
        samples_file = tmp_path / "samples.txt"
        samples_file.write_text("25D1001-HOL-Hay\n")
        proc = _make_processor(tmp_path, gtf_file, mode="samples",
                               samples_file=samples_file)
        proc.samples = ["25D1001-HOL-Hay"]
        result = proc._filter_data(outrider_df.copy(), "OUTRIDER")
        assert set(result["sampleID"].unique()) == {"25D1001-HOL-Hay"}


# =============================================================================
# RNASeqProcessorPerSample — _get_matched_samples
# =============================================================================

class TestGetMatchedSamples:
    def test_exact_match_returns_only_listed(self, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, mode="samples")
        proc.samples = ["25D1001-HOL"]
        data_samples = ["25D1001-HOL", "25D1002-HOL"]
        result = proc._get_matched_samples(data_samples)
        assert result == ["25D1001-HOL"]

    def test_partial_match_mode(self, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, mode="samples",
                               partial_match=True)
        proc.samples = ["25D1001"]
        data_samples = ["25D1001-HOL-Hay", "25D1002-HOL-Hay"]
        result = proc._get_matched_samples(data_samples)
        assert "25D1001-HOL-Hay" in result
        assert "25D1002-HOL-Hay" not in result

    def test_all_mode_returns_everything(self, tmp_path, gtf_file):
        proc = _make_processor(tmp_path, gtf_file, mode="all")
        proc.samples = None
        data_samples = ["25D1001-HOL", "25D1002-HOL"]
        result = proc._get_matched_samples(data_samples)
        assert set(result) == set(data_samples)

    def test_warning_for_missing_samples(self, tmp_path, gtf_file, caplog):
        import logging
        proc = _make_processor(tmp_path, gtf_file, mode="samples")
        proc.samples = ["25D9999"]
        with caplog.at_level(logging.WARNING):
            proc._get_matched_samples(["25D1001-HOL"])
        assert "25D9999" in caplog.text


# =============================================================================
# RNASeqProcessorPerSample — create_zip_archive
# =============================================================================

class TestCreateZipArchive:
    def test_zip_created(self, tmp_path, gtf_file):
        # Create some dummy files to zip
        files = []
        for name in ("s1.outrider.tab", "s2.fraser.tab"):
            p = tmp_path / name
            p.write_text("col1\tcol2\nA\tB\n")
            files.append(p)
        proc = _make_processor(tmp_path, gtf_file)
        zip_path = proc.create_zip_archive(files)
        assert zip_path.exists()
        assert zip_path.suffix == ".zip"

    def test_zip_contains_all_files(self, tmp_path, gtf_file):
        import zipfile
        files = []
        for name in ("s1.outrider.tab", "s2.fraser.tab"):
            p = tmp_path / name
            p.write_text("col1\tcol2\nA\tB\n")
            files.append(p)
        proc = _make_processor(tmp_path, gtf_file)
        zip_path = proc.create_zip_archive(files)
        with zipfile.ZipFile(zip_path) as zf:
            names = zf.namelist()
        assert "s1.outrider.tab" in names
        assert "s2.fraser.tab" in names


# =============================================================================
# Helper
# =============================================================================

def _make_processor(tmp_path, gtf_file, fraser_file=None, outrider_file=None,
                    mendeliome_file=None, samples_file=None,
                    mode="all", pvalue_filter=None, partial_match=False):
    """Build a minimal processor without requiring real data files."""
    return RNASeqProcessorPerSample(
        fraser_file=fraser_file,
        outrider_file=outrider_file,
        samples_file=samples_file,
        gtf_file=gtf_file,
        output_dir=tmp_path / "output",
        mendeliome_file=mendeliome_file,
        mode=mode,
        pvalue_filter=pvalue_filter,
        create_zip=False,
        partial_match=partial_match,
        workers=1,
    )
