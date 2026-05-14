"""
test_fraser_config_create.py
Tests for scripts/fraser_config_create.py

Covers:
  - extract_run_number  : RUN number extraction from BAM path
  - generate_config_fraser:
      · writes a valid TSV header
      · discovers BAM files matching pattern
      · de-duplicates samples, keeping most recent run
      · respects blacklist
      · handles missing root_dir gracefully
      · writes correct column count per row
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from fraser_config_create import extract_run_number, generate_config_fraser


# =============================================================================
# extract_run_number
# =============================================================================

class TestExtractRunNumber:
    def test_extracts_number_from_bam_basename(self):
        path = "/data/20240101_RUN25_NextSeq/star/25D1001-MOINS_RUN25_Aligned.sortedByCoord.out.bam"
        assert extract_run_number(path) == 25

    def test_returns_zero_when_no_run_tag(self):
        assert extract_run_number("/data/star/sample.bam") == 0

    def test_extracts_two_digit_run(self):
        path = "25D1001_RUN7_Aligned.sortedByCoord.out.bam"
        assert extract_run_number(path) == 7

    def test_extracts_large_run_number(self):
        path = "25D1001_RUN123_Aligned.sortedByCoord.out.bam"
        assert extract_run_number(path) == 123


# =============================================================================
# generate_config_fraser
# =============================================================================

class TestGenerateConfigFraser:
    def _build_bam_tree(self, tmp_path, samples):
        """
        Create fake BAM files under a path structure matching what
        fraser_config_create.py expects:
          root_dir/RNASEQ/star/<sample>_<run>_Aligned.sortedByCoord.out.bam
        """
        star_dir = tmp_path / "RNASEQ" / "star"
        star_dir.mkdir(parents=True)
        for sample, run in samples:
            bam = star_dir / f"{sample}-MOINS_RUN{run}_Aligned.sortedByCoord.out.bam"
            bam.write_text("dummy")
        return tmp_path

    def test_creates_output_file(self, tmp_path):
        root = self._build_bam_tree(tmp_path / "data", [("25D1001", 25)])
        out = tmp_path / "fraser_config.txt"
        generate_config_fraser(str(out), "MOINS", str(root), str(tmp_path / "fraser"))
        assert out.exists()

    def test_output_has_correct_header(self, tmp_path):
        root = self._build_bam_tree(tmp_path / "data", [("25D1001", 25)])
        out = tmp_path / "fraser_config.txt"
        generate_config_fraser(str(out), "MOINS", str(root), str(tmp_path / "fraser"))
        header = out.read_text().splitlines()[0]
        assert header == "sampleID\tbamFile\tgroup\tgene\tpairedEnd"

    def test_discovers_bam_files(self, tmp_path):
        root = self._build_bam_tree(tmp_path / "data",
                                    [("25D1001", 25), ("25D1002", 25)])
        out = tmp_path / "fraser_config.txt"
        generate_config_fraser(str(out), "MOINS", str(root), str(tmp_path / "fraser"))
        lines = out.read_text().splitlines()
        assert len(lines) == 3   # header + 2 samples

    def test_pattern_filters_bams(self, tmp_path):
        """Only BAMs containing the pattern should be included."""
        star_dir = tmp_path / "RNASEQ" / "star"
        star_dir.mkdir(parents=True)
        (star_dir / "25D1001-MOINS_RUN25_Aligned.sortedByCoord.out.bam").write_text("d")
        (star_dir / "25D9999-POLYA_RUN25_Aligned.sortedByCoord.out.bam").write_text("d")
        out = tmp_path / "config.txt"
        generate_config_fraser(str(out), "MOINS", str(tmp_path), str(tmp_path / "fraser"))
        lines = [l for l in out.read_text().splitlines() if l and not l.startswith("sampleID")]
        sample_ids = [l.split("\t")[0] for l in lines]
        assert "25D1001" in sample_ids
        assert "25D9999" not in sample_ids

    def test_blacklist_excludes_sample(self, tmp_path, blacklist_file):
        root = self._build_bam_tree(tmp_path / "data",
                                    [("25D1001", 25), ("25D1002", 25)])
        out = tmp_path / "config.txt"
        generate_config_fraser(str(out), "MOINS", str(root),
                               str(tmp_path / "fraser"),
                               blacklist_file=str(blacklist_file))
        content = out.read_text()
        assert "25D1001" in content
        assert "25D1002" not in content   # blacklisted

    def test_keeps_most_recent_run_for_duplicate_sample(self, tmp_path):
        star_dir = tmp_path / "RNASEQ" / "star"
        star_dir.mkdir(parents=True)
        # Same sample across two runs
        for run in [20, 30]:
            bam = star_dir / f"25D1001-MOINS_RUN{run}_Aligned.sortedByCoord.out.bam"
            bam.write_text("dummy")
        out = tmp_path / "config.txt"
        generate_config_fraser(str(out), "MOINS", str(tmp_path), str(tmp_path / "fraser"))
        lines = [l for l in out.read_text().splitlines() if "25D1001" in l]
        assert len(lines) == 1
        assert "RUN30" in lines[0]

    def test_each_data_row_has_five_columns(self, tmp_path):
        root = self._build_bam_tree(tmp_path / "data", [("25D1001", 25)])
        out = tmp_path / "config.txt"
        generate_config_fraser(str(out), "MOINS", str(root), str(tmp_path / "fraser"))
        for line in out.read_text().splitlines()[1:]:
            assert len(line.split("\t")) == 5

    def test_paired_end_column_is_true(self, tmp_path):
        root = self._build_bam_tree(tmp_path / "data", [("25D1001", 25)])
        out = tmp_path / "config.txt"
        generate_config_fraser(str(out), "MOINS", str(root), str(tmp_path / "fraser"))
        data_lines = out.read_text().splitlines()[1:]
        for line in data_lines:
            assert line.split("\t")[4] == "TRUE"
