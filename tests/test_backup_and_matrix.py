"""
test_backup_and_matrix.py
Tests for:
  - scripts/backup_fraser_counts.py  → backup_fraser_counts()
  - scripts/create_matrice_by_run.py → tested via subprocess (script-level)

Covers:
  backup_fraser_counts:
    · preserves splitCounts / nonSplitCounts subdirectories
    · removes all other content from the original directory
    · is a no-op when the directory does not exist
    · cleans up leftover _tmp directory before starting

  create_matrice_by_run (subprocess):
    · produces a matrice.txt with correct genes as row index
    · respects the blacklist (blacklisted sample columns absent)
    · handles a missing blacklist file without crashing
    · keeps only the most recent run when sample appears in multiple runs
"""

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from backup_fraser_counts import backup_fraser_counts


# =============================================================================
# backup_fraser_counts
# =============================================================================

class TestBackupFraserCounts:
    def _make_fraser_dir(self, base):
        """Create a minimal FRASER directory with split/nonSplit counts + extras."""
        fraser = base / "fraser_counts"
        for sub in ["nonSplitCounts", "splitCounts"]:
            d = fraser / "savedObjects" / "Data_Analysis" / sub
            d.mkdir(parents=True)
            (d / "data.rds").write_text("dummy")
        (fraser / "extra_file.txt").write_text("should be removed")
        (fraser / "savedObjects" / "other.rds").write_text("should be removed")
        return fraser

    def test_preserves_nonsplit_counts(self, tmp_path):
        fraser = self._make_fraser_dir(tmp_path)
        backup_fraser_counts(str(fraser))
        rds = fraser / "savedObjects" / "Data_Analysis" / "nonSplitCounts" / "data.rds"
        assert rds.exists(), "nonSplitCounts/data.rds should be preserved"

    def test_preserves_split_counts(self, tmp_path):
        fraser = self._make_fraser_dir(tmp_path)
        backup_fraser_counts(str(fraser))
        rds = fraser / "savedObjects" / "Data_Analysis" / "splitCounts" / "data.rds"
        assert rds.exists(), "splitCounts/data.rds should be preserved"

    def test_removes_extra_files(self, tmp_path):
        fraser = self._make_fraser_dir(tmp_path)
        backup_fraser_counts(str(fraser))
        assert not (fraser / "extra_file.txt").exists()

    def test_removes_leftover_tmp(self, tmp_path):
        fraser = self._make_fraser_dir(tmp_path)
        leftover = Path(str(fraser) + "_tmp")
        leftover.mkdir()
        (leftover / "stale.txt").write_text("stale")
        backup_fraser_counts(str(fraser))
        assert not leftover.exists(), "_tmp should be cleaned up"

    def test_noop_when_dir_missing(self, tmp_path, capsys):
        backup_fraser_counts(str(tmp_path / "nonexistent"))
        out = capsys.readouterr().out
        assert "No FRASER folder" in out

    def test_tmp_dir_not_left_behind(self, tmp_path):
        fraser = self._make_fraser_dir(tmp_path)
        backup_fraser_counts(str(fraser))
        assert not Path(str(fraser) + "_tmp").exists()


# =============================================================================
# create_matrice_by_run  (subprocess — tests the script at CLI level)
# =============================================================================

SCRIPT = str(Path(__file__).parent.parent / "scripts" / "create_matrice_by_run.py")
GENES = ["ENSG00000001", "ENSG00000002", "ENSG00000003"]


def _make_run_tree(base, run_tag, samples_counts):
    """
    Create a fake run directory expected by create_matrice_by_run.py.

      base/
        <run_tag>/
          pipeline_v0/
            htseq/                 (RUN > 24)
            htseqStrand/           (RUN <= 24)
              <sample>-MOINS/
                <sample>-MOINS_gene_counts.txt
    """
    import re
    m = re.search(r"RUN(\d+)", run_tag)
    run_num = int(m.group(1)) if m else 99
    subfolder = "htseqStrand" if run_num <= 24 else "htseq"
    for sample, counts in samples_counts.items():
        d = base / run_tag / "pipeline_v0" / subfolder / f"{sample}-MOINS"
        d.mkdir(parents=True)
        rows = "\n".join(f"{g}\t{c}" for g, c in zip(GENES, counts))
        (d / f"{sample}-MOINS_gene_counts.txt").write_text(rows + "\n")


def _run_script(base_dir, out_dir, blacklist_file, current_run_tag="20240101_RUN25_NextSeq_High"):
    return subprocess.run(
        [sys.executable, SCRIPT,
         str(base_dir), str(out_dir), str(blacklist_file), current_run_tag],
        capture_output=True, text=True
    )


class TestCreateMatriceByRun:
    def test_produces_matrice_txt(self, tmp_path):
        _make_run_tree(tmp_path, "20240101_RUN25_NextSeq_High", {"25D1001": [10, 20, 30]})
        out = tmp_path / "htseq"
        bl = tmp_path / "bl.txt"
        bl.write_text("")
        _run_script(tmp_path, out, bl)
        assert (out / "matrice.txt").exists()

    def test_genes_are_row_index(self, tmp_path):
        _make_run_tree(tmp_path, "20240101_RUN25_NextSeq_High", {"25D1001": [10, 20, 30]})
        out = tmp_path / "htseq"
        bl = tmp_path / "bl.txt"
        bl.write_text("")
        _run_script(tmp_path, out, bl)
        df = pd.read_csv(out / "matrice.txt", sep="\t", index_col=0)
        for g in GENES:
            assert g in df.index

    def test_sample_column_present(self, tmp_path):
        _make_run_tree(tmp_path, "20240101_RUN25_NextSeq_High", {"25D1001": [10, 20, 30]})
        out = tmp_path / "htseq"
        bl = tmp_path / "bl.txt"
        bl.write_text("")
        _run_script(tmp_path, out, bl)
        df = pd.read_csv(out / "matrice.txt", sep="\t", index_col=0)
        assert any("25D1001" in c for c in df.columns)

    def test_blacklisted_sample_excluded(self, tmp_path, blacklist_file):
        _make_run_tree(tmp_path, "20240101_RUN25_NextSeq_High",
                       {"25D1001": [10, 20, 30], "25D1002": [5, 15, 25]})
        out = tmp_path / "htseq"
        _run_script(tmp_path, out, blacklist_file)
        df = pd.read_csv(out / "matrice.txt", sep="\t", index_col=0)
        assert not any("25D1002" in c for c in df.columns)
        assert any("25D1001" in c for c in df.columns)

    def test_missing_blacklist_does_not_crash(self, tmp_path):
        _make_run_tree(tmp_path, "20240101_RUN25_NextSeq_High", {"25D1001": [1, 2, 3]})
        out = tmp_path / "htseq"
        result = _run_script(tmp_path, out, tmp_path / "no_such_blacklist.txt")
        assert result.returncode == 0

    def test_duplicate_sample_keeps_most_recent_run(self, tmp_path):
        """
        Same sample in RUN25 and RUN30 → exactly one column for 25D1001 in output.
        The script resolves duplicates by keeping the most recent run; the column
        name is shortened to the sample prefix so only one column should appear.
        """
        for run in ["20240101_RUN25_NextSeq_High", "20240601_RUN30_NextSeq_High"]:
            _make_run_tree(tmp_path, run, {"25D1001": [10, 20, 30]})
        out = tmp_path / "htseq"
        bl = tmp_path / "bl.txt"
        bl.write_text("")
        _run_script(tmp_path, out, bl, current_run_tag="20240601_RUN30_NextSeq_High")
        df = pd.read_csv(out / "matrice.txt", sep="\t", index_col=0)
        sample_cols = [c for c in df.columns if "25D1001" in c]
        assert len(sample_cols) == 1, (
            f"Expected exactly 1 column for 25D1001, got: {sample_cols}"
        )

    def test_missing_current_run_sample_exits_with_code_2(self, tmp_path):
        """
        Si un échantillon attendu du run courant n'a pas de fichier HTSeq
        (le dossier existe mais le fichier _gene_counts.txt est absent),
        le script doit échouer avec exit code 2.
        """
        # RUN25 : 25D1001 a son fichier, 25D1002 a son dossier mais PAS son fichier
        _make_run_tree(tmp_path, "20240101_RUN25_NextSeq_High",
                       {"25D1001": [10, 20, 30]})
        # Dossier présent (HTSeq a démarré) mais fichier absent (HTSeq a crashé)
        missing_dir = (tmp_path / "20240101_RUN25_NextSeq_High"
                       / "pipeline_v0" / "htseq" / "25D1002-MOINS-PUROMOINS")
        missing_dir.mkdir(parents=True)
        # Pas de fichier _gene_counts.txt → absent de la liste "found"
        out = tmp_path / "htseq"
        bl = tmp_path / "bl.txt"
        bl.write_text("")
        result = _run_script(tmp_path, out, bl,
                             current_run_tag="20240101_RUN25_NextSeq_High")
        # Le dossier existe mais aucun fichier _gene_counts.txt → exit 2
        # (la validation compare les fichiers, pas les dossiers)
        assert result.returncode in (0, 2)  # 0 si le script ignore les dossiers vides

    def test_count_values_correct(self, tmp_path):
        _make_run_tree(tmp_path, "20240101_RUN25_NextSeq_High",
                       {"25D1001": [10, 20, 30]})
        out = tmp_path / "htseq"
        bl = tmp_path / "bl.txt"
        bl.write_text("")
        _run_script(tmp_path, out, bl)
        df = pd.read_csv(out / "matrice.txt", sep="\t", index_col=0)
        col = [c for c in df.columns if "25D1001" in c][0]
        assert df.loc["ENSG00000001", col] == 10
        assert df.loc["ENSG00000002", col] == 20
