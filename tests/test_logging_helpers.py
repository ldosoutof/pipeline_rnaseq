"""
test_logging_helpers.py
Tests for the Python functions defined in rules/logging.smk:
  - log_start()           : generates a valid bash snippet
  - collect_failed_logs() : scans log/ for failed rules

These are pure-Python helpers; no Snakemake execution required.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Import logging.smk as a Python module
# (it is valid Python even though it lives in rules/ with a .smk extension)
# ---------------------------------------------------------------------------

def _load_logging_smk():
    import importlib.machinery
    smk_path = Path(__file__).parent.parent / "rules" / "logging.smk"
    loader = importlib.machinery.SourceFileLoader("logging_smk", str(smk_path))
    spec = importlib.util.spec_from_loader("logging_smk", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


logging_smk = _load_logging_smk()
log_start = logging_smk.log_start
collect_failed_logs = logging_smk.collect_failed_logs


# =============================================================================
# log_start
# =============================================================================

class TestLogStart:
    def _wildcards(self, **kwargs):
        """Simple namespace that behaves like Snakemake wildcards."""
        class _WC(dict):
            pass
        return _WC(kwargs)

    def test_returns_string(self):
        result = log_start("fastqc_report", self._wildcards(sample="25D1001"), threads=8)
        assert isinstance(result, str)

    def test_contains_rule_name(self):
        result = log_start("alignment_star", self._wildcards(sample="25D1001"), threads=16)
        assert "alignment_star" in result

    def test_contains_mkdir(self):
        """Shell snippet must create the log directory before writing."""
        result = log_start("fastp", self._wildcards(sample="25D1001"))
        assert "mkdir" in result

    def test_contains_echo_or_printf(self):
        """Shell snippet must write something to the log file."""
        result = log_start("htseq_gene", self._wildcards(sample="25D1001"))
        assert "echo" in result or "printf" in result

    def test_wildcard_values_appear_in_snippet(self):
        """Sample wildcard value should be embedded so logs are self-contained."""
        result = log_start("outrider", self._wildcards(sample="25D9999"))
        assert "25D9999" in result

    def test_none_wildcards_dont_crash(self):
        result = log_start("fraser", None)
        assert isinstance(result, str)

    def test_threads_embedded(self):
        result = log_start("alignment_star", self._wildcards(sample="X"), threads=32)
        assert "32" in result


# =============================================================================
# collect_failed_logs
# =============================================================================

class TestCollectFailedLogs:
    def _write_log(self, log_dir, rule_name, exit_code):
        """Write a synthetic structured log (start + end JSON lines)."""
        d = log_dir / rule_name
        d.mkdir(parents=True, exist_ok=True)
        p = d / "log.txt"
        lines = [
            json.dumps({"event": "start", "rule": rule_name,
                        "wildcards": {}, "threads": 1,
                        "timestamp": "2026-01-01T00:00:00Z"}),
            "some tool output line",
            json.dumps({"event": "end", "rule": rule_name,
                        "exit_code": exit_code,
                        "timestamp": "2026-01-01T00:01:00Z"}),
        ]
        p.write_text("\n".join(lines))
        return p

    def test_returns_empty_list_when_all_succeeded(self, tmp_path):
        log_dir = tmp_path / "log"
        self._write_log(log_dir, "fastqc_report", exit_code=0)
        self._write_log(log_dir, "fastp", exit_code=0)
        result = collect_failed_logs(log_dir=str(log_dir))
        assert result == []

    def test_detects_single_failure(self, tmp_path):
        log_dir = tmp_path / "log"
        self._write_log(log_dir, "alignment_star", exit_code=1)
        result = collect_failed_logs(log_dir=str(log_dir))
        assert len(result) == 1
        assert result[0]["exit_code"] == 1

    def test_detects_multiple_failures(self, tmp_path):
        log_dir = tmp_path / "log"
        self._write_log(log_dir, "outrider", exit_code=1)
        self._write_log(log_dir, "fraser", exit_code=2)
        self._write_log(log_dir, "fastqc_report", exit_code=0)
        result = collect_failed_logs(log_dir=str(log_dir))
        assert len(result) == 2

    def test_failure_record_has_expected_keys(self, tmp_path):
        log_dir = tmp_path / "log"
        self._write_log(log_dir, "htseq_gene", exit_code=1)
        result = collect_failed_logs(log_dir=str(log_dir))
        rec = result[0]
        assert "rule" in rec
        assert "log_path" in rec
        assert "exit_code" in rec
        assert "tail" in rec

    def test_tail_contains_log_content(self, tmp_path):
        log_dir = tmp_path / "log"
        self._write_log(log_dir, "matrix", exit_code=127)
        result = collect_failed_logs(log_dir=str(log_dir))
        assert "some tool output line" in result[0]["tail"]

    def test_returns_empty_when_log_dir_missing(self, tmp_path):
        result = collect_failed_logs(log_dir=str(tmp_path / "no_such_dir"))
        assert result == []

    def test_ignores_log_with_no_json_end_record(self, tmp_path):
        """A log file with only plain text (tool crashed before LOG_END) is skipped."""
        log_dir = tmp_path / "log"
        d = log_dir / "bad_rule"
        d.mkdir(parents=True)
        (d / "log.txt").write_text("Segfault\nAborted (core dumped)\n")
        result = collect_failed_logs(log_dir=str(log_dir))
        assert result == []

    def test_exit_code_zero_is_not_failure(self, tmp_path):
        log_dir = tmp_path / "log"
        self._write_log(log_dir, "multiqc", exit_code=0)
        result = collect_failed_logs(log_dir=str(log_dir))
        assert result == []

    def test_non_zero_exit_codes_all_captured(self, tmp_path):
        log_dir = tmp_path / "log"
        for code in [1, 2, 127, 139]:
            self._write_log(log_dir, f"rule_exit_{code}", exit_code=code)
        result = collect_failed_logs(log_dir=str(log_dir))
        codes = {r["exit_code"] for r in result}
        assert codes == {1, 2, 127, 139}
