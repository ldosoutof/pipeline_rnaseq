#!/usr/bin/env python3
"""
push_metrics.py
===============
Collect RNA-seq pipeline QC metrics from all runs under PROD_ROOT
and push them to a Prometheus Pushgateway.

Called automatically by the pipeline after each run (rule generate_metrics),
or manually for historical backfill:

    python push_metrics.py \\
        --prod_root /datawork2/genetique/RNASeq/diag/prod \\
        --gateway   http://localhost:9091 \\
        --backfill                         # scan all existing runs

Metrics pushed
--------------
Per-sample QC (from qc_summary.tsv):
  rnaseq_nb_reads              Total read pairs
  rnaseq_mapped_reads          Mapped read pairs
  rnaseq_dup_pct               Duplication percentage
  rnaseq_on_target             On-target read count
  rnaseq_nb_outrider_events    OUTRIDER aberration count
  rnaseq_nb_fraser_events      FRASER aberration count
  rnaseq_hba_total_tpm         HBA1+HBA2+HBB summed TPM
  rnaseq_hba1_tpm / hba2 / hbb Individual haemoglobin TPM
  rnaseq_di_green_pct          % DI-green genes with TPM > 10
  rnaseq_di_green_expressed    Count of expressed DI-green genes

Per-rule performance (from benchmarks/*.tsv):
  rnaseq_rule_runtime_seconds  Wall-clock time per rule
  rnaseq_rule_cpu_seconds      CPU time per rule
  rnaseq_rule_max_rss_mb       Peak RSS memory per rule

Pipeline status (from structured logs):
  rnaseq_rule_exit_code        Last exit code per rule
"""

import argparse
import glob
import json
import os
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests


# ---------------------------------------------------------------------------
# Pushgateway helpers
# ---------------------------------------------------------------------------

def _push(gateway, job, labels, metrics, run_name=""):
    """
    Push metrics to Pushgateway.
    Uses run_name as part of the job path so each run has its own slot
    and never overwrites another run's data.
    """
    safe_run  = re.sub(r"[^a-zA-Z0-9_-]", "_", run_name) if run_name else ""
    full_job  = f"{job}/{safe_run}" if safe_run else job
    label_path = "/".join(f"{k}/{v}" for k, v in sorted(labels.items()))
    url = f"{gateway}/metrics/job/{full_job}/{label_path}"

    label_str = ",".join(f'{k}="{v}"' for k, v in sorted(labels.items()))
    lines = []
    for name, (value, help_text, mtype) in metrics.items():
        if value is None:
            continue
        lines.append(f"# HELP {name} {help_text}")
        lines.append(f"# TYPE {name} {mtype}")
        lines.append(f"{name}{{{label_str}}} {float(value)}")

    body = "\n".join(lines) + "\n"
    try:
        r = requests.put(url, data=body.encode(),
                         headers={"Content-Type": "text/plain"},
                         timeout=10)
        r.raise_for_status()
    except Exception as e:
        print(f"[WARN] push failed for {labels}: {e}", file=sys.stderr)


# ---------------------------------------------------------------------------
# QC summary metrics
# ---------------------------------------------------------------------------

def push_qc_summary(run_path, run_name, run_id, gateway):
    """Parse qc_summary.tsv and push per-sample metrics."""
    qc_file = Path(run_path) / "pipeline_v0" / "metrics" / "qc_summary.tsv"
    if not qc_file.exists():
        print(f"[SKIP] No qc_summary.tsv in {run_name}", file=sys.stderr)
        return

    try:
        df = pd.read_csv(qc_file, sep="\t")
    except Exception as e:
        print(f"[WARN] Could not read {qc_file}: {e}", file=sys.stderr)
        return

    for _, row in df.iterrows():
        sample = str(row.get("sample id", "unknown"))
        labels = {"run": run_id, "run_name": run_name, "sample": sample}

        metrics = {
            "rnaseq_nb_reads": (
                row.get("nb read"), "Total read pairs", "gauge"),
            "rnaseq_mapped_reads": (
                row.get("mapped"), "Mapped read pairs", "gauge"),
            "rnaseq_dup_pct": (
                row.get("dup"), "Duplication percentage", "gauge"),
            "rnaseq_on_target": (
                row.get("on target"), "On-target read count", "gauge"),
            "rnaseq_nb_outrider_events": (
                row.get("nb event out"), "OUTRIDER aberration event count", "gauge"),
            "rnaseq_nb_fraser_events": (
                row.get("nb event fraser"), "FRASER aberration event count", "gauge"),
            "rnaseq_hba1_tpm": (
                row.get("HBA1"), "HBA1 TPM expression", "gauge"),
            "rnaseq_hba2_tpm": (
                row.get("HBA2"), "HBA2 TPM expression", "gauge"),
            "rnaseq_hbb_tpm": (
                row.get("HBB"), "HBB TPM expression", "gauge"),
            "rnaseq_hba_total_tpm": (
                row.get("HBA_total"), "HBA1+HBA2+HBB summed TPM", "gauge"),
            "rnaseq_di_green_expressed": (
                row.get("nb_DI_green_expressed"), "DI-green genes expressed (TPM>10)", "gauge"),
            "rnaseq_di_green_pct": (
                row.get("pct_DI_green_TPM>10"), "Percent DI-green genes with TPM>10", "gauge"),
        }
        _push(gateway, "rnaseq_pipeline", labels, metrics, run_name)

    print(f"[OK] QC metrics pushed for {run_name} ({len(df)} samples)")


# ---------------------------------------------------------------------------
# Benchmark metrics
# ---------------------------------------------------------------------------

def push_benchmarks(run_path, run_name, run_id, gateway):
    """
    Parse Snakemake benchmark TSV files and push per-rule performance metrics.
    Benchmark columns: s  h:m:s  max_rss  max_vms  max_uss  max_pss
                       io_in  io_out  mean_load  cpu_time
    """
    bench_dir = Path(run_path) / "benchmarks"
    if not bench_dir.exists():
        return

    for tsv in bench_dir.rglob("*.tsv"):
        # derive rule name from path: benchmarks/<rule>/<sample>.tsv
        parts = tsv.relative_to(bench_dir).parts
        rule = parts[0] if len(parts) >= 1 else tsv.stem
        sample = tsv.stem if len(parts) >= 2 else "all"

        try:
            df = pd.read_csv(tsv, sep="\t")
            if df.empty:
                continue
            row = df.iloc[-1]   # last attempt
        except Exception:
            continue

        labels = {"run": run_id, "run_name": run_name,
                  "rule": rule, "sample": sample}
        metrics = {
            "rnaseq_rule_runtime_seconds": (
                row.get("s"), "Rule wall-clock runtime in seconds", "gauge"),
            "rnaseq_rule_cpu_seconds": (
                row.get("cpu_time"), "Rule CPU time in seconds", "gauge"),
            "rnaseq_rule_max_rss_mb": (
                row.get("max_rss"), "Rule peak RSS memory in MB", "gauge"),
        }
        _push(gateway, "rnaseq_benchmarks", labels, metrics, run_name)

    print(f"[OK] Benchmark metrics pushed for {run_name}")


# ---------------------------------------------------------------------------
# Structured log metrics (exit codes + durations)
# ---------------------------------------------------------------------------

def _parse_log_duration(log_file):
    """Extract start/end timestamps and exit_code from a structured log file."""
    start_ts = end_ts = exit_code = rule = None
    try:
        for line in Path(log_file).read_text(errors="replace").splitlines():
            line = line.strip()
            if not line.startswith("{"):
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("event") == "start":
                start_ts = rec.get("timestamp")
                rule = rec.get("rule", rule)
            elif rec.get("event") == "end":
                end_ts = rec.get("timestamp")
                exit_code = rec.get("exit_code", 0)
    except OSError:
        pass
    return rule, start_ts, end_ts, exit_code


def _iso_to_epoch(ts):
    if not ts:
        return None
    from datetime import datetime, timezone
    try:
        return datetime.fromisoformat(ts.replace("Z", "+00:00")).timestamp()
    except Exception:
        return None


def push_log_metrics(run_path, run_name, run_id, gateway):
    """Push exit codes and durations derived from structured JSON logs."""
    log_dir = Path(run_path) / "log"
    if not log_dir.exists():
        return

    for log_file in sorted(log_dir.rglob("log.txt")):
        rule, start_ts, end_ts, exit_code = _parse_log_duration(log_file)
        if exit_code is None:
            continue

        # derive sample from path: log/<rule>/<sample>/log.txt
        parts = log_file.relative_to(log_dir).parts
        rule_name = parts[0] if len(parts) >= 1 else "unknown"
        sample = parts[1] if len(parts) >= 3 else "all"

        duration = None
        t0, t1 = _iso_to_epoch(start_ts), _iso_to_epoch(end_ts)
        if t0 and t1:
            duration = t1 - t0

        labels = {"run": run_id, "run_name": run_name,
                  "rule": rule_name, "sample": sample}
        metrics = {
            "rnaseq_rule_exit_code": (
                exit_code, "Last exit code of rule execution", "gauge"),
            "rnaseq_rule_duration_seconds": (
                duration, "Rule wall-clock duration from structured logs", "gauge"),
        }
        _push(gateway, "rnaseq_logs", labels, metrics, run_name)

    print(f"[OK] Log metrics pushed for {run_name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_run(run_path, gateway):
    run_name = Path(run_path).name
    m = re.search(r"RUN(\d+)", run_name, re.IGNORECASE)
    run_id = m.group(0) if m else run_name

    push_qc_summary(run_path, run_name, run_id, gateway)
    push_benchmarks(run_path, run_name, run_id, gateway)
    push_log_metrics(run_path, run_name, run_id, gateway)


def main():
    p = argparse.ArgumentParser(description="Push pipeline metrics to Prometheus Pushgateway")
    p.add_argument("--run_path",  help="Path to a single run folder")
    p.add_argument("--prod_root", help="Root production directory (for --backfill)")
    p.add_argument("--gateway",   default="http://localhost:9091",
                   help="Prometheus Pushgateway URL (default: http://localhost:9091)")
    p.add_argument("--backfill",  action="store_true",
                   help="Scan prod_root and push metrics for all existing runs")
    args = p.parse_args()

    if args.backfill:
        if not args.prod_root:
            print("--prod_root required with --backfill", file=sys.stderr)
            sys.exit(1)
        runs = sorted(glob.glob(os.path.join(args.prod_root, "20*_RUN*")))
        print(f"[INFO] Backfilling {len(runs)} runs from {args.prod_root}")
        for run in runs:
            process_run(run, args.gateway)
    elif args.run_path:
        process_run(args.run_path, args.gateway)
    else:
        print("Provide --run_path or --backfill + --prod_root", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
