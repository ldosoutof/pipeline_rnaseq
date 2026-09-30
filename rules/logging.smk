# =============================================================================
# logging.smk  —  shared structured logging helpers
# Include ONCE in pipeline.smk before any rule file:
#   include: '../rules/logging.smk'
#
# Every rule shell block should start / end with the two macros:
#
#   {params.log_start}   — writes a JSON header to {log.run_info}
#   {params.log_end}     — appends a JSON footer (exit_code, duration, outputs)
#
# On failure Snakemake calls onerror; the last log entry written is the
# JSON footer with exit_code != 0, making grep / jq trivial.
# =============================================================================

import json
import os
from pathlib import Path
from datetime import datetime


def _log_header(rule, wildcards_dict, input_dict, params_dict, threads, log_path):
    """
    Return a bash snippet that writes a JSON header line to log_path.
    Called at rule parse-time to build the params.log_start string.
    The actual timestamp is evaluated at shell execution time via $(date …).
    """
    meta = {
        "event":     "start",
        "rule":      rule,
        "wildcards": wildcards_dict,
        "threads":   threads,
    }
    # Embed static metadata; timestamp injected at runtime
    meta_json = json.dumps(meta, separators=(',', ':'))
    # Strip the closing } so we can splice in the runtime timestamp
    meta_open = meta_json[:-1]  # everything except trailing }
    snippet = (
        'mkdir -p "$(dirname {log.run_info})" && '
        f'echo \'{meta_open},"timestamp":"\'$(date -u +"%Y-%m-%dT%H:%M:%SZ")\'"\' >> {{log.run_info}}'
    )
    return snippet


def _log_footer():
    """
    Return a bash snippet appended after the tool command.
    Captures $? so the exit code is recorded even when the rule is about
    to fail — Snakemake sees the non-zero exit AFTER this snippet runs.
    """
    return (
        '_EXIT=$? ; '
        '_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ") ; '
        'echo \'{"event":"end","exit_code":\'$_EXIT\',"timestamp":"\'$_END\'"}\' >> {log.run_info} ; '
        '(exit $_EXIT)'          # re-raise so Snakemake marks the rule failed
    )


# ---------------------------------------------------------------------------
# Shell preamble / postamble strings — injected via params in every rule
# ---------------------------------------------------------------------------
# Usage in a rule:
#
#   params:
#       ...,
#       log_start = lambda wc, input, threads: log_start(RULENAME, wc, threads),
#       log_end   = LOG_END,
#
#   shell:
#       """
#       set -euo pipefail
#       {params.log_start}
#       <tool command> >> {log.run_info} 2>&1
#       {params.log_end}
#       """
#
# LOG_END is a module-level constant; log_start() is a helper function.
# ---------------------------------------------------------------------------

LOG_END = (
    # Prefer $SNAKEMAKE_LOG (set by log_start); fall back to {log.run_info}
    # so LOG_END is safe even if log_start was not called first.
    # Second fallback to /dev/null guards against rules with no log: block.
    ': "${SNAKEMAKE_LOG:={log.run_info}}" ; '
    ': "${SNAKEMAKE_LOG:=/dev/null}" ; '
    '_EXIT=$? ; '
    '_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ") ; '
    'echo \'{"event":"end","exit_code":\'$_EXIT\',"timestamp":"\'$_END\'"}\' >> "$SNAKEMAKE_LOG" ; '
    '(exit $_EXIT)'
)


def log_start(rule_name, wildcards, threads=1):
    """
    Returns a bash one-liner that appends a JSON 'start' record to {log.run_info}.
    Wildcard values are serialised so the log is self-contained.

    Usage inside a rule params lambda:
        log_start = lambda wc, input, threads: log_start("my_rule", wc, threads),
    """
    wc_dict = dict(wildcards) if wildcards else {}
    meta = json.dumps({
        "event":     "start",
        "rule":      rule_name,
        "wildcards": wc_dict,
        "threads":   threads,
    }, separators=(',', ':'))
    # Drop the trailing } so we can splice in the runtime timestamp.
    # We use plain string concatenation to avoid f-string / .format() collisions
    # with the Snakemake {log.run_info} placeholder.
    meta_open = meta[:-1]   # e.g. '{"event":"start","rule":"fastqc",...'
    ts_suffix  = ',"timestamp":"\'$(date -u +"%Y-%m-%dT%H:%M:%SZ")\'"}'
    log_ref    = "$SNAKEMAKE_LOG"
    return (
        'SNAKEMAKE_LOG="{log.run_info}" ; '
        + 'mkdir -p "$(dirname $SNAKEMAKE_LOG)" && '
        + "echo '" + meta_open + ts_suffix + "' >> " + log_ref
    )


# ---------------------------------------------------------------------------
# onerror hook helper — call from pipeline.smk onerror block
# ---------------------------------------------------------------------------

def collect_failed_logs(log_dir="log", tail_lines=30):
    """
    Scan log_dir for run_info logs whose last JSON line has exit_code != 0.
    Returns a list of dicts  { rule, log_path, last_lines, exit_code }.
    Used in the onerror email to give an instant summary of what broke.
    """
    failures = []
    log_root = Path(log_dir)
    if not log_root.exists():
        return failures

    for log_file in sorted(log_root.rglob("log.txt")) :
        try:
            lines = log_file.read_text(errors="replace").splitlines()
        except OSError:
            continue

        # Walk backwards looking for a JSON end record
        for line in reversed(lines):
            line = line.strip()
            if not line.startswith("{"):
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("event") == "end" and rec.get("exit_code", 0) != 0:
                failures.append({
                    "rule":       rec.get("rule", str(log_file.parent.name)),
                    "log_path":   str(log_file),
                    "exit_code":  rec["exit_code"],
                    "timestamp":  rec.get("timestamp", "?"),
                    "tail":       "\n".join(lines[-tail_lines:]),
                })
            break   # only inspect the last end record per file

    return failures


# ---------------------------------------------------------------------------
# Lint helper — call from CI or a pre-run check to detect rules that have
# one of log_start / log_end but not both.
# ---------------------------------------------------------------------------

def check_logging_symmetry(smk_dir="rules"):
    """
    Parse all .smk files under smk_dir and return a list of (file, rule_name)
    tuples where log_start is present without log_end, or vice-versa.

    Usage (from pipeline.smk or a standalone script):
        from rules.logging import check_logging_symmetry
        issues = check_logging_symmetry("rules")
        if issues:
            raise WorkflowError("Logging asymmetry: " + str(issues))
    """
    import re
    from pathlib import Path

    issues = []
    for smk in sorted(Path(smk_dir).rglob("*.smk")):
        content = smk.read_text(errors="replace")
        # Find each rule block (naive split — good enough for lint)
        for m in re.finditer(r'^rule\s+(\w+)\s*:', content, re.MULTILINE):
            rule_name = m.group(1)
            # Extract the rule body up to the next rule or EOF
            start = m.start()
            next_rule = re.search(r'^rule\s+\w+\s*:', content[start + 1:], re.MULTILINE)
            end = start + 1 + next_rule.start() if next_rule else len(content)
            body = content[start:end]
            has_start = bool(re.search(r'log_start\s*=\s*lambda', body))
            has_end   = bool(re.search(r'log_end\s*=\s*LOG_END', body))
            if has_start != has_end:
                issues.append((str(smk), rule_name, "log_start" if has_start else "log_end"))
    return issues
