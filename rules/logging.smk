# =============================================================================
# logging.smk  —  shared structured logging helpers
# Include ONCE in pipeline.smk before any rule file:
#   include: '../rules/logging.smk'
#
# Every rule shell block starts with {params.log_start} and ends with
# {params.log_end}.
#
#   {params.log_start}  appends a JSON "start" record to log/events/<rule>.jsonl
#                       and installs an EXIT trap that appends the JSON "end"
#                       record with the REAL exit code — on success AND on failure.
#   {params.log_end}    kept for compatibility; no-op (the trap writes the end).
#
# The events file path is computed in Python (rule name), so it never depends on
# a Snakemake placeholder: a params value is inserted verbatim into the shell
# command and is NOT formatted a second time (a literal "{log.run_info}" in a
# params value used to create a file literally named "{log.run_info}").
# =============================================================================

import json
import os
from pathlib import Path
from datetime import datetime


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

LOG_END = ":"   # no-op : the EXIT trap installed by log_start writes the end record


def log_start(rule_name, wildcards, threads=1):
    """
    Bash snippet: append a JSON 'start' record to log/events/<rule>.jsonl and
    install an EXIT trap appending the JSON 'end' record with the real exit code
    (the trap also fires when `set -e` aborts the shell on a failing command).

    Usage inside a rule params lambda:
        log_start = lambda wc, input, threads: log_start("my_rule", wc, threads),
    """
    wc_dict = dict(wildcards) if wildcards else {}
    meta = json.dumps({"rule": rule_name, "wildcards": wc_dict, "threads": threads},
                      separators=(',', ':'))[1:-1]          # inner part, no braces
    meta = meta.replace("'", "'\\''")                         # safe inside '...'
    log_file = f"log/events/{rule_name}.jsonl"
    ts = '$(date -u +%Y-%m-%dT%H:%M:%SZ)'
    start = ('echo "{\\"event\\":\\"start\\",$_LOG_META,\\"timestamp\\":\\"' + ts
             + '\\"}" >> "$_LOG_FILE"')
    end = ('_rc=$? ; echo "{\\"event\\":\\"end\\",$_LOG_META,\\"exit_code\\":$_rc,'
           '\\"timestamp\\":\\"' + ts + '\\"}" >> "$_LOG_FILE"')
    return (f"_LOG_FILE='{log_file}' ; _LOG_META='{meta}' ; "
            'mkdir -p "$(dirname "$_LOG_FILE")" ; '
            + start + " ; "
            + f"trap '{end}' EXIT")


# ---------------------------------------------------------------------------
# onerror hook helper — call from pipeline.smk onerror block
# ---------------------------------------------------------------------------

def collect_failed_logs(log_dir="log", tail_lines=30):
    """
    Read log/events/<rule>.jsonl and return one dict per failed or interrupted
    job: { rule, log_path, exit_code, timestamp, tail }.
    A job is identified by (rule, wildcards). exit_code is None when a start
    record has no matching end record (job killed or still running).
    """
    failures = []
    ev_root = Path(log_dir) / "events"
    if not ev_root.exists():
        return failures
    for ev_file in sorted(ev_root.glob("*.jsonl")):
        try:
            lines = ev_file.read_text(errors="replace").splitlines()
        except OSError:
            continue
        last, raw = {}, {}             # job key -> last record / its raw lines
        for line in lines:
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            key = json.dumps(rec.get("wildcards", {}), sort_keys=True)
            last[key] = rec
            raw.setdefault(key, []).append(line)
        for key, rec in last.items():
            if rec.get("event") == "start" or rec.get("exit_code", 0) != 0:
                rule = rec.get("rule", ev_file.stem)
                text_log = _find_rule_log(Path(log_dir), rule, rec.get("wildcards", {}))
                tail = raw[key][-tail_lines:]
                if text_log is not None:
                    try:
                        tail += ["--- " + str(text_log) + " ---"] + \
                                text_log.read_text(errors="replace").splitlines()[-tail_lines:]
                    except OSError:
                        pass
                failures.append({
                    "rule":      rule,
                    "log_path":  str(text_log or ev_file),
                    "exit_code": rec.get("exit_code"),
                    "timestamp": rec.get("timestamp", "?"),
                    "tail":      "\n".join(tail),
                })
    return failures


def _find_rule_log(log_root, rule, wildcards):
    """
    Best effort : log.txt de la règle (sortie de l'outil). Les règles écrivent
    en général dans log/<règle>/[<wildcard>/]log.txt ; on retient le fichier dont
    le chemin contient le nom de la règle et toutes les valeurs de wildcards.
    """
    values = [str(v) for v in (wildcards or {}).values()]
    best = None
    for p in log_root.rglob("log.txt"):
        parts = p.parts
        if rule in parts and all(v in str(p) for v in values):
            if best is None or len(parts) > len(best.parts):
                best = p
    return best


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
