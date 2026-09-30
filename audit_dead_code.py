#!/usr/bin/env python3
"""
audit_dead_code.py — inventaire du code mort dans le pipeline RNA-seq
======================================================================

Lit pipeline.smk, calcule la fermeture des includes, et signale :
  - les .smk non atteignables (candidats archive/)
  - les scripts non appelés par aucune règle active
  - les fichiers non-.smk dans rules/ (archives, .zip, etc.)

Usage :
    python audit_dead_code.py [--pipeline snakemake/pipeline.smk]
                              [--rules-dir rules/]
                              [--scripts-dir scripts/]
                              [--archive-dir rules/archive/]
                              [--move]   # déplace les fichiers au lieu d'afficher
"""

import argparse
import re
import shutil
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Audit dead code in Snakemake pipeline")
    p.add_argument("--pipeline",    default="snakemake/pipeline.smk")
    p.add_argument("--rules-dir",   default="rules")
    p.add_argument("--scripts-dir", default="scripts")
    p.add_argument("--archive-dir", default="rules/archive")
    p.add_argument("--move",        action="store_true",
                   help="Move archive-candidate files to --archive-dir")
    return p.parse_args()


def resolve_includes(pipeline_smk: Path, rules_dir: Path) -> set:
    """Walk include: directives recursively and return the set of included filenames."""
    included = set()
    to_visit = [pipeline_smk]
    while to_visit:
        current = to_visit.pop()
        if not current.exists():
            continue
        content = current.read_text(errors="replace")
        for m in re.finditer(r"include:\s*['\"]([^'\"]+)['\"]", content):
            raw = m.group(1)
            # Resolve relative to rules_dir
            resolved = (rules_dir / Path(raw).name).resolve()
            fname = Path(raw).name
            if fname not in included:
                included.add(fname)
                candidate = rules_dir / fname
                if candidate.exists():
                    to_visit.append(candidate)
                # Also check experimental/
                exp = rules_dir / "experimental" / fname
                if exp.exists():
                    to_visit.append(exp)
    return included


def find_called_scripts(included_smks: list, pipeline_smk: Path) -> set:
    """Return the set of script basenames called by any active rule file."""
    called = set()
    for p in included_smks + [pipeline_smk]:
        if not Path(p).exists():
            continue
        content = Path(p).read_text(errors="replace")
        for m in re.findall(r"scripts/([A-Za-z0-9_]+\.(?:py|R))", content):
            called.add(m)
    return called


def main():
    args = parse_args()
    rules_dir   = Path(args.rules_dir).resolve()
    scripts_dir = Path(args.scripts_dir).resolve()
    pipeline    = Path(args.pipeline).resolve()
    archive_dir = Path(args.archive_dir).resolve()

    print(f"Pipeline  : {pipeline}")
    print(f"Rules dir : {rules_dir}")
    print(f"Scripts   : {scripts_dir}")
    print()

    # ── 1. Included smk files ──────────────────────────────────────────────
    included_names = resolve_includes(pipeline, rules_dir)
    included_names.add("pipeline.smk")   # pipeline itself

    # Collect all .smk on disk (direct children only — experimental/ handled separately)
    all_smks = {f.name: f for f in rules_dir.glob("*.smk")}

    active_paths  = [rules_dir / n for n in included_names if (rules_dir / n).exists()]
    archive_smks  = [f for name, f in sorted(all_smks.items())
                     if name not in included_names]

    print(f"Active rules ({len(included_names)}):")
    for n in sorted(included_names):
        status = "✓" if (rules_dir / n).exists() else "✗ NOT ON DISK"
        print(f"  {status}  {n}")

    print(f"\nArchive candidates ({len(archive_smks)}) — not reachable from pipeline.smk:")
    for f in archive_smks:
        print(f"  {f.name}")

    # ── 2. Non-smk files in rules/ ─────────────────────────────────────────
    odd_files = [f for f in rules_dir.iterdir()
                 if not f.is_dir() and f.suffix != ".smk"
                 and not re.search(r'\.bak\d*$', f.name)]
    if odd_files:
        print(f"\nNon-.smk files in {rules_dir.name}/ ({len(odd_files)}):")
        for f in sorted(odd_files):
            print(f"  {f.name}  ({f.stat().st_size:,} bytes)")

    # ── 3. Orphaned scripts ────────────────────────────────────────────────
    called = find_called_scripts(active_paths, pipeline)
    all_scripts = {f.name for f in scripts_dir.iterdir()
                   if f.suffix in (".py", ".R")}
    orphaned = sorted(all_scripts - called)

    # Known intentional non-rule scripts (CLI tools, ops scripts)
    cli_tools = {
        "backfill_db.py", "backfill_qc_summary.py", "preprocess.py",
        "annotation_coverage.py",   # added by F24, called manually
    }
    # Called by experimental module (not included by default)
    experimental = {"mean_chrY_expression.R", "plot_vaf_violin_by_sample.py"}

    true_orphans = [s for s in orphaned
                    if s not in cli_tools and s not in experimental]

    if experimental & set(orphaned):
        print(f"\nScripts for experimental module (not in active rules by design):")
        for s in sorted(experimental & set(orphaned)):
            print(f"  {s}")

    if cli_tools & set(orphaned):
        print(f"\nCLI/ops scripts (called manually, not from rules):")
        for s in sorted(cli_tools & set(orphaned)):
            print(f"  {s}")

    if true_orphans:
        print(f"\nUnexplained orphaned scripts ({len(true_orphans)}):")
        for s in true_orphans:
            print(f"  {s}")
    else:
        print(f"\n✓ No unexplained orphaned scripts")

    # ── 4. Move if requested ───────────────────────────────────────────────
    if args.move and archive_smks:
        archive_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nMoving {len(archive_smks)} files to {archive_dir} …")
        for f in archive_smks:
            dest = archive_dir / f.name
            shutil.move(str(f), str(dest))
            print(f"  mv {f.name} → archive/")
        print("Done. Use 'git mv' in the actual repo to preserve history.")

    print()
    return 1 if true_orphans else 0


if __name__ == "__main__":
    sys.exit(main())
