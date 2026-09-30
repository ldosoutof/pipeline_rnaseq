#!/usr/bin/env python3
"""
merge_blacklists.py

Fusionne les blacklists FRASER et OUTRIDER (ancien format "sampleID raison",
séparateur espace, sans en-tête) en UN fichier unifié au format attendu par le
pipeline (_load_blacklist) :

    sample_id<TAB>tool<TAB>reason
    23D2296<TAB>fraser<TAB>RNU
    23D2544<TAB>outrider<TAB>gtf2i

- Chaque fichier source est tagué avec son tool (fraser / outrider).
- Un même échantillon peut apparaître dans les deux -> deux lignes (une par tool),
  ce qui est correct (exclu de fraser ET d'outrider).
- Déduplique au sein d'un même tool (garde la première raison).
- Ignore lignes vides ou malformées (ex. "BAB" sans raison), en les signalant.

Usage :
    python merge_blacklists.py \
        --fraser   blacklist_fraser.txt \
        --outrider blacklist_outrider.txt \
        --output   blacklist.tsv
"""
import argparse
import sys


def read_old_format(path, tool):
    """Lit un fichier ancien format -> liste (sample_id, tool, reason), dédupliqué."""
    seen = {}
    order = []
    skipped = []
    dup = 0
    with open(path) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                skipped.append((lineno, line))
                continue
            sid = parts[0]
            reason = "_".join(parts[1:])
            if sid in seen:
                dup += 1
                continue
            seen[sid] = reason
            order.append(sid)
    entries = [(sid, tool, seen[sid]) for sid in order]
    return entries, dup, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fraser",   required=True, help="blacklist_fraser.txt (ancien format)")
    ap.add_argument("--outrider", required=True, help="blacklist_outrider.txt (ancien format)")
    ap.add_argument("--output",   required=True, help="blacklist.tsv unifiée (sortie)")
    args = ap.parse_args()

    all_entries = []
    for path, tool in [(args.fraser, "fraser"), (args.outrider, "outrider")]:
        entries, dup, skipped = read_old_format(path, tool)
        all_entries.extend(entries)
        sys.stderr.write(
            f"[{tool}] {len(entries)} échantillons uniques | {dup} doublons ignorés "
            f"| {len(skipped)} lignes malformées ignorées.\n"
        )
        for lineno, line in skipped:
            sys.stderr.write(f"  [{tool}] ligne {lineno} ignorée : '{line}'\n")

    with open(args.output, "w") as out:
        out.write("sample_id\ttool\treason\n")
        for sid, tool, reason in all_entries:
            out.write(f"{sid}\t{tool}\t{reason}\n")

    sys.stderr.write(f"[merge] {len(all_entries)} lignes écrites dans {args.output}\n")


if __name__ == "__main__":
    main()
