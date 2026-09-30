#!/usr/bin/env python3
"""
convert_blacklist.py

Convertit une blacklist FRASER de l'ancien format (une ligne "sampleID raison",
séparateur espace, sans en-tête) vers le format TSV attendu par le pipeline :

    sample_id<TAB>tool<TAB>reason
    23D2296<TAB>fraser<TAB>RNU

- Déduplique les échantillons (garde la première raison rencontrée).
- Ignore les lignes vides ou malformées (ex. une ligne isolée "BAB" sans raison).
- Assigne le tool passé en argument (défaut : fraser).

Usage :
    python convert_blacklist.py --input blacklist_fraser.txt \
        --output blacklist.tsv --tool fraser
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  required=True, help="blacklist ancien format")
    ap.add_argument("--output", required=True, help="blacklist TSV de sortie")
    ap.add_argument("--tool",   default="fraser",
                    help="outil concerné (fraser/outrider/all) [défaut: fraser]")
    args = ap.parse_args()

    seen = {}           # sample_id -> reason (première rencontrée)
    order = []          # préserve l'ordre d'apparition
    skipped = []        # lignes ignorées (malformées)
    dup = 0

    with open(args.input) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()          # séparateur = espaces (1 ou plusieurs)
            if len(parts) < 2:
                # ligne sans raison (ex. "BAB" seul) -> ignorée et signalée
                skipped.append((lineno, line))
                continue
            sample_id = parts[0]
            reason = "_".join(parts[1:])  # raison = tout après le sample (jointe)
            if sample_id in seen:
                dup += 1
                continue
            seen[sample_id] = reason
            order.append(sample_id)

    # écriture TSV avec en-tête
    with open(args.output, "w") as out:
        out.write("sample_id\ttool\treason\n")
        for sid in order:
            out.write(f"{sid}\t{args.tool}\t{seen[sid]}\n")

    # rapport
    sys.stderr.write(
        f"[convert_blacklist] {len(order)} échantillons uniques écrits "
        f"(tool={args.tool}) | {dup} doublons ignorés | "
        f"{len(skipped)} lignes malformées ignorées.\n"
    )
    for lineno, line in skipped:
        sys.stderr.write(f"  ligne {lineno} ignorée (pas de raison) : '{line}'\n")


if __name__ == "__main__":
    main()
