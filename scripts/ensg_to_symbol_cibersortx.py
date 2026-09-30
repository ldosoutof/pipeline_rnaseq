#!/usr/bin/env python3
"""
ensg_to_symbol_cibersortx.py

Convertit une matrice TPM gène-level indexée par identifiant Ensembl (ENSG…)
en matrice indexée par symbole de gène HGNC, au format attendu par CIBERSORTx
(mixture file) : 1re colonne = symbole, colonnes suivantes = échantillons, TSV.

- Correspondance ENSG -> symbole lue depuis le GTF Ensembl du pipeline
  (mêmes gene_id / gene_name que le reste du pipeline : pas de dépendance externe).
- Agrégation des ENSG multiples partageant un même symbole : SOMME des TPM
  (standard pour la déconvolution ; préserve l'abondance totale du gène).
- Les ENSG sans symbole (gene_name absent) sont ignorés (loggés).

Usage :
    python ensg_to_symbol_cibersortx.py \
        --matrix   matrice_gene_tpm.tsv \
        --gtf      Homo_sapiens.GRCh38.106.gtf \
        --out      mixture_cibersortx.tsv
"""
import argparse
import sys
import re
from collections import defaultdict


def parse_gtf_id_to_symbol(gtf_path):
    """
    Lit le GTF et renvoie un dict {gene_id_sans_version: gene_name}.
    On lit uniquement les lignes 'gene' (suffit pour la correspondance) et on
    retire la version du gene_id (ENSG00000000003.14 -> ENSG00000000003) pour
    matcher une matrice qui peut être versionnée ou non.
    """
    id_re   = re.compile(r'gene_id "([^"]+)"')
    name_re = re.compile(r'gene_name "([^"]+)"')
    mapping = {}
    n_lines = 0
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = cols[8]
            mid = id_re.search(attrs)
            mname = name_re.search(attrs)
            if not mid:
                continue
            gene_id = mid.group(1).split(".")[0]   # retire la version
            gene_name = mname.group(1) if mname else None
            if gene_name:
                mapping[gene_id] = gene_name
            n_lines += 1
    sys.stderr.write(f"[ensg_to_symbol] {len(mapping)} correspondances "
                     f"ENSG->symbole lues depuis le GTF ({n_lines} lignes gene).\n")
    return mapping


def convert(matrix_path, mapping, out_path):
    """
    Lit la matrice TPM (ENSG en 1re colonne), convertit en symboles, somme les
    doublons, écrit la matrice CIBERSORTx.
    """
    with open(matrix_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # La 1re cellule de l'en-tête peut être vide (index R) : on la remplace.
        sample_cols = header[1:]
        n_samples = len(sample_cols)

        # somme des TPM par symbole
        summed = defaultdict(lambda: [0.0] * n_samples)
        n_rows = n_nomap = n_mapped = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            n_rows += 1
            ensg = parts[0].split(".")[0]          # retire version éventuelle
            symbol = mapping.get(ensg)
            if symbol is None:
                n_nomap += 1
                continue
            n_mapped += 1
            vals = summed[symbol]
            for i in range(n_samples):
                try:
                    vals[i] += float(parts[i + 1])
                except (ValueError, IndexError):
                    pass  # cellule vide/non numérique -> ignorée (reste inchangée)

    # écriture au format CIBERSORTx : 'GeneSymbol' + échantillons
    with open(out_path, "w") as out:
        out.write("GeneSymbol\t" + "\t".join(sample_cols) + "\n")
        for symbol in sorted(summed):
            vals = summed[symbol]
            out.write(symbol + "\t" + "\t".join(f"{v:.6g}" for v in vals) + "\n")

    sys.stderr.write(
        f"[ensg_to_symbol] {n_rows} lignes lues | {n_mapped} mappées | "
        f"{n_nomap} sans symbole (ignorées) | {len(summed)} symboles uniques écrits.\n"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, help="matrice TPM (ENSG en 1re colonne)")
    ap.add_argument("--gtf",    required=True, help="GTF Ensembl (gene_id + gene_name)")
    ap.add_argument("--out",    required=True, help="matrice de sortie (mixture CIBERSORTx)")
    args = ap.parse_args()

    mapping = parse_gtf_id_to_symbol(args.gtf)
    if not mapping:
        sys.exit("[ensg_to_symbol] ERREUR : aucune correspondance lue du GTF.")
    convert(args.matrix, mapping, args.out)


if __name__ == "__main__":
    main()
