#!/usr/bin/env python3
"""
validate_known_events.py — retrouve-t-on les événements connus ?

Contrôle de performance (sensibilité sur des positifs connus) et de
non-régression entre deux versions du pipeline. Pour chaque événement attendu
(échantillon, outil, gène), le script cherche le résultat dans les tables
COMPLÈTES d'un run (et non dans les livrables filtrés), ce qui distingue
« non retrouvé » (testé, non significatif) de « non testé » (filtré avant le
test) et d'« exclu » (échantillon absent ou blacklisté).

Usage :
  python3 tools/validate_known_events.py --truth evenements_attendus.tsv \\
      --run /datawork2/.../prod/<run>/pipeline_v0 [--branch both] --out rapport_validation
  # non-régression : référence (production) contre candidat (dev)
  python3 tools/validate_known_events.py --truth ... --run <candidat>/pipeline_v0 \\
      --compare <reference>/pipeline_v0 --out ... --strict

Fichier des événements attendus (TSV, en-tête ; noms de colonnes insensibles à
la casse) : Sample, Tool (OUTRIDER|FRASER), Gene, variant, effet
(« Effet total » / « Effet partiel »). Une ligne sans Sample est rattachée à la
précédente (complément de variant, ex. coordonnée GRCh38).
CE FICHIER CONTIENT DES DONNÉES DE SANTÉ : ne pas le versionner.

Critères (modifiables) :
  OUTRIDER retrouvé : padj < --padj et |log2FC| >= --l2fc  (critère « aberrant » du pipeline)
  FRASER   retrouvé : padj < --padj et |deltaPsi| >= --dpsi
  --max-rank N      : exige en plus un rang <= N dans l'échantillon
  FRASER : jonction cherchée par symbole de gène ; à défaut, à ± --window pb de
  la position GRCh38 du variant (résultat marqué « par position »).
Statuts : RETROUVE, SIGNIFICATIF_EFFET_FAIBLE, NOMINAL_SEULEMENT, NON_RETROUVE,
          NON_TESTE, ABSENT (échantillon absent des résultats), BLACKLISTE.
Code de sortie avec --strict : 1 si un événement « Effet total » n'est pas
RETROUVE, ou (avec --compare) si un événement retrouvé par la référence ne
l'est plus par le candidat.
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

CORE_RE = r"(\d{2,}[A-Z]\d+)"          # identifiant court : 24D1112, X24D1112… -> 24D1112
RANK_ORDER = ["RETROUVE", "SIGNIFICATIF_EFFET_FAIBLE", "NOMINAL_SEULEMENT",
              "NON_RETROUVE", "NON_TESTE", "ABSENT", "BLACKLISTE", "FICHIER_ABSENT"]
FILES = {
    "normal": {"OUTRIDER": "outrider/outrider_htseq_all.tsv", "FRASER": "fraser/fraser.tab"},
    "hyper":  {"OUTRIDER": "outrider_hyper/outrider_htseq_all.tsv", "FRASER": "fraser_hyper/fraser_results_all.tsv"},
}


# ── Entrées ───────────────────────────────────────────────────────────────────
def read_truth(path):
    rows, cols = [], None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = [p.strip() for p in line.split("\t")]
        if cols is None:
            cols = [c.lower() for c in parts]; continue
        rec = dict(zip(cols, parts + [""] * (len(cols) - len(parts))))
        if not rec.get("sample"):                      # ligne de continuation
            if rows:
                rows[-1]["variant"] += " | " + " ".join(p for p in parts if p)
            continue
        rows.append({"sample": re.search(CORE_RE, rec["sample"]).group(1) if re.search(CORE_RE, rec["sample"]) else rec["sample"],
                     "tool": rec.get("tool", "").upper(), "gene": rec.get("gene", ""),
                     "variant": rec.get("variant", ""), "effet": rec.get("effet", "")})
    return pd.DataFrame(rows)


def grch38_positions(variant):
    """Positions génomiques explicitement GRCh38/hg38 (les autres builds sont ignorés)."""
    out = []
    for m in re.finditer(r"chr([0-9XYM]+)\((GRCh38|hg38)\):g\.(\d+)", variant, flags=re.I):
        out.append((m.group(1).upper(), int(m.group(3))))
    return out


def load_blacklist(path):
    bl = {}
    p = Path(path) if path else None
    if not p or not p.is_file():
        return bl
    for line in p.read_text(errors="replace").splitlines():
        f = line.split("\t")
        if len(f) >= 2 and re.search(CORE_RE, f[0]) and not line.startswith("#"):
            bl.setdefault(re.search(CORE_RE, f[0]).group(1), set()).add(f[1].strip().lower())
    return bl


def gene_ids_from_gtf(gtfs, genes):
    """Symbole -> ensemble d'ENSG (sans version), lu sur les lignes 'gene' des GTF."""
    want, out = set(genes), {g: set() for g in genes}
    for gtf in gtfs:
        if not gtf or not Path(gtf).is_file():
            continue
        with open(gtf) as fh:
            for l in fh:
                if l.startswith("#"):
                    continue
                f = l.split("\t", 9)
                if len(f) < 9 or f[2] != "gene":
                    continue
                n = re.search(r'gene_name "([^"]+)"', f[8]); i = re.search(r'gene_id "([^"]+)"', f[8])
                if n and i and n.group(1) in want:
                    out[n.group(1)].add(i.group(1).split(".")[0])
    return out


def read_filtered(path, samples, usecols_wanted):
    """Lecture par blocs, en ne gardant que les échantillons de la liste."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    cols = [c for c in header if c in usecols_wanted]
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=cols, chunksize=2_000_000, low_memory=False):
        chunk["core"] = chunk["sampleID"].fillna("").astype(str).str.extract(CORE_RE, expand=False)
        parts.append(chunk[chunk["core"].isin(samples)])
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols + ["core"])


# ── Évaluation ────────────────────────────────────────────────────────────────
def num(df, col):
    """Colonne numérique, ou NaN si absente."""
    return pd.to_numeric(df[col], errors="coerce") if col in df.columns else pd.Series(float("nan"), index=df.index)


def md_table(df, floatfmt="{:.2e}"):
    """Tableau Markdown sans dépendance (pas de tabulate)."""
    def f(v):
        if isinstance(v, float):
            return "" if pd.isna(v) else floatfmt.format(v)
        return str(v)
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "| " + " | ".join("---" for _ in df.columns) + " |"
    body = ["| " + " | ".join(f(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([head, sep] + body)


def classify(padj, stat, eff, thr_eff, p, a):
    if pd.isna(p):
        return "NON_TESTE"
    if pd.notna(padj) and padj < a.padj:
        return "RETROUVE" if abs(eff) >= thr_eff else "SIGNIFICATIF_EFFET_FAIBLE"
    return "NOMINAL_SEULEMENT" if p < 0.05 else "NON_RETROUVE"


def eval_outrider(df, ev, ensg, a):
    sub = df[df["core"] == ev.sample]
    if sub.empty:
        return {"statut": "ABSENT"}
    padj_col = "padjValue" if "padjValue" in sub.columns else "padjust"
    sub = sub.assign(_p=num(sub, "pValue"), _q=num(sub, padj_col), _l=num(sub, "l2fc"), _z=num(sub, "zScore"))
    if "hgnc_symbol" in sub.columns and sub["hgnc_symbol"].eq(ev.gene).any():
        is_gene = sub["hgnc_symbol"].fillna("").eq(ev.gene)
    else:
        is_gene = sub["geneID"].fillna("").astype(str).str.split(".").str[0].isin(ensg)
    best = None
    for sid, s in sub.groupby("sampleID"):
        tested = s[s["_p"].notna()]
        g = s[is_gene.loc[s.index]]
        if g.empty:
            r = {"statut": "NON_TESTE", "echantillon": sid}
        else:
            row = g.sort_values("_p", na_position="last").iloc[0]
            rank = int((tested["_p"] < row["_p"]).sum() + 1) if pd.notna(row["_p"]) else None
            st = classify(row["_q"], None, row["_l"], a.l2fc, row["_p"], a)
            if st == "RETROUVE" and a.max_rank and rank and rank > a.max_rank:
                st = "SIGNIFICATIF_EFFET_FAIBLE"
            r = {"statut": st, "echantillon": sid, "pValue": row["_p"], "padj": row["_q"],
                 "effet_mesure": f"log2FC={row['_l']:.2f} z={row['_z']:.1f}" if pd.notna(row["_l"]) else "",
                 "sens": ("baisse" if row["_l"] < 0 else "hausse") if pd.notna(row["_l"]) else "",
                 "rang": rank, "n_testes": len(tested)}
        if best is None or RANK_ORDER.index(r["statut"]) < RANK_ORDER.index(best["statut"]):
            best = r
    return best


def eval_fraser(df, ev, a):
    sub = df[df["core"] == ev.sample]
    if sub.empty:
        return {"statut": "ABSENT"}
    padj_col = "padjust" if "padjust" in sub.columns else "padjValue"
    sub = sub.assign(_p=num(sub, "pValue"), _q=num(sub, padj_col), _d=num(sub, "deltaPsi"))
    sym = sub["hgncSymbol"].fillna("").astype(str) if "hgncSymbol" in sub.columns else pd.Series("", index=sub.index)
    is_gene = sym.str.split(r"[;,]").apply(lambda xs: ev.gene in [x.strip() for x in xs])
    pos = grch38_positions(ev.variant)
    best = None
    for sid, s in sub.groupby("sampleID"):
        tested = s[s["_p"].notna()]
        g = s[is_gene.loc[s.index] & s["_p"].notna()]
        by_pos = False
        if g.empty and pos and {"seqnames", "start", "end"} <= set(s.columns):
            # repli : jonctions proches du variant (annotation de gène absente ou différente)
            chrom_s = tested["seqnames"].fillna("").astype(str).str.replace("chr", "", regex=False).str.upper()
            near = pd.Series(False, index=tested.index)
            for c, p in pos:
                st_, en_ = pd.to_numeric(tested["start"], errors="coerce"), pd.to_numeric(tested["end"], errors="coerce")
                near |= (chrom_s == c) & ((st_ - p).abs().le(a.window) | (en_ - p).abs().le(a.window)
                                          | ((st_ <= p) & (en_ >= p)))
            g, by_pos = tested[near], True
        if g.empty:
            r = {"statut": "NON_TESTE", "echantillon": sid}
        else:
            row = g.sort_values("_p").iloc[0]
            rank = int((tested["_p"] < row["_p"]).sum() + 1)
            st = classify(row["_q"], None, row["_d"], a.dpsi, row["_p"], a)
            if st == "RETROUVE" and a.max_rank and rank > a.max_rank:
                st = "SIGNIFICATIF_EFFET_FAIBLE"
            n_sig = int(((g["_q"] < a.padj) & (g["_d"].abs() >= a.dpsi)).sum())
            dist = ""
            if pos and "start" in g.columns:
                chrom = str(row.get("seqnames", "")).replace("chr", "").upper()
                d = [min(abs(p - int(row["start"])), abs(p - int(row["end"]))) for c, p in pos if c == chrom]
                dist = min(d) if d else ""
            r = {"statut": st, "echantillon": sid, "pValue": row["_p"], "padj": row["_q"],
                 "effet_mesure": f"dPsi={row['_d']:.2f} {row.get('type', '')} "
                                 f"{row.get('seqnames', '')}:{row.get('start', '')}-{row.get('end', '')}"
                                 + (f" (par position, ±{a.window} pb)" if by_pos else ""),
                 "sens": ("baisse" if row["_d"] < 0 else "hausse") if pd.notna(row["_d"]) else "",
                 "rang": rank, "n_testes": len(tested), "jonctions_sig_gene": n_sig,
                 "distance_variant_pb": dist}
        if best is None or RANK_ORDER.index(r["statut"]) < RANK_ORDER.index(best["statut"]):
            best = r
    return best


def evaluate(run_dir, truth, branches, a, blacklist, ensg_map):
    out = []
    samples = set(truth["sample"])
    for br in branches:
        cache = {}
        for tool in ("OUTRIDER", "FRASER"):
            f = Path(run_dir) / FILES[br][tool]
            if tool in set(truth["tool"]) and f.is_file():
                want = {"sampleID", "geneID", "hgnc_symbol", "pValue", "padjValue", "padjust", "l2fc", "zScore"} \
                    if tool == "OUTRIDER" else \
                    {"sampleID", "hgncSymbol", "seqnames", "start", "end", "type", "pValue", "padjust", "deltaPsi"}
                print(f"[INFO] lecture {f} …", file=sys.stderr)
                cache[tool] = read_filtered(f, samples, want)
        for ev in truth.itertuples():
            tools = ("OUTRIDER", "FRASER") if ev.tool in ("LES DEUX", "BOTH") else (ev.tool,)
            for tool in tools:
                base = {"branche": br, "sample": ev.sample, "outil": tool, "gene": ev.gene,
                        "effet_attendu": ev.effet, "variant": ev.variant}
                if tool.lower() in blacklist.get(ev.sample, set()) or "all" in blacklist.get(ev.sample, set()):
                    out.append({**base, "statut": "BLACKLISTE"}); continue
                if tool not in cache:
                    out.append({**base, "statut": "FICHIER_ABSENT"}); continue
                r = eval_outrider(cache[tool], ev, ensg_map.get(ev.gene, set()), a) if tool == "OUTRIDER" \
                    else eval_fraser(cache[tool], ev, a)
                out.append({**base, **r})
    return pd.DataFrame(out)


# ── Rapport ───────────────────────────────────────────────────────────────────
def report(res, a, title):
    lines = [f"# {title}", "", f"Critères : padj < {a.padj} ; OUTRIDER |log2FC| ≥ {a.l2fc} ; "
             f"FRASER |ΔΨ| ≥ {a.dpsi}" + (f" ; rang ≤ {a.max_rank}" if a.max_rank else ""), ""]
    for br, d in res.groupby("branche"):
        n = len(d); ok = (d["statut"] == "RETROUVE").sum()
        lines.append(f"## Branche {br} : {ok}/{n} événements retrouvés ({100 * ok / n:.0f} %)")
        lines.append("")
        lines.append(md_table(d.groupby(["outil", "effet_attendu"])["statut"].value_counts()
                              .unstack(fill_value=0).reset_index()))
        lines.append("")
    cols = [c for c in ["branche", "sample", "outil", "gene", "effet_attendu", "statut", "statut_reference",
                        "padj", "effet_mesure", "sens", "rang", "n_testes", "distance_variant_pb"] if c in res.columns]
    lines += ["## Détail", "", md_table(res[cols])]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--run", required=True, help="dossier pipeline_v0 évalué (candidat)")
    ap.add_argument("--compare", help="dossier pipeline_v0 de référence (non-régression)")
    ap.add_argument("--branch", choices=["normal", "hyper", "both"], default="both")
    ap.add_argument("--site-paths", default=str(Path(__file__).resolve().parent.parent / "template" / "site_paths.yml"))
    ap.add_argument("--gtf", nargs="*", help="GTF pour symbole -> ENSG (défaut : gtfFile et gtf_refseq de site_paths)")
    ap.add_argument("--blacklist", help="blacklist unifiée (défaut : clé blacklist de site_paths)")
    ap.add_argument("--padj", type=float, default=0.05)
    ap.add_argument("--l2fc", type=float, default=1.0)
    ap.add_argument("--dpsi", type=float, default=0.1)
    ap.add_argument("--max-rank", type=int, default=None)
    ap.add_argument("--window", type=int, default=10000,
                    help="FRASER : si aucune jonction n'est annotée au gène, chercher à ± N pb du variant GRCh38 (défaut 10 000)")
    ap.add_argument("--out", default="rapport_validation")
    ap.add_argument("--strict", action="store_true")
    a = ap.parse_args()

    sp = yaml.safe_load(open(a.site_paths)) if Path(a.site_paths).is_file() else {}
    truth = read_truth(a.truth)
    branches = ["normal", "hyper"] if a.branch == "both" else [a.branch]
    blacklist = load_blacklist(a.blacklist or sp.get("blacklist", ""))
    gtfs = a.gtf if a.gtf is not None else [sp.get("gtfFile", ""), sp.get("gtf_refseq", "")]
    ensg = gene_ids_from_gtf(gtfs, truth.loc[truth["tool"].isin(["OUTRIDER", "LES DEUX", "BOTH"]), "gene"])
    print(f"[INFO] {len(truth)} événements attendus, {truth['sample'].nunique()} échantillons", file=sys.stderr)

    res = evaluate(a.run, truth, branches, a, blacklist, ensg)
    fail = False
    if a.compare:
        ref = evaluate(a.compare, truth, branches, a, blacklist, ensg)
        key = ["branche", "sample", "outil", "gene"]
        res = res.merge(ref[key + ["statut", "rang"]].rename(columns={"statut": "statut_reference", "rang": "rang_reference"}),
                        on=key, how="left")
        lost = res[(res["statut_reference"] == "RETROUVE") & (res["statut"] != "RETROUVE")]
        if len(lost):
            print(f"[ALERTE] {len(lost)} événement(s) retrouvé(s) par la référence et perdu(s) :", file=sys.stderr)
            print(lost[key + ["statut_reference", "statut"]].to_string(index=False), file=sys.stderr)
            fail = True
    total_missed = res[res["effet_attendu"].str.contains("total", case=False, na=False) & (res["statut"] != "RETROUVE")]
    if len(total_missed):
        fail = True

    res.to_csv(f"{a.out}.tsv", sep="\t", index=False)
    Path(f"{a.out}.md").write_text(report(res, a, f"Validation sur événements connus — {a.run}"), encoding="utf-8")
    print(f"[INFO] résultats : {a.out}.tsv, {a.out}.md", file=sys.stderr)
    print(res.groupby("branche")["statut"].value_counts().to_string())
    return 1 if (a.strict and fail) else 0


if __name__ == "__main__":
    sys.exit(main())
