#!/usr/bin/env python3
"""
check_run_outputs.py — le run est-il allé jusqu'au bout, et chaque fonction
a-t-elle produit sa sortie ? (test de fumée, run pilote, ou run de production)

    python3 tools/check_run_outputs.py --run-dir <dossier du run (contient pipeline_v0/)> \\
        [--workdir <dossier de travail Snakemake : <pipeline>/<run>>] [--snakemake-log <log>]

Contrôles (OK / KO / À VOIR) :
  fin          Snakemake terminé sans erreur (si --snakemake-log)
  versions     pipeline_versions.tsv : ligne « pipeline », aucun outil « not found »
  evenements   log/events/*.jsonl : toutes les tâches terminées avec le code 0 ;
               pas de fichier parasite « {log.run_info} »
  metriques    qc_summary.tsv et qc_warnings.tsv présents et non vides
  annotation   tables OUTRIDER annotées sans gène non testé (p-value vide)
  gnomad       livrables par échantillon : colonne loeuf présente et renseignée
  nettoyage    .cleanup_done présent, analysis_input/ supprimé
  deepvariant  un VCF chrX par BAM (si le module a tourné)
  version      même version dans run_info.json, pipeline_versions.tsv, qc_summary.tsv
               et les zips de livrables ; ni « -dirty » ni inconnue
  retour       marqueur .results_synced (information ; désactivé en développement)
Code de sortie : 1 si au moins un KO.
"""
import argparse
import io
import json
import sys
import zipfile
from pathlib import Path

RES = []


def check(name, status, detail):
    RES.append((name, status, detail))


def c_fin(log):
    if not log:
        return check("fin", "À VOIR", "--snakemake-log non fourni")
    p = Path(log)
    if not p.is_file():
        return check("fin", "KO", f"{p} introuvable")
    txt = p.read_text(errors="replace")
    if "Error in rule" in txt or "Exiting because a job execution failed" in txt or "MissingInputException" in txt:
        err = [l for l in txt.splitlines() if "Error in rule" in l or "Exception" in l][:3]
        return check("fin", "KO", "erreur Snakemake : " + " ; ".join(err))
    if "steps (100%) done" in txt or "Nothing to be done" in txt:
        return check("fin", "OK", "terminé sans erreur")
    check("fin", "À VOIR", "ni erreur ni fin détectée (run en cours ?)")


def c_versions(pv):
    f = pv / "metrics" / "pipeline_versions.tsv"
    if not f.is_file():
        return check("versions", "KO", f"{f} absent")
    lines = f.read_text(errors="replace").splitlines()
    pipe = [l for l in lines if l.split("\t")[0].strip().lower() == "pipeline"]
    nf = [l.split("\t")[0] for l in lines if "not found" in l.lower()]
    if not pipe:
        return check("versions", "KO", "ligne « pipeline » absente (version git)")
    v = pipe[0].split("\t")[-1]
    if nf:
        return check("versions", "KO", f"pipeline={v} ; not found : {', '.join(nf)}")
    check("versions", "OK" if "inconnue" not in v else "À VOIR", f"pipeline={v} ; {len(lines) - 1} outils")


def c_events(wd):
    if not wd:
        return check("evenements", "À VOIR", "--workdir non fourni")
    wd = Path(wd)
    if (wd / "{log.run_info}").exists():
        return check("evenements", "KO", "fichier parasite « {log.run_info} » présent (ancien journal)")
    files = sorted((wd / "log" / "events").glob("*.jsonl"))
    if not files:
        return check("evenements", "KO", f"aucun fichier dans {wd}/log/events/")
    open_jobs, failed, n = {}, [], 0
    for f in files:
        for l in f.read_text(errors="replace").splitlines():
            try:
                r = json.loads(l)
            except json.JSONDecodeError:
                continue
            k = (r.get("rule"), json.dumps(r.get("wildcards", {}), sort_keys=True))
            if r.get("event") == "start":
                open_jobs[k] = r; n += 1
            elif r.get("event") == "end":
                open_jobs.pop(k, None)
                if r.get("exit_code", 0) != 0:
                    failed.append(f"{k[0]}{k[1]} code={r['exit_code']}")
    if failed:
        return check("evenements", "KO", f"{len(failed)} tâche(s) en échec : " + " ; ".join(failed[:3]))
    if open_jobs:
        return check("evenements", "À VOIR", f"{len(open_jobs)} tâche(s) sans fin : " +
                     " ; ".join(k[0] for k in list(open_jobs)[:3]))
    check("evenements", "OK", f"{n} tâches, {len(files)} règles, toutes terminées avec le code 0")


def c_metrics(pv):
    miss = [n for n in ("qc_summary.tsv", "qc_warnings.tsv")
            if not (pv / "metrics" / n).is_file() or (pv / "metrics" / n).stat().st_size == 0]
    check("metriques", "KO" if miss else "OK", "absent/vide : " + ", ".join(miss) if miss else "qc_summary et qc_warnings présents")


def c_annot(pv):
    files = sorted(pv.glob("outrider*/outrider_htseq_annot.tsv"))
    if not files:
        return check("annotation", "KO", "aucune table outrider*_annot")
    bad = []
    for f in files:
        with open(f, errors="replace") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            if "pValue" not in head:
                bad.append(f"{f.parent.name}: pas de colonne pValue"); continue
            i = head.index("pValue")
            n = sum(1 for l in fh if (l.rstrip("\n").split("\t") + [""] * (i + 1))[i] in ("", "NA"))
            if n:
                bad.append(f"{f.parent.name}: {n} lignes non testées")
    check("annotation", "KO" if bad else "OK", " ; ".join(bad) if bad else f"{len(files)} table(s) sans gène non testé")


def _headers_in_zip(z):
    """(membre, en-tête, 1re ligne) pour les membres tabulés ; xlsx via openpyxl si disponible."""
    out = []
    for m in z.namelist():
        low = m.lower()
        if low.endswith((".tsv", ".txt", ".csv")):
            lines = z.read(m).decode(errors="replace").splitlines()
            sep = "," if low.endswith(".csv") else "\t"
            if lines:
                out.append((m, lines[0].split(sep), lines[1].split(sep) if len(lines) > 1 else []))
        elif low.endswith(".xlsx"):
            try:
                import openpyxl
                wb = openpyxl.load_workbook(io.BytesIO(z.read(m)), read_only=True)
                rows = list(wb.active.iter_rows(max_row=2, values_only=True))
                if rows:
                    out.append((m, [str(c) for c in rows[0]], [str(c) for c in rows[1]] if len(rows) > 1 else []))
            except ImportError:
                pass
    return out


def c_gnomad(pv):
    zips = sorted(pv.glob("per_sample_hyper/*.zip"))
    if not zips:
        return check("gnomad", "KO", "aucun livrable per_sample_hyper/*.zip")
    seen = with_col = filled = 0
    for zp in zips:
        with zipfile.ZipFile(zp) as z:
            for m, head, row in _headers_in_zip(z):
                seen += 1
                if "loeuf" in head:
                    with_col += 1
                    v = row[head.index("loeuf")] if len(row) > head.index("loeuf") else ""
                    filled += v not in ("", "None", "nan", "NA")
    if not seen:
        return check("gnomad", "À VOIR", "livrables non lisibles (format non tabulé)")
    if not with_col:
        return check("gnomad", "KO", f"colonne loeuf absente des {seen} tables lues (gnomAD v4.1 non actif ?)")
    check("gnomad", "OK" if filled else "À VOIR",
          f"loeuf présent dans {with_col}/{seen} tables ; renseigné en 1re ligne dans {filled}")


def c_cleanup(pv):
    ok = (pv / ".cleanup_done").is_file() and not (pv / "analysis_input").exists()
    check("nettoyage", "OK" if ok else "KO",
          "nettoyage fait" if ok else ".cleanup_done absent ou analysis_input/ encore présent")


def c_deepvariant(pv):
    d = pv / "deepvariant_WES"
    if not d.is_dir():
        return check("deepvariant", "À VOIR", "module non exécuté (run_variant_calling désactivé ?)")
    nb = len(list((pv / "star").glob("*_Aligned.sortedByCoord.out.bam")))
    nv = len(list(d.glob("*/*_chrX.vcf.gz")))
    check("deepvariant", "OK" if nv and nv == nb else "KO", f"{nv} VCF chrX pour {nb} BAM")


def c_version(pv):
    found = {}
    ri = pv / "run_info.json"
    if ri.is_file():
        try:
            found["run_info.json"] = json.loads(ri.read_text()).get("pipeline_version", "")
        except json.JSONDecodeError:
            found["run_info.json"] = "illisible"
    f = pv / "metrics" / "pipeline_versions.tsv"
    if f.is_file():
        for l in f.read_text(errors="replace").splitlines():
            if l.split("\t")[0].strip().lower() == "pipeline":
                found["pipeline_versions.tsv"] = l.split("\t")[-1].strip()
    q = pv / "metrics" / "qc_summary.tsv"
    if q.is_file():
        lines = q.read_text(errors="replace").splitlines()
        head = lines[0].split("\t") if lines else []
        if "pipeline_version" in head:
            vals = {l.split("\t")[head.index("pipeline_version")] for l in lines[1:] if l.strip()}
            found["qc_summary.tsv"] = ";".join(sorted(vals))
    for zp in sorted(pv.glob("per_sample_hyper/*.zip")):
        with zipfile.ZipFile(zp) as z:
            if "PIPELINE_VERSION.txt" in z.namelist():
                found[f"{zp.name}"] = z.read("PIPELINE_VERSION.txt").decode().strip()
    if not found:
        return check("version", "KO", "version du pipeline absente de toutes les sorties")
    vals = set(found.values())
    detail = " ; ".join(f"{k}={v}" for k, v in found.items())
    if len(vals) > 1:
        return check("version", "KO", "versions différentes : " + detail)
    v = vals.pop()
    missing = [n for n in ("run_info.json", "pipeline_versions.tsv", "qc_summary.tsv") if n not in found]
    if v.endswith("-dirty") or "inconnue" in v or not v:
        return check("version", "KO", f"version non identifiable ({v})")
    check("version", "À VOIR" if missing else "OK",
          f"{v} dans {len(found)} sortie(s)" + (f" ; absente de : {', '.join(missing)}" if missing else ""))


def c_retour(rd):
    m = rd / ".results_synced"
    check("retour", "OK" if m.is_file() else "À VOIR",
          m.read_text().strip() if m.is_file() else "pas de .results_synced (normal en développement)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--workdir")
    ap.add_argument("--snakemake-log")
    a = ap.parse_args()
    rd = Path(a.run_dir).resolve(); pv = rd / "pipeline_v0"
    if not pv.is_dir():
        print(f"[ERROR] {pv} introuvable"); return 1
    c_fin(a.snakemake_log); c_versions(pv); c_events(a.workdir); c_metrics(pv); c_annot(pv)
    c_gnomad(pv); c_cleanup(pv); c_deepvariant(pv); c_version(pv); c_retour(rd)
    w = max(len(n) for n, _, _ in RES)
    print(f"Contrôle de fin de run — {rd.name}")
    for n, s, d in RES:
        print(f"  [{s:6s}] {n:<{w}}  {d}")
    ko = sum(s == "KO" for _, s, _ in RES)
    print(f"Bilan : {sum(s == 'OK' for _, s, _ in RES)} OK, {ko} KO, {sum(s == 'À VOIR' for _, s, _ in RES)} à voir")
    return 1 if ko else 0


if __name__ == "__main__":
    sys.exit(main())
