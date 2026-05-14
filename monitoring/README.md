# RNA-seq Pipeline Monitoring — Grafana + Prometheus

## Architecture

```
Pipeline run
    └── rule generate_metrics
            ├── recup_metrics.py  →  qc_summary.tsv
            └── push_metrics.py   →  Pushgateway → Prometheus → Grafana
```

## Quick start

### 1. Start the stack

```bash
cd monitoring/
docker-compose up -d
```

Services:
| Service      | URL                    |
|-------------|------------------------|
| Grafana      | http://localhost:3001  |
| Prometheus   | http://localhost:9090  |
| Pushgateway  | http://localhost:9091  |

Login: `admin` / `rnaseq_admin`

The RNA-seq dashboard loads automatically.

### 2. Set the gateway URL in your config

```yaml
# config.yml
prometheus_gateway: "http://localhost:9091"
```

Metrics are pushed automatically at the end of each pipeline run (rule `generate_metrics`). If the Pushgateway is unavailable, the pipeline still completes — a warning is logged but no error is raised.

### 3. Backfill historical runs

To push metrics for all past runs at once:

```bash
python monitoring/push_metrics.py \
    --prod_root /datawork2/genetique/RNASeq/diag/prod \
    --gateway   http://localhost:9091 \
    --backfill
```

## Dashboard panels

### Per-Sample QC
| Panel | Metric | Thresholds |
|---|---|---|
| Total Reads | `rnaseq_nb_reads` | 🟡 20M  🟢 40M |
| Duplication Rate | `rnaseq_dup_pct` | 🟢 <30%  🟡 <50%  🔴 >50% |
| OUTRIDER / FRASER Events | `rnaseq_nb_outrider_events` / `_fraser_events` | — |
| DI-Green Expressed % | `rnaseq_di_green_pct` | 🔴 <70%  🟡 <85%  🟢 >85% |
| Haemoglobin TPM | `rnaseq_hba1/2_tpm`, `rnaseq_hbb_tpm`, `rnaseq_hba_total_tpm` | — |

### Run-Level Overview
Trend plots of average reads, average duplication rate, and average DI-green expression across all runs — useful for detecting batch effects or reagent lot changes.

### Pipeline Performance
| Panel | Metric |
|---|---|
| Rule Runtime | `rnaseq_rule_runtime_seconds` (from benchmarks) |
| Peak Memory | `rnaseq_rule_max_rss_mb` (from benchmarks) |
| Failed Rules | `rnaseq_rule_exit_code != 0` (from structured logs) |

### Sample Table
Full interactive table of all samples across all selected runs, with conditional colouring on duplication rate and DI-green expression.

## Variables (top of dashboard)
- **Run** — filter by one or more run IDs (e.g. `RUN48`)
- **Sample** — filter by one or more sample IDs within the selected run(s)

## Adding alerts

In Grafana → Alerting → Alert rules, create rules such as:
- Duplication > 50% → notify
- DI-green < 70% → notify
- Any failed rule → notify

Configure notification channels (email / Slack) under Alerting → Contact points.
