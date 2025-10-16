# pipeline_rnaseq

##Pré-requis


""# 🧬 RNA-Seq Analysis Pipeline (Snakemake)

**Author:** Laura Do Souto Ferreira  
**Affiliation:** CHU Nantes  
**Pipeline type:** Automated RNA-Seq analysis with Snakemake  
**Last major revision:** 2025-10-14  

---

## 📖 Overview

This repository contains a **Snakemake-based RNA-Seq pipeline** used for clinical and research sequencing analyses.  
It performs alignment, quantification, quality control, outlier detection (OUTRIDER / FRASER), and metrics summarization.

The workflow is fully reproducible using **Conda environments** and modularized into independent **Snakemake rules** and **scripts** for post-processing and visualization.

---

## 🗂️ Repository Structure

## ⚙️ Installation & Prerequisites

### 1️⃣ Requirements
- **Linux** or HPC environment with Slurm/Singularity  
- **Conda (miniconda or mamba)**  
- **Snakemake ≥ 7.30**  
- **Graphviz** (for rulegraph & filegraph PNGs)

```bash
conda install -c bioconda snakemake
conda install -c conda-forge graphviz


## ⚙️ Template Files — Manual Adjustments Required (:x)

Before running the pipeline, **two template files must be checked and adapted** to match your environment.

### 🧩 `template/config_template.yml`

This YAML file serves as the **base configuration** for the pipeline.  
The `scripts/preprocess.py` script automatically fills in:

| Field | Replaced Automatically | Description |
|--------|-------------------------|-------------|
| `WORKDIR` | ✅ | Path to the working/output directory defined by `--workDir`. |
| `{SAMPLES}` | ✅ | List of FASTQ sample IDs detected in the `dataDir`. |
| `FASTQDIR` | ✅ | Absolute path to the directory containing FASTQ files. |
| `pipeline_dir` | ✅ | Absolute path to the root pipeline directory. |

All other parameters (below the `# Paramètres indépendants du JOB` section)  
**use absolute paths** specific to your local HPC environment or cluster setup.  

> ⚠️ **Reflection:**  
> These paths (`/dataref/...`, `/datawork/...`) are environment-specific.  
> If you deploy the pipeline elsewhere, you **must replace them manually** with the correct paths  
> to genome references, BED files, blacklist files, and annotation data.  
> They are marked conceptually as **`:x`** values — to be **edited for your infrastructure**.

---

### 🚀 `template/launch_template.sh`

This bash script is used to launch Snakemake and visualize pipeline graphs.

The placeholders inside the file are dynamically replaced by `preprocess.py`:

| Placeholder | Description |
|--------------|-------------|
| `ROOT_PIPELINE` | The absolute path to the pipeline root. |
| `PATH_TO_CONFIG` | The path to the generated config file (e.g. `launch_folder/config.yml`). |

However, a few elements **must be checked manually (:x):**

| Variable | Example Value | Action |
|-----------|----------------|--------|
| `conda_prefix` | `~/pipeline/RNASEQ/routine/conda_env` | :x: Update this to match the location of your Snakemake Conda environments. |
| `conda activate snakemake` | (top of script) | :x: Ensure the environment name matches your installed Snakemake environment. |

The rest of the file is managed automatically by Snakemake when `preprocess.py` replaces the placeholders.

---

### ✅ Summary

| Template | Automatically Modified | Manual Edits Needed |
|-----------|------------------------|---------------------|
| `config_template.yml` | Yes — `WORKDIR`, `{SAMPLES}`, `FASTQDIR`, `pipeline_dir` | Yes — reference and annotation paths (`:x`) |
| `launch_template.sh` | Yes — `ROOT_PIPELINE`, `PATH_TO_CONFIG` | Yes — conda environment path and activation line (`:x`) |

---

> 🧠 **Tip:**  
> After running `preprocess.py`, you can inspect the generated files in the `launch_folder`  
> to confirm that all replacements were made correctly before executing:
>
> ```bash
> bash launch_folder/launch.sh
> ```



## 🧩 Script: `scripts/preprocess.py`

### 📝 Description

This script prepares the working environment and configuration files required to launch the Snakemake pipeline.  
It helps structure input data (FASTQ files), set up sample directories, and generate the configuration file (`config.yml`) dynamically before running the workflow.

It is designed to make the pipeline launch reproducible and parameterized from the command line.

---

### ⚙️ Usage

```bash
python scripts/preprocess.py \
  --path /path/to/pipeline \
  --workDir /path/to/working_directory \
  --dataDir /path/to/fastq_dir \
  --samples samples.txt \


#### 📁 Expected Input Structure

dataDir/
├── SAMPLE_01_R1.fastq.gz
├── SAMPLE_01_R2.fastq.gz
├── SAMPLE_02_R1.fastq.gz
└── SAMPLE_02_R2.fastq.gz

#### 🧱 Output

After execution, the script generates:

File	Description
workDir/config.yml	Pipeline configuration file adapted from template/config_template.yml



