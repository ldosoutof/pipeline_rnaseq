# --- helpers for version detection ------------------------------------------
import subprocess
from pathlib import Path

def get_version_from_env(env_yml, cmd):
    """
    Detect a tool version by searching for its binary inside the local conda_env folder.
    """
    print(f"[DEBUG] get_version_from_env called with: {env_yml}, {cmd}")
    bin_name = cmd.split()[0]

    # Ensure PIPELINE_DIR is a Path object
    base_dir = Path(PIPELINE_DIR) / "conda_env"
    print(f"[DEBUG] Searching envs in {base_dir}")

    found_binary = None
    for env_dir in sorted(base_dir.glob("*")):
        if not env_dir.is_dir():
            continue
        for binary in env_dir.glob(f"**/bin/{bin_name}"):
            if binary.exists():
                print(f"[DEBUG] ✅ Found binary {binary}")
                found_binary = binary
                break
        if found_binary:
            break

    if not found_binary:
        print(f"[WARN] ❌ Could not find {bin_name} in any conda_env directory")
        return f"/bin/sh: 1: {bin_name}: not found"

    # Try to run version command
    try:
        result = subprocess.run(
            f"{found_binary} {cmd[len(bin_name):]}",
            shell=True,
            capture_output=True,
            text=True,
            executable="/bin/bash"
        )
        output = (result.stdout + result.stderr).strip()
        version = output.split("\n")[0]
        print(f"[DEBUG] ✅ {bin_name} version detected: {version}")
        return version if version else f"{bin_name}: version not found"
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to run {cmd} for {bin_name}: {e}")
        return f"/bin/sh: 1: {bin_name}: not found"


def get_python_version():
    """Get Python version in current Snakemake environment"""
    import sys
    return ".".join(map(str, sys.version_info[:3]))


# --- rule to collect and save versions --------------------------------------
rule save_pipeline_versions:
    """
    Collect tool versions used in the pipeline and save them to TSV.
    """
    output:
        "benchmarks/versions/pipeline_versions.tsv"
    run:
        versions = {
            "fastqc": get_version_from_env("envs/fastqc_env.yml", "fastqc --version | head -1 | cut -d' ' -f2"),
            "fastp": get_version_from_env("envs/fastp_env.yml", "fastp --version"),
            "STAR": get_version_from_env("envs/star_env.yml", "STAR --version | head -1 | cut -d' ' -f2"),
            "samtools": get_version_from_env("envs/samtools_env.yml", "samtools --version | head -1 | cut -d' ' -f2"),
            "htseq-count": get_version_from_env("envs/htseq_env.yml", "htseq-count --version | head -1 | cut -d' ' -f2"),
            "kallisto": get_version_from_env("envs/count_env.yml", "kallisto version"),
            "Rscript": get_version_from_env("envs/htseq_env.yml", "Rscript --version | head -1"),
            "mosdepth": get_version_from_env("envs/metrics_env.yml", "mosdepth --version"),
            "multiqc": get_version_from_env("envs/metrics_env.yml", "multiqc --version"),
            "picard": get_version_from_env("envs/metrics_env.yml", "picard MarkDuplicates --version"),
            "python": get_python_version(),
        }

        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as f:
            f.write("tool\tversion\n")
            for tool, version in versions.items():
                f.write(f"{tool}\t{version}\n")

        print(f"[INFO] Wrote {len(versions)} tool versions to {output[0]}")

