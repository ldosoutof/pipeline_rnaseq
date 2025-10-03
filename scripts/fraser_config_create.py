#!/usr/bin/env python3
import os
import sys
import shutil
from collections import defaultdict

def extract_run_number(bam_path):
    # Extract RUNXX from path
    basename = os.path.basename(bam_path)
    parts = basename.split('_')
    for p in parts:
        if p.startswith("RUN"):
            try:
                return int(p.replace("RUN", ""))
            except ValueError:
                continue
    return 0

#def backup_fraser_counts(fraser_output_dir):
#    if not os.path.exists(fraser_output_dir):
#        print(f"No FRASER folder at {fraser_output_dir} to backup.")
#        return
#
#    tmp_dir = fraser_output_dir + "_tmp"
#
#    # ✅ Remove any leftover tmp directory from previous runs
#    if os.path.exists(tmp_dir):
#        shutil.rmtree(tmp_dir)
#
#    # Copy the full FRASER directory to a temporary backup
#    shutil.copytree(fraser_output_dir, tmp_dir)
#
#    # Remove the original and recreate the minimal folder structure
#    shutil.rmtree(fraser_output_dir)
#    os.makedirs(os.path.join(fraser_output_dir, "savedObjects", "Data_Analysis"),
#                exist_ok=True)
#    for sub in ["nonSplitCounts", "splitCounts"]:
#        src = f"{tmp_dir}/savedObjects/Data_Analysis/{sub}"
#        dst = f"{fraser_output_dir}/savedObjects/Data_Analysis/{sub}"
#        if os.path.exists(src):
#            shutil.copytree(src, dst)
#    shutil.rmtree(tmp_dir)
#    print("✅ Saved counts backup restored.")

def generate_config_fraser(output_filename, pattern, root_dir, fraser_output_dir, blacklist_file=None):
    # Backup existing FRASER counts
    #backup_fraser_counts(fraser_output_dir)

    # Load blacklist
    blacklist = set()
    if blacklist_file:
        with open(blacklist_file, 'r') as bf:
            for line in bf:
                parts = line.strip().split()
                if parts:
                    blacklist.add(parts[0])

    # Collect BAMs by sample_id
    sample_to_bams = defaultdict(list)
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if "RNASEQ" not in dirpath or "star" not in dirpath:
            continue
        for filename in filenames:
            if filename.endswith('Aligned.sortedByCoord.out.bam') and pattern in filename:
                #sample_id = filename.split('-')[0]
                if "-" in filename:
                    sample_id = filename.split("-")[0]
                elif "_" in filename:
                    sample_id = filename.split("_")[0]
                else:
                    sample_id=filename  # no dash or underscore
                if sample_id in blacklist or 'RUN17' in dirpath:
                    continue
                bam_path = os.path.join(dirpath, filename)
                run_num = extract_run_number(bam_path)
                sample_to_bams[sample_id].append((run_num, bam_path))

    # Keep only the most recent BAM per sample
    final_samples = {}
    for sample_id, bams in sample_to_bams.items():
        bams_sorted = sorted(bams, key=lambda x: x[0], reverse=True)
        final_samples[sample_id] = bams_sorted[0][1]  # most recent run

    # Write config
    with open(output_filename, 'w') as out:
        out.write("sampleID\tbamFile\tgroup\tgene\tpairedEnd\n")
        group_increment = 1
        import string
        import itertools
        gene_sequence = (f"{letter}{number}" for number in itertools.cycle(range(1, 1000)) for letter in string.ascii_uppercase)
        for sample_id, bam_path in sorted(final_samples.items()):
            gene = next(gene_sequence)
            paired_end = "TRUE"
            out.write(f"{sample_id}\t{bam_path}\t{group_increment}\t{gene}\t{paired_end}\n")
            group_increment += 1

    print(f"✅ FRASER config written to {output_filename}, {len(final_samples)} samples included.")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python fraser_config_create.py <output_file> <pattern> <root_dir> <fraser_output_dir> [blacklist_file]")
        sys.exit(1)

    output_file = sys.argv[1]
    pattern = sys.argv[2]
    root_dir = sys.argv[3]
    fraser_output_dir = sys.argv[4]
    blacklist_file = sys.argv[5] if len(sys.argv) > 5 else None

    generate_config_fraser(output_file, pattern, root_dir, fraser_output_dir, blacklist_file)

