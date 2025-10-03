#!/usr/bin/env python3
import os, shutil, sys

def backup_fraser_counts(fraser_output_dir):
    if not os.path.exists(fraser_output_dir):
        print(f"No FRASER folder at {fraser_output_dir} to backup.")
        return

    tmp_dir = fraser_output_dir.rstrip("/") + "_tmp"
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    shutil.copytree(fraser_output_dir, tmp_dir)
    shutil.rmtree(fraser_output_dir)
    os.makedirs(os.path.join(fraser_output_dir, "savedObjects", "Data_Analysis"),
                exist_ok=True)
    for sub in ["nonSplitCounts", "splitCounts"]:
        src = f"{tmp_dir}/savedObjects/Data_Analysis/{sub}"
        dst = f"{fraser_output_dir}/savedObjects/Data_Analysis/{sub}"
        if os.path.exists(src):
            shutil.copytree(src, dst)
    shutil.rmtree(tmp_dir)
    print("✅ Saved counts backup restored.")

if __name__ == "__main__":
    backup_fraser_counts(sys.argv[1])

