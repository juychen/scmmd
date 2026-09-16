import os
import subprocess
import multiprocessing
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# === PATHS ===
INPUT_DIR = "/data1st1/junyi/output/sn0827L1/"
GRN_DIR = "/data1st1/junyi/output/sn0827L1/scenic_grn"
OUTPUT_DIR = "/data1st1/junyi/output/sn0827L1/scenic_ctx"
ANNOTATIONS = "/data1st2/yejun/scenic/mouse/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
DATABASE = "/data1st2/yejun/scenic/mouse/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Dynamic CPU detection
workers = min(32, max(1, multiprocessing.cpu_count() - 4))

def fix_grn_columns(grn_path):
    """Ensure GRN file column names are correct for PySCENIC."""
    try:
        with open(grn_path, 'r') as f:
            header = f.readline().strip()
            columns = header.split('\t')

        target_columns = ['TF', 'target', 'importance']
        if columns != target_columns:
            df = pd.read_csv(grn_path, sep='\t')
            df.columns = [col.lower() for col in df.columns]
            df = df.rename(columns={'tf': 'TF'})
            df = df[['TF', 'target', 'importance']]
            df.to_csv(grn_path, sep='\t', index=False)
            print(f"✅ Fixed columns: {os.path.basename(grn_path)}")

        return True
    except Exception as e:
        print(f"⚠️ Column fix failed for {grn_path}: {e}")
        return False

def run_pyscenic_ctx(loom_file):
    """Run PySCENIC ctx on one input loom with concurrency-safe flagging."""
    expr_path = os.path.join(INPUT_DIR, loom_file)
    base_name = os.path.splitext(loom_file)[0]
    grn_path = os.path.join(GRN_DIR, f"adj_{base_name}.tsv")
    output_path = os.path.join(OUTPUT_DIR, f"ctx_{base_name}.csv")
    flag_file = os.path.join(OUTPUT_DIR, f"{base_name}.processing")

    # Skip if already processed or being processed
    if os.path.exists(flag_file):
        print(f"⏩ Skipping {loom_file}: already processing.")
        return f"Skipped: {loom_file}"
    if os.path.exists(output_path):
        print(f"✅ Skipping {loom_file}: output exists.")
        return f"Skipped (exists): {loom_file}"
    if not os.path.exists(grn_path):
        print(f"🚫 Skipping {loom_file}: GRN missing.")
        return f"Missing GRN: {loom_file}"

    # Write processing flag
    with open(flag_file, "w") as f:
        f.write(f"Processing {loom_file}\n")

    # Ensure valid GRN file
    if not fix_grn_columns(grn_path):
        print(f"🚫 Invalid GRN format: {loom_file}")
        os.remove(flag_file)
        return f"Invalid GRN: {loom_file}"

    # Construct command
    cmd = [
        "pyscenic", "ctx",
        grn_path,
        DATABASE,
        "--annotations_fname", ANNOTATIONS,
        "--expression_mtx_fname", expr_path,
        "--output", output_path,
        "--num_workers", str(workers)
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"✅ Completed: {loom_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running PySCENIC on {loom_file}: {e}")
    except Exception as e:
        print(f"⚠️ Exception during {loom_file}: {e}")
    finally:
        if os.path.exists(flag_file):
            os.remove(flag_file)

    return loom_file

# === PARALLEL EXECUTION ===
max_parallel_jobs = 4

loom_files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.loom')]
files_to_run = []

# Pre-check files
for loom_file in loom_files:
    base_name = os.path.splitext(loom_file)[0]
    output_file = os.path.join(OUTPUT_DIR, f"ctx_{base_name}.csv")
    flag_file = os.path.join(OUTPUT_DIR, f"{base_name}.processing")

    if os.path.exists(output_file):
        print(f"Skip: {loom_file} (output exists)")
    elif os.path.exists(flag_file):
        print(f"Skip: {loom_file} (flag exists)")
    else:
        files_to_run.append(loom_file)

print(f"🚀 Submitting {len(files_to_run)} tasks for PySCENIC ctx...")

with ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
    future_to_file = {executor.submit(run_pyscenic_ctx, f): f for f in files_to_run}

    for future in as_completed(future_to_file):
        loom_file = future_to_file[future]
        try:
            result = future.result()
            print(f"✅ Done: {result}")
        except Exception as exc:
            print(f"⚠️ {loom_file} raised an exception: {exc}")
