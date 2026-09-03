import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging

# 设置pyscenic模块的日志级别为WARNING或更高
logging.getLogger('pyscenic').setLevel(logging.WARNING)
logging.basicConfig(level=logging.WARNING)

# Directory containing the loom files
input_directory = "/data1st1/junyi/output/sn0827L1"

# Define the output directory for saving results
output_directory = "/data1st1/junyi/output/sn0827L1/scenic_grn"

# Path to the TFs file
tfs_file = "/data1st2/yejun/scenic/mouse/mm_mgi_tfs.txt"

os.makedirs(output_directory, exist_ok=True)

# List all loom files in the specified directory
loom_files = [f for f in os.listdir(input_directory) if f.endswith(".loom")]

def run_pyscenic(loom_file: str):
    """Run pyscenic grn on one loom file, with a .processing flag file."""

    input_file = os.path.join(input_directory, loom_file)
    output_file = os.path.join(
        output_directory,
        f"adj_{loom_file.replace('.loom', '.tsv')}"
    )
    flag_file = os.path.join(
        output_directory,
        f"{os.path.splitext(loom_file)[0]}.processing"
    )

    # ---- if processing file exists, skip ----
    if os.path.exists(flag_file):
        print(f"Skipping {loom_file}: processing flag exists.")
        return f"Skipped: {loom_file}"

    # ---- else proceed ----
    with open(flag_file, "w") as f:
        f.write(f"Processing {loom_file}\n")

    # Construct the pyscenic command
    command = [
        "pyscenic", "grn",
        "--num_workers", "32",
        "-o", output_file,
        input_file,
        tfs_file,
    ]

    try:
        subprocess.run(command, check=True)
        print(f"Finished processing {loom_file}, results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {loom_file}: {e}")
    finally:
        # Remove processing flag file regardless of success/failure
        if os.path.exists(flag_file):
            os.remove(flag_file)

    return loom_file

# Max parallel jobs
max_parallel_jobs = 8

# ---- Only submit tasks for files without .processing flag / existing output ----
files_to_run = []
for loom_file in loom_files:
    flag_file = os.path.join(
        output_directory, f"{os.path.splitext(loom_file)[0]}.processing"
    )
    output_file = os.path.join(
        output_directory,
        f"adj_{loom_file.replace('.loom', '.tsv')}"
    )
    if os.path.exists(flag_file):
        print(f"Skipping {loom_file}: {flag_file} already exists.")
    elif os.path.exists(output_file):
        print(f"Skipping {loom_file}: output file {output_file} already exists.")
    else:
        files_to_run.append(loom_file)

print(f"Submitting {len(files_to_run)} tasks...")

with ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
    future_to_file = {
        executor.submit(run_pyscenic, loom_file): loom_file
        for loom_file in files_to_run
    }

    for future in as_completed(future_to_file):
        loom_file = future_to_file[future]
        try:
            future.result()
        except Exception as exc:
            print(f"{loom_file} generated an exception: {exc}")
