import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

# Directory containing the CSV files
input_directory = "/data2st2/junyi/output/stg1028/subsetsn"

# Define the output directory for saving results
output_directory = "/data2st2/junyi/output/stg1028/scenic_grn"

# Path to the TFs file
tfs_file = "/data1st2/yejun/scenic/mouse/mm_mgi_tfs.txt"

os.makedirs(output_directory, exist_ok=True)

# List all CSV files in the specified directory
csv_files = [f for f in os.listdir(input_directory) if f.endswith(".loom")]

def run_pyscenic(csv_file: str):
    """Run pyscenic on one loom file, with a .processing flag file."""

    input_file = os.path.join(input_directory, csv_file)
    output_file = os.path.join(
        output_directory,
        f"adj_{csv_file.replace('.loom', '.tsv')}"
    )
    flag_file = os.path.join(
        output_directory,
        f"{os.path.splitext(csv_file)[0]}.processing"
    )

    # ---- NEW LOGIC: if processing file exists, skip ----
    if os.path.exists(flag_file):
        print(f"Skipping {csv_file}: processing flag exists.")
        return f"Skipped: {csv_file}"

    # ---- else proceed ----
    with open(flag_file, "w") as f:
        f.write(f"Processing {csv_file}\n")

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
        print(f"Finished processing {csv_file}, results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {csv_file}: {e}")
    finally:
        # Remove processing flag file regardless of success/failure
        if os.path.exists(flag_file):
            os.remove(flag_file)

    return csv_file

# Max parallel jobs
max_parallel_jobs = 4

# ---- Only submit tasks for files without .processing flag ----
files_to_run = []
for csv_file in csv_files:
    flag_file = os.path.join(
        output_directory, f"{os.path.splitext(csv_file)[0]}.processing"
    )
    if os.path.exists(flag_file):
        print(f"Skipping {csv_file}: {flag_file} already exists.")
    else:
        files_to_run.append(csv_file)

print(f"Submitting {len(files_to_run)} tasks...")

with ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
    future_to_file = {
        executor.submit(run_pyscenic, csv_file): csv_file
        for csv_file in files_to_run
    }

    for future in as_completed(future_to_file):
        csv_file = future_to_file[future]
        try:
            future.result()
        except Exception as exc:
            print(f"{csv_file} generated an exception: {exc}")