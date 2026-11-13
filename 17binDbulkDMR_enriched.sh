#!/bin/bash
# 17binDbulkDMR_enriched.sh
# TOBIAS BINDetect on all DMRs (dmr_3regions_ext.bed) for the 6 pseudobulk footprint bigwigs
# (AMY / HIP / PFC x MC / MW), using the pycisTarget v10nr_clust DMR-enriched MEME motif set.
# One BINDetect run per region (MC vs MW), output in BULK/BIND_DMR_enriched/<sample>_footprints.bw/
set -uo pipefail

source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias

# BINDetect keeps one file handle per motif (1126) when writing the per-TF bed files;
# with the default soft limit of 1024 it dies with "Too many open files".
ulimit -n 65534

cd /data2st2/junyi/output/atac1112/tobiasbam

folder1="/data2st2/junyi/output/atac1112/tobiasbam/BULK/corrected/"   # MC footprints
folder2="/data2st2/junyi/output/atac1112/tobiasbam/BULK/corrected/"   # MW footprints

output_dir="/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched"
finish_dir="$output_dir/finish_list"
log_dir="$output_dir/logs"
motif_dir="$output_dir/motifs"
mkdir -p "$output_dir" "$finish_dir" "$log_dir" "$motif_dir"

peaks=/data2st1/junyi/output/atac1112/cCRE/dmr_3regions_ext.bed
genome=/data2st1/junyi/ref/GRCm38.p6.genome.fa
cores=64

# --- motif file -------------------------------------------------------------------------------
# The source MEME file has "MOTIF <id>" lines only (no name column). TOBIAS reads column 2 as the
# motif id and column 3 as the name, so without a name the "name" column of bindetect_results
# would be empty. Add the id as its own name: MOTIF <id> <id>.
motifs_src=/data1st2/hannan_25/data/Nanopore_processV1/nanopore_08_differential/summary/dmr_sets/DMR_pycis_pipeline/pycistarget_out/motif_use/v10nr_clust.from_hdf5_DMRpycisfimoenriched.meme
motifs="$motif_dir/v10nr_clust.from_hdf5_DMRpycisfimoenriched.named.meme"
if [[ ! -s "$motifs" ]]; then
    echo "[$(date '+%F %T')] building named motif file: $motifs"
    sed -E 's/^MOTIF (\S+)$/MOTIF \1 \1/' "$motifs_src" > "$motifs.tmp" && mv "$motifs.tmp" "$motifs"
fi
echo "[$(date '+%F %T')] motifs: $motifs ($(grep -c '^MOTIF' "$motifs") motifs)"
echo "[$(date '+%F %T')] peaks : $peaks ($(wc -l < "$peaks") regions)"
echo "[$(date '+%F %T')] genome: $genome"

# --- one BINDetect run per region -------------------------------------------------------------
for region in AMY HIP PFC; do
    mc_file="${region}_MC_footprints.bw"
    mw_file="${region}_MW_footprints.bw"
    bam_file="$folder1/$mc_file"

    if [[ ! -s "$folder1/$mc_file" || ! -s "$folder2/$mw_file" ]]; then
        echo "[$(date '+%F %T')] ERROR: missing footprint file for $region ($folder1/$mc_file or $folder2/$mw_file)"
        exit 1
    fi

    filename=$(basename "$bam_file")          # e.g. AMY_MC_footprints.bw
    outdir="$output_dir/$filename"
    marker="$finish_dir/$filename.txt"
    log="$log_dir/$filename.$(date +%Y%m%d_%H%M%S).log"

    if [[ -s "$marker" ]]; then
        echo "[$(date '+%F %T')] skip $region ($marker exists)"
        continue
    fi

    # keep a previous partial run out of the way instead of mixing files into it
    if [[ -d "$outdir" ]]; then
        mv "$outdir" "${outdir}.partial_$(date +%Y%m%d_%H%M%S)"
    fi

    echo "[$(date '+%F %T')] START $region : $filename vs $mw_file  (log: $log)"
    TOBIAS BINDetect \
        --motifs "$motifs" \
        --signals "$folder1/$mc_file" "$folder2/$mw_file" \
        --outdir "$outdir" \
        --cond-names MC MW \
        --naming id \
        --cores "$cores" \
        --genome "$genome" \
        --peaks "$peaks" \
        --verbosity 2 \
        > "$log" 2>&1

    if [[ $? -eq 0 && -s "$outdir/bindetect_results.txt" ]]; then
        touch "$marker"
        echo "[$(date '+%F %T')] DONE  $region -> $outdir  (motifs scored: $(($(wc -l < "$outdir/bindetect_results.txt")-1)))"
    else
        echo "[$(date '+%F %T')] FAILED $region, see $log"
        exit 1
    fi
done

touch "$finish_dir/ALL_DONE"
echo "[$(date '+%F %T')] Bindetect complete. Output saved to $output_dir."
