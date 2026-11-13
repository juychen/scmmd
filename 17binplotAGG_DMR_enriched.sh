#!/bin/bash
# 17binplotAGG_DMR_enriched.sh
# Per-motif footprint aggregate plots (MC vs MW) for the DMR BINDetect results.
# Same idea as 17binplotAGGAMY/HIP/PFC.sh (which was applied to AMY_Jaspar26 etc.),
# here applied to BULK/BIND_DMR_enriched: for every motif bed (*_all.bed) run
#   TOBIAS PlotAggregate --signals <region>_MC_corrected.bw <region>_MW_corrected.bw
# and write <motif>_all.bed.pdf + <motif>_all.bed_signal.txt next to the bed.
set -uo pipefail

source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias
ulimit -n 65534

BULK=/data2st2/junyi/output/atac1112/tobiasbam/BULK
output_dir=$BULK/BIND_DMR_enriched
finish_dir=$output_dir/finish_list
log_dir=$output_dir/plotagg_logs
mkdir -p "$finish_dir" "$log_dir"

parallel=24

for region in AMY HIP PFC; do
    mc_sig=$BULK/corrected/${region}_MC_corrected.bw
    mw_sig=$BULK/corrected/${region}_MW_corrected.bw
    workdir=$output_dir/${region}_MC_footprints.bw
    marker=$finish_dir/plotagg_${region}.txt
    log=$log_dir/plotagg_${region}.log

    if [[ ! -s "$mc_sig" || ! -s "$mw_sig" ]]; then
        echo "[$(date '+%F %T')] ERROR: missing corrected bigwig for $region ($mc_sig / $mw_sig)"
        exit 1
    fi

    if [[ -s "$marker" ]]; then
        echo "[$(date '+%F %T')] skip $region ($marker exists)"
        continue
    fi

    expected=$(find "$workdir" -path "*/beds/*_all.bed" | wc -l)
    echo "[$(date '+%F %T')] START $region: $expected motif beds | MC=$mc_sig | MW=$mw_sig" | tee -a "$log"

    find "$workdir" -path "*/beds/*_all.bed" -print0 \
      | xargs -0 -P "$parallel" -I{} bash -c '
            bed="$1"; mc="$2"; mw="$3"
            [[ -s "${bed}.pdf" ]] && exit 0
            TOBIAS PlotAggregate --TFBS "$bed" \
                --signals "$mc" "$mw" \
                --output "${bed}.pdf" \
                --share_y both --plot_boundaries --signal-on-x \
                --output-txt "${bed}_signal.txt" \
                --verbosity 2
        ' _ {} "$mc_sig" "$mw_sig" >> "$log" 2>&1

    done_pdf=$(find "$workdir" -path "*/beds/*_all.bed.pdf" | wc -l)
    echo "[$(date '+%F %T')] DONE $region: $done_pdf / $expected pdf" | tee -a "$log"
    if [[ "$done_pdf" -eq "$expected" ]]; then
        touch "$marker"
    else
        echo "[$(date '+%F %T')] WARNING $region incomplete ($done_pdf/$expected), rerun the script to fill in the rest" | tee -a "$log"
    fi
done

touch "$finish_dir/plotagg_ALL_DONE"
echo "[$(date '+%F %T')] PlotAggregate complete. Plots are in $output_dir/<region>_MC_footprints.bw/*/beds/*.pdf"
