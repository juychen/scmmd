#!/bin/bash
ts=$(date +"%Y-%m-%d-%H-%M-%S")
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2
nohup python darmotif.py --celltype Astro --region AMY > /home/junyichen/logs/darmotif_Astro_AMY_$ts.log 2>&1 &
nohup python darmotif.py --celltype OPC --region AMY > /home/junyichen/logs/darmotif_OPC_AMY_$ts.log 2>&1 &
nohup python darmotif.py --celltype Neuron --region HIP > /home/junyichen/logs/darmotif_Neuron_HIP_$ts.log 2>&1 &
nohup python darmotif.py --celltype OPC --region HIP > /home/junyichen/logs/darmotif_OPC_HIP_$ts.log 2>&1 &
nohup python darmotif.py --celltype Oligo --region HIP > /home/junyichen/logs/darmotif_Oligo_HIP_$ts.log 2>&1 &
nohup python darmotif.py --celltype Neuron --region PFC > /home/junyichen/logs/darmotif_Neuron_PFC_$ts.log 2>&1 &
nohup python darmotif.py --celltype Oligo --region PFC > /home/junyichen/logs/darmotif_Oligo_PFC_$ts.log 2>&1 &
nohup python darmotif.py --celltype Astro --region PFC > /home/junyichen/logs/darmotif_Astro_PFC_$ts.log
