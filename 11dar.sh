#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2

# output path make if not exists
mkdir -p /data2st1/junyi/output/atac1112/darmr

# python 11darbysubset.py --input /data2st1/junyi/output/atac1112/3REGIONS_dmr.h5ads \
# --output /data2st1/junyi/output/atac1112/darmr \
# --celltype_column All \
# --method memento-binary \
# --region PFC

# python 11darbysubset.py --input /data2st1/junyi/output/atac1112/3REGIONS_dmr.h5ads \
# --celltype_column All \
# --output /data2st1/junyi/output/atac1112/darmr \
# --method memento-binary \
# --region HIP

# python 11darbysubset.py --input /data2st1/junyi/output/atac1112/3REGIONS_dmr.h5ads \
# --celltype_column All \
# --output /data2st1/junyi/output/atac1112/dar \
# --method memento-binary \
# --region AMY

python 11darbysubset.py --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_PFC_Glut.h5ad \
--output /data2st1/junyi/output/atac1112/dar/celltype.L2/wil \
--celltype_column celltype.L2 \
--method wilcoxon \
--region PFC

python 11darbysubset.py --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_Glut.h5ad \
--celltype_column celltype.L2 \
--output /data2st1/junyi/output/atac1112/dar/celltype.L2/wil \
--method wilcoxon \
--region HIP

python 11darbysubset.py --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_AMY_Glut.h5ad \
--celltype_column celltype.L2 \
--output /data2st1/junyi/output/atac1112/dar/celltype.L2/wil \
--method wilcoxon \
--region AMY

python 11darbysubset.py --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_PFC_GABA.h5ad \
--output /data2st1/junyi/output/atac1112/dar/celltype.L2/wil \
--celltype_column celltype.L2 \
--method wilcoxon \
--region PFC

python 11darbysubset.py --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_GABA.h5ad \
--celltype_column celltype.L2 \
--output /data2st1/junyi/output/atac1112/dar/celltype.L2/wil \
--method wilcoxon \
--region HIP

python 11darbysubset.py --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_AMY_GABA.h5ad \
--celltype_column celltype.L2 \
--output /data2st1/junyi/output/atac1112/dar/celltype.L2/wil \
--method wilcoxon \
--region AMY