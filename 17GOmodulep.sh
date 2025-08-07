#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snap
python testclf.py --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gosummary_module/1_gosummary.csv |\
python testclf.py --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gosummary_module/2_gosummary.csv |\
python testclf.py --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gosummary_module/3_gosummary.csv |\
python testclf.py --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gosummary_module/4_gosummary.csv ;
python testclf.py --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gosummary_module/5_gosummary.csv |\
python testclf.py --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gosummary_module/6_gosummary.csv |\
python testclf.py --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gosummary_module/7_gosummary.csv |\
python testclf.py --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gosummary_module/8_gosummary.csv ;