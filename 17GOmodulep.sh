#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snap
python testclf.py --column term_idname --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gostmodule_new/TFmodule1_gosummary.csv |\
python testclf.py --column term_idname --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gostmodule_new/TFmodule2_gosummary.csv |\
python testclf.py --column term_idname --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gostmodule_new/TFmodule3_gosummary.csv |\
python testclf.py --column term_idname --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gostmodule_new/TFmodule4_gosummary.csv ;
python testclf.py --column term_idname --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gostmodule_new/TFmodule5_gosummary.csv |\
python testclf.py --column term_idname --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gostmodule_new/TFmodule6_gosummary.csv |\
python testclf.py --column term_idname --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gostmodule_new/TFmodule7_gosummary.csv |\
python testclf.py --column term_idname --csv_test_path /data2st1/junyi/output/atac0627/snregulation/gostmodule_new/TFmodule8_gosummary.csv ;