import pandas as pd
import re
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import h5py
import shutil
import pysam
from multiprocessing import Pool
import concurrent.futures
from sklearn.svm import SVC
from sklearn.preprocessing import LabelEncoder

def conclude_pycistargets(file_paths:list):

    """
    Concludes the PycisTargets analysis by generating a summary report and cleaning up temporary files.
    """
    dict_table = {}
    # for celltype in cell_types:
    motif_region = {}

    region_list = []
    explist = []
    motif_list = []

    for file_path in file_paths:

        with h5py.File(file_path, 'r') as f:
            # Open the HDF5 file

            # get experiment name
            expname = list(f.keys())[0]
            dataset = f[expname]  # Replace with your dataset name
            dict_dataset = {}
        
            table = f[expname]['motif_enrichment']['table'][:]
            for key in dataset.keys():
                dict_dataset[key] = dataset[key]
            dict_table[expname] = table

            # get the motif region
            for key in f[expname]['motif_hits']['region_set'].keys():
                reg_set = f[expname]['motif_hits']['region_set'][key][:]
                reg_set = [item.decode('utf-8') for item in reg_set]

                for reg in reg_set:
                        region_list.append(reg)
                        motif_list.append(key)
                        explist.append(expname)
                motif_region[key] = ",".join(reg_set)

    df_motif_region = pd.DataFrame({"motif":motif_list,"region":region_list,"key":explist})

    # get the TF and scores
    list_table = []
    list_key = []
    list_TFs = []

    df_TF_celltype = pd.DataFrame()

    for key in dict_table.keys():
        data = dict_table[key]
        for row in data:
            list_table.append(row[0].decode('utf-8'))
            list_key.append(key)
            TFs_tmp=[]
            for i in range(1, len(row)-2):
                if 'img' in row[i][0].decode('utf-8'):
                    img = row[i][0].decode('utf-8')
                    continue
                elif key in row[i][0].decode('utf-8'):
                    continue
                TFs_tmp+=([ tf.decode('utf-8') for tf in row[i]])

            df_tfs = pd.DataFrame({'TF':",".join(TFs_tmp).split(',')})
            df_tfs['NES'] = row[-2][0]
            df_tfs['AUC'] = row[-2][1]
            df_tfs['Rank'] = row[-2][2]
            df_tfs['key'] = key
            df_tfs['motif'] = row[0].decode('utf-8')
            df_tfs['hit'] = row[-1][0]
            df_tfs['logo'] = img

            list_TFs.append(set(TFs_tmp))
            df_TF_celltype = pd.concat([df_TF_celltype,df_tfs],axis=0)
            
    df_TF = pd.DataFrame({'motif':list_table,'key':list_key,'TFs':list_TFs})

    df_merged_motif = df_TF.merge(df_motif_region,left_on=['motif','key'],right_on=['motif','key'])
    df_score_unique=df_TF_celltype.drop(columns=['TF']).drop_duplicates()
    return df_merged_motif.merge(df_score_unique,left_on=['motif','key'],right_on=['motif','key'])

def annotate_region(df_input,bedfile,region_col='region', temp_dir='./tmp') -> pd.DataFrame:
    """
    Annotates the regions in the input DataFrame using a bed file.
    region_col: str
        The name of the column in df_input that contains the regions to be annotated.
        Format: chr:start-end
    bedfile: str
        The path to the bed file used for annotation.
    temp_dir: str
        The directory where temporary files will be stored.
    Returns:
        pd.DataFrame: A DataFrame with the annotated regions.

    """
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    df_region = df_input[region_col].str.split(r'[ ,!\-;:|]',expand=True).drop_duplicates(subset=[0,1,2])
    df_region.to_csv(os.path.join(temp_dir,'tmp_region.bed'), sep='\t', header=False, index=False)
    # Get the bedfile name and sort it
    bedbasename = os.path.basename(bedfile)
    bedbasename = bedbasename.replace(".bed","_sorted.bed")
    destbedname = os.path.join(temp_dir,bedbasename)
    #!sort -k1,1 -k2,2n {filename} > {destname}
    command_str = f"sort -k1,1 -k2,2n {bedfile} > {destbedname}"
    subprocess.run(command_str, shell=True)
    #!sort -k1,1 -k2,2n {filename} > {destname}
    destregionname = os.path.join(temp_dir,'tmp_region_sorted.bed')
    command_str = f"sort -k1,1 -k2,2n {os.path.join(temp_dir,'tmp_region.bed')} > {destregionname}"
    subprocess.run(command_str, shell=True)
    closest_name = os.path.join(temp_dir,'tmp_cloeset.bed')
    # !bedtools closest -a {filename} -b {destname} -D ref > {destname}
    command_str = f"bedtools closest -a {destregionname} -b {destbedname} -D ref > {closest_name}"
    subprocess.run(command_str, shell=True)

    df_closest = pd.read_csv(closest_name, sep='\t', header=None)
    df_closest.columns = ['chr','start','end','gchr','gstart','gend','score','strand','gene_name','gene_id','annotation','distance']
    # Annotate the regions based on the distance and strand
    # df_closest.loc[(df_closest.distance>2000) & (df_closest.strand=='+'),'annotation'] = 'distal'
    # df_closest.loc[(df_closest.distance>0) & (df_closest.distance<2000) & (df_closest.strand=='+'),'annotation'] = 'promoter'
    # df_closest.loc[(df_closest.distance<0) & (df_closest.strand=='+'),'annotation'] = 'downstream'
    # df_closest.loc[(df_closest.distance<-2000) & (df_closest.strand=='-'),'annotation']= 'distal'
    # df_closest.loc[(df_closest.distance<0) & (df_closest.distance>-2000) & (df_closest.strand=='-'),'annotation']= 'promoter'
    # df_closest.loc[(df_closest.distance>0) & (df_closest.strand=='-'),'annotation'] = 'downstream'

    
    shutil.rmtree(temp_dir)
    df_closest[region_col]= df_closest['chr'] + ':' + df_closest['start'].astype(str) + '-' + df_closest['end'].astype(str)
    # Only keep cloest
    df_sorted = df_closest.sort_values(by='distance', key=lambda x: abs(x))
    df_sorted = df_sorted.drop_duplicates(subset=[region_col], keep='first')

    df_result = df_input.merge(df_sorted[[region_col, 'gene_name','gene_id', 'distance','gstart','gend']], on=region_col, how='left')
    return df_result

def process_celltype(allcfile,regions):
    tabix_file = pysam.TabixFile(allcfile)
    list_mc_covs = []
    for region in regions:
        mc = 0
        cov = 0
        chr = region[0]
        start = region[1]
        end = region[2]
        for row in tabix_file.fetch(chr, start, end):
            columns = row.split('\t')
            mc += int(columns[4])
            cov += int(columns[5])
        list_mc_covs.append([mc, cov])

    return list_mc_covs

def mc_coverage_by_region(allcfiles:list,regions:list,ncpus=4):

    # Create a pool of workers
    with Pool(processes=ncpus) as pool:
        results = pool.starmap(process_celltype, [(allcfile, regions) for allcfile in allcfiles])
    return dict(zip(allcfiles, results))

# test_region = [['chr1', 6526617, 6526633],
#  ['chr1', 6527054, 6527089],
#  ['chr1', 12577161, 12577177],
#  ['chr1', 16878736, 16878786],
#  ['chr1', 18098940, 18099022]]

# allcs=['/data2st1/junyi/merged_allcs/PFC/Neuron.tsv.gz',
#  '/data2st1/junyi/merged_allcs/PFC/Astro-Epen.tsv.gz',
#  '/data2st1/junyi/merged_allcs/PFC/Vascular.tsv.gz',
#  '/data2st1/junyi/merged_allcs/PFC/Immune.tsv.gz',
#  '/data2st1/junyi/merged_allcs/PFC/OPC-Oligo.tsv.gz']

# results = mc_coverage_by_region(allcs,test_region,ncpus=4)

# print(results)

def label_transfer_svc(adata_sc_train, adata_sc_test,label_col='celltype.L1',embedding_col='X_pca',svc_kwargs={'kernel': 'linear', 'C': 1}):
    """
    Transfer labels using SVC
    """
    y_train = adata_sc_train.obs[label_col]
    le = LabelEncoder()
    y_train = le.fit_transform(y_train)
    X_train = adata_sc_train.obsm[embedding_col]
    clf = SVC(**svc_kwargs)
    clf.fit(X_train, y_train)
    X_test = adata_sc_test.obsm[embedding_col]
    y_test = clf.predict(X_test)
    y_test = le.inverse_transform(y_test)
    adata_sc_test.obs[label_col] = y_test
    return adata_sc_test
