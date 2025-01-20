# %%
import pandas as pd
import pyBigWig
import pandas as pd# %%
from wmb import cemba, mm10
import glob as glob
import numpy as np
import os
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Calculate the methlyation level for a given region")
parser.add_argument("--region", type=str, required=False, help="Input path of the region bed file",
                    default="/data2st1/junyi/generegion_vM33/promoter.bed")
parser.add_argument("--input", type=str, required=False, help="Input path of the bw file folder",
                    default="/data2st1/junyi/methlyatlas/mCseq/data.nemoarchive.org/biccn/grant/u19_cemba/ecker/epigenome/cellgroup/mCseq3/mouse/processed/other/")
parser.add_argument("--output", type=str, required=False, help="Output path of the methylation level",
                    default="/data2st1/junyi/output/")
parser.add_argument("--meth_type", type=str, required=False, help="Genome path",default="CGN")
parser.add_argument("--genome", type=str, required=False, help="Genome path",default="/data2st1/junyi/ref/GRCm38.p6.genome.fa")

args = parser.parse_args()

bw_file_path = args.input
methtype = args.meth_type

bw_files_mseq = glob.glob(bw_file_path+"*"+methtype+"*bw")

cell_groups = []
methtypes =  []

genome = SeqIO.index(args.genome, "fasta")
use_gene_meta = pd.read_csv(args.region,sep='\t',header=None)
annotation = use_gene_meta.iloc[:,-1].mode()[0]

for bw_file in bw_files_mseq:
    cell_group = os.path.basename(bw_file).split('.')[0]
    cell_groups.append(cell_group)
    methtype = os.path.basename(bw_file).split('.')[1]
    methtypes.append(methtype)
    bw = pyBigWig.open(bw_file)

    # %%
    bw.isBigWig()

    # %%
    use_gene_meta

    # %%
    stats =[]
    CGN_values = []
    for index, row in use_gene_meta.iterrows():
        chrom = row[0]
        start = row[1]
        end = row[2]
        
        values = bw.values(chrom, start, end)
        chromosome = genome[chrom]
        sequence = chromosome.seq[start:end]

        arr_values = np.array(values)
        arr_values = np.nan_to_num(arr_values)
        try:
            beta_vale = arr_values.sum()/sequence.count('C')
        except ZeroDivisionError:
            beta_vale = 0
        stats.append(beta_vale)
        CGN_values.append(values)


    if not os.path.exists(f'output/{cell_group}'):
        os.makedirs(f'output/{cell_group}')
    plt.hist(stats, bins=100, color='skyblue', edgecolor='black')
    # Add labels and title
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of Random Data')
    plt.savefig(f'output/{cell_group}/{methtype}_{annotation}.png')
    plt.close()

    value_array = np.array(CGN_values)
    np.save(f'output/{cell_group}/{methtype}_{annotation}.npy',value_array)
    stats = np.array(stats)
    np.save(f'output/{cell_group}/{methtype}_{annotation}_beta.npy',stats)