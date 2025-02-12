# %%
import pandas as pd
import pyBigWig
import pandas as pd# %%
import glob as glob
import numpy as np
import os
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse
import copy

parser = argparse.ArgumentParser(description="Calculate the methlyation level for a given region")
parser.add_argument("--region", type=str, required=False, help="Input path of the region bed file",
                    default="/data2st1/junyi/generegion_vM23/genebody_selected.bed")
parser.add_argument("--input", type=str, required=False, help="Input path of the bw file folder",
                    default="/data2st1/junyi/methlyatlas/mCseq/data.nemoarchive.org/biccn/grant/u19_cemba/ecker/epigenome/cellgroup/mCseq3/mouse/processed/other/")
parser.add_argument("--output", type=str, required=False, help="Output path of the methylation level",
                    default="/data2st1/junyi/output/")
parser.add_argument("--meth_type", type=str, required=False, help="Methlylation type, can be CGN,CHN, based on your definetion of bw file",default="CGN")
parser.add_argument("--genome", type=str, required=False, help="Genome path",default="/data2st1/junyi/ref/GRCm38.p6.genome.fa")
parser.add_argument("--tile_length", type=int, required=False, help="Size of each tile",default=500)
parser.add_argument("--flanking", type=int, required=False, help="Flanking base base",default=250000)


args = parser.parse_args()

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(seq))

bw_file_path = args.input
methtype = args.meth_type
tile_length = args.tile_length
flanking = args.flanking

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
    meta_entries = []
    for index, row in use_gene_meta.iterrows():
        chrom = row[0]
        start = row[1]
        end = row[2]

        try:
        
            values = bw.values(chrom, start, end)
            chromosome = genome[chrom]
            if row[4] == '-':
                sequence = chromosome.seq[start:end]
                sequence = reverse_complement(sequence)
            else:
                sequence = chromosome.seq[start:end]

        except:
            print(f"Error in {chrom} {start} {end}")
            arr_values = np.zeros(end-start)
            stats.append(0)
            CGN_values.append(arr_values)
            meta_entries.append(row)
            continue

        arr_values = np.array(values)
        arr_values = np.nan_to_num(arr_values)
        if sequence.count('C') > 0:
            beta_vale = arr_values.sum()/sequence.count('C')
        else:
            beta_vale = 0
        stats.append(beta_vale)
        CGN_values.append(values)
        meta_entries.append(row)

    

        if flanking >0 and tile_length >0:
            extra_rows = []
            l_starts = np.arange(start - flanking, start, tile_length)
            l_ends = l_starts+tile_length
            r_starts = np.arange(end, end+flanking, tile_length)
            r_ends = r_starts+tile_length

            starts= np.concatenate([l_starts,r_starts])
            ends = np.concatenate([l_ends,r_ends])

            df_l = pd.DataFrame({'start':l_starts,'end':l_ends,'annotation':'upstream'})
            df_l['annotation'] = 'upstream_tile' 
            df_r = pd.DataFrame({'start':r_starts,'end':r_ends,'annotation':'downstream'})
            df_r['annotation'] = 'downstream_tile'
            df_extra = pd.concat([df_l,df_r])

            for index, extra_row in df_extra.iterrows():
                start = extra_row[0]
                end = extra_row[1]

                meta_row = copy.deepcopy(row)
                meta_row[2] = end
                meta_row[1] = start
                meta_row[7] = extra_row["annotation"]

                try:
                
                    values = bw.values(chrom, start, end)
                    chromosome = genome[chrom]
                    if meta_row[4] == '-':
                        sequence = chromosome.seq[start:end]
                        sequence = reverse_complement(sequence)
                    else:
                        sequence = chromosome.seq[start:end]

                except:
                    print(f"Error in {chrom} {start} {end}")
                    arr_values = np.zeros(end-start)
                    stats.append(0)
                    #CGN_values.append(arr_values)
                    meta_entries.append(row)
                    continue

                arr_values = np.array(values)
                arr_values = np.nan_to_num(arr_values)

                if sequence.count('C') > 0:
                    beta_vale = arr_values.sum()/sequence.count('C')
                else:
                    beta_vale = 0
                stats.append(beta_vale)
                #CGN_values.append(values)

                meta_entries.append(meta_row)




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
    df_meta = pd.DataFrame(meta_entries)
    df_meta.to_csv(f'output/{cell_group}/{methtype}_{annotation}_meta.csv',index=False)