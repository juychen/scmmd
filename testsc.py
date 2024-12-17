import matplotlib.pyplot as plt
import anndata
import scanpy as sc

from ALLCools.clustering import significant_pc_test
from ALLCools.plot import *


adata = anndata.read_h5ad('../input/10X.Neuron.h5ad')
