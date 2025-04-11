import pandas as pd
import re
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
import h5py
import shutil
import numpy as np

def predict_methylation_contributions(coverage_A, methylated_A, total_methylated_B, global_rate_B, gene_A_expression=None):
    # Calculate methylation rates
    rates_A = np.array(methylated_A) / np.array(coverage_A)
    
    # Expected methylated reads in Data A
    expected_A = rates_A * coverage_A
    
    # Adjust for Data B's global rate (if rates differ)
    total_expected_A = expected_A.sum()
    adjustment_factor = global_rate_B / (total_expected_A / sum(coverage_A))
    adjusted_expected = expected_A * adjustment_factor
    
    # Incorporate Gene A expression (if provided)
    if gene_A_expression is not None:
        # Assume gene_A_expression is a list of values per cell type
        # Normalize to weights (e.g., linear/logistic relationship)
        gene_weights = np.array(gene_A_expression) / sum(gene_A_expression)
        adjusted_expected *= gene_weights  # Scale contributions by gene expression
    
    # Rescale to total_methylated_B
    contributions = (adjusted_expected / adjusted_expected.sum()) * total_methylated_B
    return contributions
