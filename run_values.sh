#!/bin/bash
ts=$(date +"%Y-%m-%d-%H-%M-%S")
source activate allcools
conda activate allcools
python beta_byregion.py