#!/bin/bash
#SBATCH --partition=rohani_p
#SBATCH --ntasks=1
#SBATCH --job-name=construct_tree
#SBATCH --output=construct_tree_%A_%a.out
#SBATCH --error=construct_tree_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G

# Load the module
module load IQ-TREE/2.2.2.6-gompi-2022a

# Define variables
INPUT_FILE=$1
OUTPUT_DIR="/home/oa36725/scripts/Model1/phylogeny/trees/hiv/tree_output_$(basename $INPUT_FILE .fasta)"
OUTPUT_PREFIX="${OUTPUT_DIR}/tree"

# Create output directory
mkdir -p $OUTPUT_DIR

# Run IQ-TREE
iqtree2 -s $INPUT_FILE -m MFP -bb 1000 -alrt 1000 -nt AUTO -pre $OUTPUT_PREFIX
