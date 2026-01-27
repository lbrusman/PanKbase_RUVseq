#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=4:00:00
#SBATCH --account=csd854
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH -o outs/job_%j.out # File to which STDOUT will be written
#SBATCH -e outs/job_%j.err # File to which STDERR will be written
#SBATCH --mail-type ALL           # Optional, Send mail when the job ends
#SBATCH --mail-user lebrusman@health.ucsd.edu   # Optional, Send mail to this address

# Point to your conda
source /tscc/nfs/home/lebrusman/miniconda3/etc/profile.d/conda.sh
conda activate pankbase_degs # Activate your specific conda env

Rscript --vanilla RUVseq_with_fGSEA.R \
    contrast_id="$contrast_id" contrast_var="$contrast_var" sample_col="$sample_col" control_grp="$control_grp" experimental_grp="$experimental_grp" \
    subset_on1="$subset_on1" subset_on2="$subset_on2" subset1="$subset1" subset2="$subset2" \
    covariates="$covariates" correlate_vars="$correlate_vars" outdir="$outdir" metadata="$metadata" \
    scale="$scale" cell_types="$cell_types" file_pattern="$file_pattern" 