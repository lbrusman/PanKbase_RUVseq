# Make sure file with contrast specifications can be read in properly
input_vars="pankbase_new_contrasts.csv"

dos2unix "$input_vars" # Change to unix

sed -i -e '$a\' "$input_vars" # Add blank line at end of csv if not already there


# Now run script for each contrast
while IFS=',' read -r contrast_id contrast_var sample_col control_grp experimental_grp subset_on1 subset_on2 subset1 subset2 covariates correlate_vars scale cell_types file_pattern
do

echo $contrast_id
echo $contrast_var

export contrast_id
export contrast_var
export sample_col
export control_grp
export experimental_grp
export subset_on1
export subset_on2
export subset1
export subset2
export covariates
export correlate_vars
export scale
export cell_types
export file_pattern


outdir=DE_outputs/

mkdir -p $outdir
export outdir

metadata=sc_object_cellcounts.csv
export metadata

sleep 2

sbatch run_RUVseq_with_fGSEA.sh
done < <(tail -n +2 pankbase_new_contrasts.csv)
