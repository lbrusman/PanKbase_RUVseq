# Make sure file with contrast specifications can be read in properly
input_vars="pankbase_new_contrasts.csv"

dos2unix "$input_vars" # Change to unix

sed -i -e '$a\' "$input_vars" # Add blank line at end of csv if not already there

# Download donor metadata from PanKbase
mkdir -p metadata_files

wget -O metadata_files/pankbase_donor.tar.gz 'https://pankbase-data-v1.s3.us-west-2.amazonaws.com/analysis_resources/meta_data_table/human-donor_update_251210_v3.tar.gz'
tar -xzvf metadata_files/pankbase_donor.tar.gz -C metadata_files

# This is where the file was downloaded to/unzipped
meta_donor_path=metadata_files/human-donor_update_251210_v3/pankbase_human_donor.txt
export meta_donor_path

# Get biosample metadata url from PanKbase
biosample_url=https://pankbase-data-v1.s3.us-west-2.amazonaws.com/analysis_resources/meta_data_table/biosamples.2_2025-11-06.txt
export biosample_url

# Get cell counts table
cell_counts_file=metadata_files/sc_object_cellcounts.csv
export cell_counts_file

# Get chemistry metadata file
chem_meta_file=metadata_files/chemistry_metadata.csv
export chem_meta_file



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


# Create outdir to save all files - must have slash at end!
outdir=DE_outputs/

mkdir -p $outdir
export outdir


# Wait 2 seconds between submitting each job
sleep 2

# Run bash script that submits a separate job for each contrast (line in pankbase_new_contrasts.csv)
sbatch run_RUVseq_with_fGSEA.sh
done < <(tail -n +2 pankbase_new_contrasts.csv)
