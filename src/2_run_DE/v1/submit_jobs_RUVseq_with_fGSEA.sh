#make sure file can be read in properly
input_vars="pankbase_new_contrasts.csv"

dos2unix "$input_vars" #change to unix

sed -i -e '$a\' "$input_vars" #add blank line at end of csv if not already there


#now run script for each contrast
while IFS=, read -r contrast_id contrast_var donor_col control_grp experimental_grp subset_on1 subset_on2 subset1 subset2 covariates correlate_vars scale file_pattern assay
do

echo $contrast_id
echo $contrast_var

export contrast_id
export contrast_var
export donor_col
export control_grp
export experimental_grp
export subset_on1
export subset_on2
export subset1
export subset2
export covariates
export correlate_vars
export scale
export file_pattern
export assay

indir=/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/outputs/250903_outs/pseudobulk_counts/
export indir

outdir=/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/outputs_251118/
mkdir -p $outdir
export outdir

metadata=/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/pankbase_sc_metadada_withcellcounts_wAAB.tsv
export metadata

sleep 2

sbatch run_RUVseq_with_fGSEA.sh
done < <(tail -n +2 pankbase_new_contrasts.csv)
