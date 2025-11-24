# PanKbase RUVseq
This repo houses the modular RUVseq pipeline for PanKbase.

Here are some additional details about the analysis, which was performed on a per-cell type basis.

## General filtering
Donors were filtered based on the following criteria:
1. Must have > 20 cells in that cell type
2. Cells must not be treated with any reagents (`treatments == "no_treatment"`)

Genes were filtered based on the following criteria:
1. Must have $\geq$ 5 counts in $\geq$ 25% of cells in that cell type

For a cell type/contrast to be attempted: 
1. There must be > 2 donors total left after filtering
2. If contrast variable is categorical, there must be $\geq$ 2 samples left in each group after filtering
3. If the contrast variable is categorical, there must > 1 unique group after filtering

Reasons pipeline could fail/exit early:
1. Initial DESeq failed (error message is returned)
2. Max number of latent variables (after accounting for number of known covariates and sample size) is zero
3. RUVseq failed (error message is returned)
4. k chosen is null or zero


## Contrast-specific parameters
All parameters are summarized in [pankbase_new_contrasts.csv](src/2_run_DE/pankbase_new_contrasts.csv)
- No diabetes vs. type 1 diabetes (ND_vs_T1D): 
    - Additional donor filtering:
        - None
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + diabetes.status.description`

- No diabetes vs. type 2 diabetes (ND_vs_T2D):
    - Additional donor filtering:
        - None
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + diabetes.status.description`

- Age (age):
    - Additional donor filtering:
        - Non-diabetic donors **only**
    - Initial DESeq formula:
        - `~ sex + bmi + ethnicity + chemistry + age`

- Sex (sex):
    - Additional donor filtering:
        - Non-diabetic donors **only**
    - Initial DESeq formula:
        - `~ age + bmi + ethnicity + chemistry + sex`

- BMI (BMI):
    - Additional donor filtering:
        - Non-diabetic donors **only**
    - Initial DESeq formula:
        - `~ age + sex + ethnicity + chemistry + bmi`

- AAB+ status (AAB+_status):
    - Additional donor filtering:
        - Only donors with AAB data
        - Non-diabetic donors **only**, regardless of AAB status
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + aab`

- AAB GADA+ (AAB_GADA):
    - Additional donor filtering:
        - Only donors with AAB data
        - Non-diabetic donors **only**
        - No donors positive for AAB *other than* GADA (must be AAB- or GADA+ only)
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + only.aab.gada`

- AAB IA2+ (AAB _IA2):
    - Additional donor filtering:
        - Only donors with AAB data
        - Non-diabetic donors **only**
        - No donors positive for AAB *other than* IA2 (must be AAB- or IA2+ only)
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + only.aab.ia.2`

- AAB IAA+ (AAB_IAA):
    - Additional donor filtering:
        - Only donors with AAB data
        - Non-diabetic donors **only**
        - No donors positive for AAB *other than* IAA (must be AAB- or IAA+ only)
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + only.aab.iaa`

- AAB ZNT8+ (AAB_ZNT8):
    - Additional donor filtering:
        - Only donors with AAB data
        - Non-diabetic donors **only**
        - No donors positive for AAB *other than* ZNT8 (must be AAB- or ZNT8+ only)
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + only.aab.znt8`

- Multiple AABs (Multi_AAB)
    - Additional donor filtering:
        - Only donors with AAB data
        - Non-diabetic donors **only**
        - Must be either AAB- or positive for $\geq$ 2 AABs
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + multi.aab`

- C. peptide (C_peptide):
    - Additional donor filtering:
        - Non-diabetic donors **only**
        - AAB- donors **only**
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + c.peptide`

- HbA1c (HbA1c):
    - Additional donor filtering:
        - Non-diabetic donors **only**
        - AAB- donors **only**
    - Initial DESeq formula:
        - `~ age + sex + bmi + ethnicity + chemistry + hb.a1c`



## Output files
- `<celltype>_<contrast_id>_results_initial_DESeq.tsv`: DESeq results from initial DESeq run with known co-variates (prior to RUVseq)
- `<celltype>_<contrast_id>_results_final_RUVseq.tsv`: DESeq results from final DESeq run, which uses known co-variates *and* latent variables from RUVseq
- `<celltype>_<contrast_id>_fGSEA_res_all.tsv`: Results from fGSEA for *all* pathways
- `<celltype>_<contrast_id>_fGSEA_res_signif.tsv`: Results from fGSEA for *significant* pathways (p < 0.1)
- `<celltype>_<contrast_id>_all_plots.pdf`: Concatenated PDF of all plots generated during analysis. Excludes plots that could not be generated due to pipeline stopping early.

## To run this pipeline yourself:

### Step 1:
**Important note:** metadata column names must be free of underscores `("_")`. Please change to periods `(".")` or eliminate if necessary before running this pipeline.

Modify [`pankbase_new_contrasts.csv`](src/2_run_DE/pankbase_new_contrasts.csv) to designate the contrasts you would like to run.
Here are descriptions of each column:

- `contrast_id`: what you would like to name the contrast. Generally used as output filename prefix.
- `contrast_var`: the name of the column you would like to perform the contrast on.
- `donor_col`: the name of the column with sample (donor) names.
- `control_grp`: if this is a binary contrast, whichever group is your control group (i.e. if expression is higher in control_grp, log2FoldChange will be negative).
- `experimental_grp`: if this is a binary contrast, whichever group is your experimental group (i.e. if expression is higher in experimetal group, log2FoldChange will be positive).
- `subset_on1`: would you like to subset data on a variable before running differential expression? Put that column name here.
- `subset_on2`: would you like to subset data on a second variable before running differential expression? Put that column name here.
- `subset1`: the groups you would like to *keep* from the column in `subset_on1`. The groups must be underscore-separated.
- `subset2`: the groups you would like to *keep* from the column in `subset_on2`. The groups must be underscore-separated.
- `covariates`: names of columns (covariates) you would like to include in the differential expression formula. If more than one, must be underscore-separated.
- `correlate_vars`: names of columns you would like to correlate with latent variables in a correlation matrix.
- `scale`: do you want to scale continuous variables? Must be TRUE or FALSE.
- `file_pattern`: the file ending of all of your pseudobulk matrix files. This should include **everything** other than the cell type.
- `assay`: must be "RNA" (no quotes).

### Step 2:

Modify [`submit_jobs_RUVseq_with_fGSEA.sh`](src/2_run_DE/submit_jobs_RUVseq_with_fGSEA.sh) to point to files and directories you want. 

**Note:** `pankbase_new_contrasts.csv` path must be changed at the beginning AND end of the file.

### Step 3:

Modify [`run_RUVseq_with_fGSEA.sh`](src/2_run_DE/run_RUVseq_with_fGSEA.sh) to point to your conda environment (with DESeq2 and RUVseq installed) and change all the SLURM info to your info.

### Step 4:

To run the pipeline:

`bash submit_jobs_RUVseq_with_fGSEA.sh`

This will submit a separate job for each contrast (row in your input csv).
