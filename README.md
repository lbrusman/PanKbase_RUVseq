# PanKbase RUVseq
This repo houses the modular RUVseq pipeline for PanKbase.

Here are some additional details about the analysis, which was performed on a per-cell type basis.

## The most recent version (v2) of the code found [here](src/2_run_DE/v2) runs my updated RUVseq differential expression pipeline.
## Change log for v2:
- Updated the RUVseq R script to pull pseudobulk matrices directly from the [PanKbase Data Library](https://data.pankbase.org) so that the pipeline no longer requires Seurat or the .rds object to run. This improves runtime and decreases memory requirements.
- Updated the bash and R scripts to pull donor- and biosample-level metadata directly from the PanKbase Data Library.
    - Link to Data Library [Donor Metadata](https://data.pankbase.org/tabular-files/PKBFI3142YFIU/) used.
    - Link to Data Library [Biosample Metadata](https://data.pankbase.org/tabular-files/PKBFI5333IJJK/) used.
- Created new script [make_sc_metadata.R](src/1_make_inputs/make_sc_metadata.R) to pull chemistry metadata and get cell counts from the single cell object.
    - This uses the single cell object as a **local** file. If you need to download the .rds object, you can get it [here](https://zenodo.org/records/15596314).
- Updated the method to filter the genes using a new `fancy_filter()` function.
- Got rid of redundant code and optimized more resource-intensive steps. This *drastically* improves runtime (from ~45 mins/contrast to ~12 mins/contrast).
- Updated the method used to optimize the number of latent variables (k). Briefly, the new method works as follows:
    - Find background non-DE genes by running an initial DESeq2 with variable of interest and co-variates. Background genes are genes with uncorrected p-value > 0.5
    - Use RUVg to normalize count data separately for every 1:max_k.
    - Calculate the relative log expression (RLE) per-sample from the normalized counts returned by RUVseq. More details about how RUVseq normalizes the data can be found [here](https://www.bioconductor.org/packages//release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html).
    - Run an ANOVA to determine how well each k normalizes the data (we want lowest F-statistic).
    - Plot k vs. F-statistic - the goal is to find the elbow of the plot (first k that minimizes F-statistic before diminishing returns)
    - In this version of the pipeline, I found the elbow of the plot by finding the point with the maximum second derivative (positive change in slope, or when the negative slope flattens out).
    - Remove any 1:best_k latent variables that are correlated with the variable of interest.
    - For final results, re-run DESeq2 with initial co-variates *and* final latent variables.
- The pipeline now requires two new input files.
    - [chemistry_metadata.csv](src/2_run_DE/v2/metadata_files/chemistry_metadata.csv) contains the 10x chemistry used for each sample.
    - [sc_object_cellcounts.csv](src/2_run_DE/v2/metadata_files/sc_object_cellcounts.csv) contains the number of cells in each cell type per-sample.
    - These should eventually be added to the Data Library and pulled directly from there.
- The parameter input .csv [pankbase_new_contrasts.csv](src/2_run_DE/v2/pankbase_new_contrasts.csv) has slightly different columns.
- I have included an [environment.yml](environment.yml) file. You can re-create my environment with:
    `conda env create -f environment.yml`. It may be missing some R packages - if you get an error about a package loading, make sure to check/install!

## Here is how I ran the contrasts that I have run so far using v2. Code for this version of the pipeline can be found [here](src/2_run_DE/v2).

### General filtering
Donors were filtered based on the following criteria:
1. Must have $\geq$ 20 cells in that cell type
2. Sample must not be treated with any reagents (`treatment == "no_treatment"`)

Genes were filtered based on the following criteria (`fancy_filter()` method):
1. Filter based on expression per-donor:
    - If the contrast variable is *categorical*, the gene must have $\geq$ 5 counts in at least half of the samples in the control group and/or half of the samples in experimental group.
    - If the contrast variable is *continuous*, the gene must have $\geq$ 5 counts in at least half of all samples.

For a cell type/contrast to be attempted: 
1. There must be $\geq$ 3 donors total left after filtering
2. If contrast variable is categorical, there must be $\geq$ 3 samples left in each group after filtering
3. If the contrast variable is categorical, there must > 1 unique group after filtering

Reasons pipeline could fail/exit early:
1. Not enough samples after filtering for the number of co-variates
2. Initial DESeq failed (error message is returned)
3. max_k = 0
4. RUVseq failed (error message is returned)

### Contrast-specific parameters
All parameters are summarized in [pankbase_new_contrasts.csv](src/2_run_DE/v2/pankbase_new_contrasts.csv). So far I have only used this pipeline on T1D and Age (to test a categorical and continuous variable).
- Control Without Diabetes vs. Type 1 Diabetes (ND_vs_T1D): 
    - Additional donor filtering:
        - None
    - Initial DESeq formula:
        - `~ Age..years. + Gender + BMI + Ethnicities + chemistry + Description.of.diabetes.status`
- Age (years)
    - Additional donor filtering:
        - Control Without Diabetes only
    - Initial DESeq formula:
        - `~ Gender + BMI + Ethnicities + chemistry + Age..years.`

## To run this pipeline yourself:

### Step 1:
**Important note:** metadata column names must be free of underscores `("_")`. Please change to periods `(".")` or eliminate if necessary before running this pipeline.

Modify [`pankbase_new_contrasts.csv`](src/2_run_DE/v2/pankbase_new_contrasts.csv) to designate the contrasts you would like to run.
Here are descriptions of each column:

- `contrast_id`: what you would like to name the contrast. Generally used as output filename prefix.
- `contrast_var`: the name of the column you would like to perform the contrast on.
- `sample_col`: the name of the column with sample names.
- `control_grp`: if this is a binary contrast, whichever group is your control group (i.e. if expression is higher in control_grp, log2FoldChange will be negative).
- `experimental_grp`: if this is a binary contrast, whichever group is your experimental group (i.e. if expression is higher in experimetal group, log2FoldChange will be positive).
- `subset_on1`: would you like to subset data on a variable before running differential expression? Put that column name here.
- `subset_on2`: would you like to subset data on a second variable before running differential expression? Put that column name here.
- `subset1`: the groups you would like to *keep* from the column in `subset_on1`. The groups must be underscore-separated.
- `subset2`: the groups you would like to *keep* from the column in `subset_on2`. The groups must be underscore-separated.
- `covariates`: names of columns (covariates) you would like to include in the differential expression formula. If more than one, must be underscore-separated.
- `correlate_vars`: names of columns you would like to correlate with latent variables in a correlation matrix.
- `scale`: do you want to scale continuous variables? Must be TRUE or FALSE.
- `cell_types`: cell types you want to perform contrasts on.
- `file_pattern`: the file ending of all of your pseudobulk matrix files. This should include **everything** other than the cell type.

### Step 2:

Modify [`submit_jobs_RUVseq_with_fGSEA.sh`](src/2_run_DE/v2/submit_jobs_RUVseq_with_fGSEA.sh) to point to files and directories you want. This script now pulls the donor- and biosample-level metadata directly from PanKbase (linked above), so if you want to change these files, this is where it should be done.

**Note:** `pankbase_new_contrasts.csv` path must be changed at the beginning AND end of the file.

### Step 3:

Modify [`run_RUVseq_with_fGSEA.sh`](src/2_run_DE/v2/run_RUVseq_with_fGSEA.sh) to point to your conda environment (can be re-created using provided `environment.yml`) and change all the SLURM info to your info.

### Step 4:

To run the pipeline:

`bash submit_jobs_RUVseq_with_fGSEA.sh`

This will submit a separate job for each contrast (row in your input csv).



## Here is how I ran the contrasts for v1 (results uploaded to PanKbase Data Library Dec 2025). Code for this version of the pipeline can be found [here](src/2_run_DE/v1).
<details>
<summary>Details here in drop-down </summary>


### General filtering
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


### Contrast-specific parameters
All parameters are summarized in [pankbase_new_contrasts.csv](src/2_run_DE/v1/pankbase_new_contrasts.csv)
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


## To run this pipeline yourself:

### Step 1:
**Important note:** metadata column names must be free of underscores `("_")`. Please change to periods `(".")` or eliminate if necessary before running this pipeline.

Modify [`pankbase_new_contrasts.csv`](src/2_run_DE/v1/pankbase_new_contrasts.csv) to designate the contrasts you would like to run.
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

Modify [`submit_jobs_RUVseq_with_fGSEA.sh`](src/2_run_DE/v1/submit_jobs_RUVseq_with_fGSEA.sh) to point to files and directories you want. 

**Note:** `pankbase_new_contrasts.csv` path must be changed at the beginning AND end of the file.

### Step 3:

Modify [`run_RUVseq_with_fGSEA.sh`](src/2_run_DE/v1/run_RUVseq_with_fGSEA.sh) to point to your conda environment (with DESeq2 and RUVseq installed) and change all the SLURM info to your info.

### Step 4:

To run the pipeline:

`bash submit_jobs_RUVseq_with_fGSEA.sh`

This will submit a separate job for each contrast (row in your input csv).

</details>

## Output files
- `<celltype>_<contrast_id>_meta_for_contrast.csv`: Final metadata for samples used in the contrast
- `<celltype>_<contrast_id>_results_initial_DESeq.tsv`: DESeq results from initial DESeq run with known co-variates (prior to RUVseq)
- `<celltype>_<contrast_id>_results_final_RUVseq.tsv`: DESeq results from final DESeq run, which uses known co-variates *and* latent variables from RUVseq
- `<celltype>_<contrast_id>_fGSEA_res_all.tsv`: Results from fGSEA for *all* pathways
- `<celltype>_<contrast_id>_fGSEA_res_signif.tsv`: Results from fGSEA for *significant* pathways (p < 0.1)
- `<celltype>_<contrast_id>_all_plots.pdf`: Concatenated PDF of all plots generated during analysis. May be an empty/corrupted pdf if contrast could not be attempted.
- `All_DE_stats_<contrast_id>.csv`: Stats and info for every the contrast in each cell type (number of samples, number of DE features, error messages, etc.)
- **If you want to save the best RUV-normalized counts:** There is code to do this in [RUVseq_with_fGSEA.R](https://github.com/lbrusman/PanKbase_RUVseq/blob/main/src/2_run_DE/v2/RUVseq_with_fGSEA.R). It is currently commented out, but it's this part:

```
# To save best RUV-normalized counts
norm_cts <- celltype_ruvseq[[best_k]]$normCounts
fname <- paste0(outdir, cell.type, "_", contrast_id, "_RUV_norm_cts.tsv")
write.table(norm_cts, fname, sep = "\t", quote = FALSE, row.names = TRUE)
```

## If you need to re-generate the pseudobulk matrices or the chemistry/cell count metadata from the single cell object
- Re-run code in [src/1_make_inputs](src/1_make_inputs).


