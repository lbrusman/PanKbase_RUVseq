suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(parallel))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(GenomicRanges))
suppressMessages(library(patchwork))
suppressMessages(library(biomaRt))
suppressMessages(library(reticulate))
suppressMessages(library(RUVSeq))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(friendlyeval))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(tidyHeatmap))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(fgsea))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))

source("0_utils_pankbase.R")


#get arguments from command line
args = commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]

keys <- attachLocally(args)
cat("Command-line arguments attached to global environment:\n");
# print(keys);
str(mget(keys, envir=globalenv()))

set.seed(123)

meta <- read.table(metadata, sep = "\t", header = TRUE)
colnames(meta) <- gsub("_", "", colnames(meta))

files <- list.files(indir, pattern = file_pattern, full.names = TRUE)
og_celltypes <- list.files(indir, pattern = file_pattern, full.names = FALSE, recursive = FALSE) %>% str_split_i(file_pattern, 1)

#for pankbase
celltypes_fix <- gsub("\\+", ".", og_celltypes)
celltypes_fix <- gsub("\\(", ".", celltypes_fix)
celltypes_fix <- gsub("\\)", ".", celltypes_fix)

#let's get list of n_cells per cell type
meta[,celltypes_fix] <- sapply(meta[,celltypes_fix], as.numeric)
all_n_cells <- colSums(meta[,celltypes_fix], na.rm = TRUE)
names(all_n_cells) <- celltypes_fix

covariates_list <- str_split(covariates, pattern = "_") %>% unlist()

deseq_formula <- paste0("~", paste0(covariates_list, collapse = "+"), "+", contrast_var)
print(deseq_formula)

#get variables for correlation with latent vars
correlate_vars <- str_split(correlate_vars, pattern = "_") %>% unlist()

# Setup additional variables
max.k.global = 10 # Roughly 1/2 samples is ok

#get RLE plot thresholds for RNA or ATAC
if (assay == "RNA") {
  f_tresh = 10
}
if (assay == "ATAC") {
  f_tresh = 100
}

#get latent vars
latent_vars <- paste0("W_", 1:max.k.global)

# Write here the covariants you want to test against the Ws and EXCLUDE from the final formula (Usually just disease)
covariants <- contrast_var

#get subsets out
subset1 <- str_split(subset1, pattern = "_") %>% unlist()
print(subset1)
subset2 <- str_split(subset2, pattern = "_") %>% unlist()

#do RUVseq for T1D vs. non-diabetic
deseq_stats = data.frame()
lm_res_all = data.frame()

fdr <- 0.05


data <- readRDS("/tscc/nfs/home/lebrusman/Gaulton_lab/data/ADA_object/250424_ADA_object_metadata_v3_3.rds")
print(data)

outdir <- outdir

fdr <- 0.05

for (s in c("SRR27326986", "SRR27326987", "SRR27326992", "SRR27326993",
            "SRR27326994", "SRR27326995", "SRR27326996", "SRR27326997")) {
    data@meta.data[data@meta.data$samples == s, "samples"] <- paste0(data@meta.data[data@meta.data$samples == s, "samples"],
                                                                    "__", data@meta.data[data@meta.data$samples == s, "treatments"])
}

data@meta.data$coarse_annot <- gsub(" ", "", data@meta.data$coarse_annot)
Idents(data) <- data@meta.data$coarse_annot
unique_cell_types <- unique(data$coarse_annot)
samples <- unique(data@meta.data$samples)

md <- data@meta.data

meta_in_sc <- md %>% group_by(samples) %>% summarise(mean_nCount_RNA = mean(nCount_RNA), mean_nFeature_RNA = mean(nFeature_RNA),
                                                  mean_rna_total_reads = mean(rna_total_reads), mean_rna_uniquely_mapped_reads = mean(rna_uniquely_mapped_reads),
                                                  mean_rna_secondary_alignments = mean(rna_secondary_alignments), mean_rna_supplementary_alignments = mean(rna_supplementary_alignments),
                                                  mean_cell_probability = mean(cell_probability), mean_post_cellbender_umis = mean(post_cellbender_umis),
                                                  mean_pct_cellbender_removed = mean(pct_cellbender_removed), mean_rna_pct_mitochondrial = mean(rna_pct_mitochondrial))
meta_in_sc <- right_join(meta_in_sc, distinct(md[, c('samples', 'study', 'rrid', 'treatments', 'chemistry',
                                                     'sex', 'age', 'bmi', 'diabetes_status_description',
                                                     'source', 'aab_gada', 'aab_ia_2', 'aab_iaa', 'aab_znt8')]))

add_meta <- read.table("/tscc/nfs/home/lebrusman/Gaulton_lab/data/ADA_object/metadata_for_ADA.csv", sep = ",", header = T, fill = T) #from Liza and Parul
add_meta <- distinct(add_meta[, c("rrid", grep("PanKbase", colnames(add_meta), value = T))])
colnames(add_meta) <- colnames(add_meta) %>% janitor::make_clean_names()
meta_in_sc <- inner_join(meta_in_sc, add_meta, by = c("rrid" = "rrid"))


# manually fix some metadata due to mismatching info across releases
meta_in_sc[meta_in_sc$rrid == "HP-21337-01", "age"] <- 26
meta_in_sc[meta_in_sc$rrid == "HP-22234-01", "age"] <- 54

meta_in_sc[meta_in_sc$rrid == "HP-21337-01", "sex"] <- "Male"
meta_in_sc[meta_in_sc$rrid == "HP-23135-01", "sex"] <- "Female"

meta_in_sc[meta_in_sc$rrid == "HP-21337-01", "bmi"] <- 30.08
meta_in_sc[meta_in_sc$rrid == "HP-22234-01", "bmi"] <- 27.86

meta_in_sc[meta_in_sc$rrid == "HP-21337-01", "pan_kbase_bmi"] <- 30.08
meta_in_sc[meta_in_sc$rrid == "HP-22234-01", "pan_kbase_bmi"] <- 27.86



# fix status to get rid of space character
meta_in_sc$pan_kbase_description_of_diabetes_status <- ifelse(meta_in_sc$pan_kbase_description_of_diabetes_status == "non-diabetic", "NonDiabetic",
                                                             meta_in_sc$pan_kbase_description_of_diabetes_status)
meta_in_sc$pan_kbase_description_of_diabetes_status <- ifelse(meta_in_sc$pan_kbase_description_of_diabetes_status == "type 1 diabetes", "T1DM",
                                                             meta_in_sc$pan_kbase_description_of_diabetes_status)
meta_in_sc$pan_kbase_description_of_diabetes_status <- ifelse(meta_in_sc$pan_kbase_description_of_diabetes_status == "type 2 diabetes", "T2DM",
                                                             meta_in_sc$pan_kbase_description_of_diabetes_status)



#fix ethnicity column where some donors have missing info
meta_in_sc$pan_kbase_ethnicities <- ifelse(meta_in_sc$pan_kbase_ethnicities == "", "na", meta_in_sc$pan_kbase_ethnicities)
meta_in_sc[meta_in_sc$rrid == "HP-21077-01", "pan_kbase_ethnicities"] <- "Asian"
meta_in_sc[meta_in_sc$rrid == "HP-23135-01", "pan_kbase_ethnicities"] <- "Hispanic"


#fix isolation center column where some donors have missing info
meta_in_sc$pan_kbase_biosample_isolation_center <- ifelse(meta_in_sc$pan_kbase_biosample_isolation_center == "",
                                                         "NA", meta_in_sc$pan_kbase_biosample_isolation_center)


# if studies have different BMI than PanKbase, use PanKbase
meta_in_sc[which(meta_in_sc$bmi != meta_in_sc$pan_kbase_bmi), "bmi"] <-
    meta_in_sc[which(meta_in_sc$bmi != meta_in_sc$pan_kbase_bmi), "pan_kbase_bmi"]

# align AAB status between HPAP website and PanKbase
meta_in_sc$aab_gada <- ifelse(meta_in_sc$aab_gada == "na", "NA", meta_in_sc$aab_gada)
meta_in_sc$aab_gada <- ifelse(meta_in_sc$aab_gada == "1", "TRUE", meta_in_sc$aab_gada)
meta_in_sc$aab_gada <- ifelse(meta_in_sc$aab_gada == "0", "FALSE", meta_in_sc$aab_gada)

meta_in_sc$aab_ia_2 <- ifelse(meta_in_sc$aab_ia_2 == "na", "NA", meta_in_sc$aab_ia_2)
meta_in_sc$aab_ia_2 <- ifelse(meta_in_sc$aab_ia_2 == "1", "TRUE", meta_in_sc$aab_ia_2)
meta_in_sc$aab_ia_2 <- ifelse(meta_in_sc$aab_ia_2 == "0", "FALSE", meta_in_sc$aab_ia_2)

meta_in_sc$aab_znt8 <- ifelse(meta_in_sc$aab_znt8 == "na", "NA", meta_in_sc$aab_znt8)
meta_in_sc$aab_znt8 <- ifelse(meta_in_sc$aab_znt8 == "1", "TRUE", meta_in_sc$aab_znt8)
meta_in_sc$aab_znt8 <- ifelse(meta_in_sc$aab_znt8 == "0", "FALSE", meta_in_sc$aab_znt8)

meta_in_sc$aab_iaa <- ifelse(meta_in_sc$aab_iaa == "na", "NA", meta_in_sc$aab_iaa)
meta_in_sc$aab_iaa <- ifelse(meta_in_sc$aab_iaa == "1", "TRUE", meta_in_sc$aab_iaa)
meta_in_sc$aab_iaa <- ifelse(meta_in_sc$aab_iaa == "0", "FALSE", meta_in_sc$aab_iaa)

meta_in_sc$pan_kbase_aab_gada_positive <- ifelse(is.na(meta_in_sc$pan_kbase_aab_gada_positive), "NA", meta_in_sc$pan_kbase_aab_gada_positive)
meta_in_sc$pan_kbase_aab_ia2_positive <- ifelse(is.na(meta_in_sc$pan_kbase_aab_ia2_positive), "NA", meta_in_sc$pan_kbase_aab_ia2_positive)
meta_in_sc$pan_kbase_aab_znt8_positive <- ifelse(is.na(meta_in_sc$pan_kbase_aab_znt8_positive), "NA", meta_in_sc$pan_kbase_aab_znt8_positive)
meta_in_sc$pan_kbase_aab_iaa_positive <- ifelse(is.na(meta_in_sc$pan_kbase_aab_iaa_positive), "NA", meta_in_sc$pan_kbase_aab_iaa_positive)

col <- c('samples', 'mean_nCount_RNA', 'mean_nFeature_RNA', 'mean_rna_total_reads', 'mean_rna_uniquely_mapped_reads', 
         'mean_rna_secondary_alignments', 'mean_rna_supplementary_alignments', 
         'mean_cell_probability', 'mean_post_cellbender_umis', 
         'mean_pct_cellbender_removed', 'mean_rna_pct_mitochondrial',
         'study', 'rrid', 'treatments', 'chemistry',
         'sex', 'age', 'bmi', 'diabetes_status_description',
         'source', 'pan_kbase_aab_gada_positive', 'pan_kbase_aab_ia2_positive', 
         'pan_kbase_aab_znt8_positive', 'pan_kbase_aab_iaa_positive', 'pan_kbase_ethnicities', 'pan_kbase_biosample_isolation_center',
          'pan_kbase_c_peptide_ng_ml', 'pan_kbase_hb_a1c_percentage')



meta_in_sc <- meta_in_sc[, col]

colnames(meta_in_sc) <- c('samples', 'mean_nCount_RNA', 'mean_nFeature_RNA', 'mean_rna_total_reads', 
                          'mean_rna_uniquely_mapped_reads', 
                          'mean_rna_secondary_alignments', 'mean_rna_supplementary_alignments', 
                          'mean_cell_probability', 'mean_post_cellbender_umis', 
                          'mean_pct_cellbender_removed', 'mean_rna_pct_mitochondrial',
                          'study', 'rrid', 'treatments', 'chemistry',
                          'sex', 'age', 'bmi', 'diabetes_status_description',
                          'tissue_source', 'aab_gada', 'aab_ia_2', 'aab_znt8', 'aab_iaa', 
                          'ethnicity', 'isolation_center',
                          'c_peptide', 'hb_a1c')

print(head(meta_in_sc))

meta_in_sc$isolation_center <- gsub("  ", "", meta_in_sc$isolation_center)
meta_in_sc$isolation_center <- gsub(" ", "", meta_in_sc$isolation_center)

#added by liza
meta_in_sc$number.aab <- rowSums(meta_in_sc[,c("aab_gada", "aab_ia_2", "aab_iaa", "aab_znt8")] == "TRUE", na.rm = FALSE) #changed from na.rm = FALSE
meta_in_sc$number.aab <- ifelse(meta_in_sc$aab_gada %in% c("TRUE", "FALSE"), meta_in_sc$number.aab, NA)

meta_in_sc$multi.aab <- ifelse(meta_in_sc$number.aab >= 2 & meta_in_sc$aab_gada %in% c("TRUE", "FALSE"), "MultiAAB", 
                              ifelse(meta_in_sc$number.aab == 1 & meta_in_sc$aab_gada %in% c("TRUE", "FALSE"), "1AAB",
                                    ifelse(meta_in_sc$number.aab == 0 & meta_in_sc$aab_gada %in% c("TRUE", "FALSE"), "ZeroAAB", NA)))

meta_in_sc$only.aab.gada <- ifelse(meta_in_sc$aab_gada == "TRUE" & meta_in_sc$aab_ia_2 == "FALSE" &
                                  meta_in_sc$aab_iaa == "FALSE" & meta_in_sc$aab_znt8 == "FALSE", "TRUE", "FALSE")
meta_in_sc$only.aab.ia.2 <- ifelse(meta_in_sc$aab_gada == "FALSE" & meta_in_sc$aab_ia_2 == "TRUE" &
                                  meta_in_sc$aab_iaa == "FALSE" & meta_in_sc$aab_znt8 == "FALSE", "TRUE", "FALSE")
meta_in_sc$only.aab.iaa <- ifelse(meta_in_sc$aab_gada == "FALSE" & meta_in_sc$aab_ia_2 == "FALSE" &
                                  meta_in_sc$aab_iaa == "TRUE" & meta_in_sc$aab_znt8 == "FALSE", "TRUE", "FALSE")
meta_in_sc$only.aab.znt8 <- ifelse(meta_in_sc$aab_gada == "FALSE" & meta_in_sc$aab_ia_2 == "FALSE" &
                                  meta_in_sc$aab_iaa == "FALSE" & meta_in_sc$aab_znt8 == "TRUE", "TRUE", "FALSE")


#get new column with all aab info
meta_in_sc$all.aab.info <- meta_in_sc$multi.aab
meta_in_sc$all.aab.info <- ifelse(meta_in_sc$only.aab.gada == "TRUE", "OnlyGADA",
                                           ifelse(meta_in_sc$only.aab.ia.2 == "TRUE", "OnlyIA2",
                                                 ifelse(meta_in_sc$only.aab.iaa == "TRUE", "OnlyIAA",
                                                       ifelse(meta_in_sc$only.aab.znt8 == "TRUE", "OnlyZNT8", meta_in_sc$all.aab.info))))



# for (cell.type in c("Acinar", "ActiveStellate", "CyclingAlpha", "Delta", "Endothelial", "Gamma+Epsilon", 
#                     "Immune(Macrophages)", "MUC5B+Ductal", "QuiescentStellate", "Alpha", "Beta", "Ductal")) {
#     cell.prop <- data.frame(table(data@meta.data[data@meta.data$coarse_annot == cell.type, "samples"]))
#     write.table(cell.prop, paste0(outdir, "cell_props/", "cell.prop_", cell.type, ".txt"),
#                 quote = F, sep = "\t", row.names = F)


# }

if (contrast_id == "ND_vs_T1D") {
    write.table(meta_in_sc, paste0(outdir, "meta_in_sc_251118.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

}


metadata <- meta_in_sc
metadata <- metadata[metadata$treatments == "no_treatment",]


# Now let's do the DE analysis - loop through all celltypes
cell.types <- c("Beta", "ActiveStellate", "CyclingAlpha", "Delta", "Endothelial", 
            "Gamma+Epsilon", "Immune(Macrophages)", "MUC5B+Ductal", "QuiescentStellate", 
            "Alpha", "Acinar", "Ductal")
cell.types <- cell.types[order(cell.types)]
print(cell.types)
ncells <- 20
minreads <- 10
minprop <- 0.25
nlatent.base <- 10

input_dir <- "/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/outputs/250903_outs/pseudobulk_counts/"


de_stats <- data.frame()
for (cell.type in cell.types) {
    print(cell.type)
    pdf(paste0(outdir, cell.type, "_", contrast_id, "_all_plots.pdf"))

    cell.prop <- read.table(paste0("/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/outputs/250903_outs/cell_props/cell.prop_", cell.type, ".txt"), header = T) # cell counts and proportions in large map
    metadata <- meta_in_sc
    metadata <- metadata[metadata$treatments == "no_treatment",]
    
    a <- data.frame(table(metadata$rrid))
    coldata <- metadata[metadata$rrid %in% a[a$Freq == 1, "Var1"],]
    for (i in a[a$Freq > 1, "Var1"]) {
        tmp <- metadata[metadata$rrid == i,]
        set.seed(1234)
        tmp <- tmp[sample(1:nrow(tmp), 1),]
        coldata <- rbind(coldata, tmp) #coldata has unique sample per donor; no donor has multiple samples
    }

    # filter samples based on the number of cells and autoantibody/diabetes status
    cell.prop <- inner_join(cell.prop, metadata, by = c("Var1" = "samples"))
    cell.prop$aab <- "0"
    cell.prop$aab <- ifelse(cell.prop$aab_gada == "TRUE", "AAB+", cell.prop$aab)
    cell.prop$aab <- ifelse(cell.prop$aab_ia_2 == "TRUE", "AAB+", cell.prop$aab)
    cell.prop$aab <- ifelse(cell.prop$aab_iaa == "TRUE", "AAB+", cell.prop$aab)
    cell.prop$aab <- ifelse(cell.prop$aab_znt8 == "TRUE", "AAB+", cell.prop$aab)
    cell.prop$diabetes_status_description <- ifelse(cell.prop$aab == "0" & cell.prop$diabetes_status_description != "T1DM",
                                                    cell.prop$diabetes_status_description,
                                                    ifelse(cell.prop$aab == "AAB+" & cell.prop$diabetes_status_description == "NonDiabetic",
                                                           cell.prop$aab, cell.prop$diabetes_status_description))
    tmp <- cell.prop[cell.prop$Freq > ncells,] #keep non-treated samples with > 20 cells
    print("here is n samps with >20 cells")
    print(nrow(tmp))
        
    colnames(tmp) <- gsub("_", ".", colnames(tmp))
    if (control_grp != "NA") {
        print(control_grp)
        print(experimental_grp)
        tmp <- tmp[which(tmp[[contrast_var]] %in% c(control_grp, experimental_grp)), ]
        print("done with this group subset")
        print(tmp[[contrast_var]])
        print(nrow(tmp))

    }

    if (subset_on1 != "NA") {
        print(subset1)
        tmp <- tmp[which(tmp[[subset_on1]] %in% subset1), ]
        print("done with this subset1")
        print(tmp[[subset_on1]])
        print(nrow(tmp))

    }

    if (subset_on2 != "NA") {
        print(subset2)
        tmp <- tmp[which(tmp[[subset_on2]] %in% subset2), ]
        print("done with this subset2")
        print(tmp[[subset_on2]])
        print(nrow(tmp))

    }
    
    tmp <- cell.prop[cell.prop$Freq > ncells & 
              cell.prop$Var1 %in% tmp$Var1,]
    print(paste("nrow tmp:", nrow(tmp)))

    
    coldata <- coldata[coldata$samples %in% tmp$Var1,]

    coldata <- inner_join(coldata, tmp[,c("Var1", "aab")], by = c("samples" = "Var1"))
    coldata$aab <- as.factor(coldata$aab)
    
    colnames(coldata) <- gsub("_", ".", colnames(coldata)) #may have to find a better place for this

    coldata <- coldata[!is.na(coldata[[contrast_var]]), ] #get rid of na values if need to

    # if (nrow(coldata) < 3) {
    #     er <- "less than 3 samples left after filtering"

    #     mini_df <- data.frame(celltype = cell.type,
    #                       contrast = contrast_id,
    #                       n_control_samps = "NA",
    #                       n_experimental_samps = "NA",
    #                       n_samples_total = nrow(coldata),
    #                       base_formula = "NA",
    #                       n_base_DE_feats = "NA",
    #                       ruv_formula = "NA",
    #                       n_ruv_DE_feats = "NA",
    #                       error_message = er
    #                      )
    #     de_stats <- rbind(de_stats, mini_df)
        
    #     dev.off()
    #     { next }
        
    # }

    coldata <- data.frame(coldata)
    rownames(coldata) <- coldata$samples

    dir <- '/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/outputs/250903_outs/pseudobulk_counts/'
    raw_mat <- read.table(paste0(dir, cell.type, "_sample_gex_total_counts.txt"), header = T)

    coldata$samples <- gsub("-", ".", coldata$samples)
    tmp$Var1 <- gsub("-", ".", tmp$Var1)

    if (is.factor(coldata[,contrast_var]) || is.character(coldata[,contrast_var])) {
        n_control_samps = nrow(coldata[coldata[,contrast_var] == control_grp, ])
        n_experimental_samps = nrow(coldata[coldata[,contrast_var] == experimental_grp, ])
    }

    else {
        n_control_samps <- "NA"
        n_experimental_samps <- "NA"
    }
    
    
    skip <- FALSE
    tryCatch({
        raw_mat <- raw_mat[, intersect(coldata$samples, tmp$Var1)] # keep only samples with > `n_cells` cells
        raw_mat <- raw_mat[, which(colSums(raw_mat) > 0)] # remove samples that do not have any cells in the population
    }, error = function(e) {
        er <<- "not enough samples left after filtering"
        print(er)
        skip <<- TRUE
    })

    if (skip) {
        mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(coldata),
                          base_formula = "NA",
                          n_base_DE_feats = "NA",
                          ruv_formula = "NA",
                          n_ruv_DE_feats = "NA",
                          error_message = er
                         )
        de_stats <- rbind(de_stats, mini_df)
        
        dev.off()
        { next }
    }

    
    # print out the number of samples per diabetes status
    # print(table(tmp[tmp$Var1 %in% colnames(raw_mat), contrast_var]))
    
    # filter genes that have some min_reads raw counts in at least min_prop ratio of samples
    basic_filter <- function (row, min_reads = minreads, min_prop = minprop) {
      mean(row >= min_reads) >= min_prop
    }
    keep <- apply(raw_mat, 1, basic_filter) #1 means apply the function to each row
    raw_mat <- raw_mat[keep, ]

    rownames(coldata) <- coldata$samples
    coldata <- coldata[colnames(raw_mat),]
    
    print("check if all rownames(coldata) is in colnames(raw_mat):")
    print(all(rownames(coldata) %in% colnames(raw_mat)))

    # make sure known covariates that are factors are indeed factors
    for (i in c('study', 'samples', 'rrid', 'treatments', 'chemistry', 'sex', #'samples', 
                'diabetes.status.description', 'tissue.source', 'ethnicity',
                'aab.gada', 'aab.ia.2', 'aab.iaa', 'aab.znt8', 
                'number.aab', 'multi.aab', 'only.aab.gada', 'only.aab.ia.2', 'only.aab.iaa', 'only.aab.znt8', 'all.aab.info')) {
        coldata[, i] <- as.factor(coldata[, i])
        
    }


    colnames(coldata)[which(colnames(coldata) == "samples")] <- "sample.id"

    # #scale age and bmi
    if (scale == "TRUE") {
        coldata$age <- scale(coldata$age)
        coldata$bmi <- scale(coldata$bmi)
    }


    # normalize for library sizes
    size.factors <- DESeq2::estimateSizeFactorsForMatrix(raw_mat)
    celltype_norm_counts <- t(apply(raw_mat, 1, function(x) x/size.factors))
    celltype_norm_counts.log <- log2(celltype_norm_counts + 1) # add a pseudocount

    # options(repr.plot.width = 12, repr.plot.height = 6, repr.plot.res = 300)

    # #make plots to explore covariates
    # ord_idx <- match(colnames(celltype_norm_counts.log), coldata$sample.id)
    # group = contrast_var
    # par(mfrow=c(2,3))
    # colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
    # EDASeq::plotRLE(celltype_norm_counts, col = colors, outline = FALSE, las = 3, cex.axis = 1, ylab = "Relative Log Expression", main = NULL, cex.main = .5)
    # EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = F, main = group)
    # group = "tissue.source"
    # colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
    # EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = T, main = group)
    # group = "isolation.center"
    # colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
    # EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = T, main = group)
    # group = "sex"
    # colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
    # EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = F, main = group)
    # group = "chemistry"
    # colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
    # EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = F, main = group)

    if (is.factor(coldata[,contrast_var])) {
        n_control_samps = nrow(coldata[coldata[,contrast_var] == control_grp, ])
        n_experimental_samps = nrow(coldata[coldata[,contrast_var] == experimental_grp, ])

        if (n_control_samps <= 1 || n_experimental_samps <= 1) {
            er <- "one or fewer samples left in one or more groups"
            print(er)

            mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(coldata),
                          base_formula = "NA",
                          n_base_DE_feats = "NA",
                          ruv_formula = "NA",
                          n_ruv_DE_feats = "NA",
                          error_message = er
                         )
            de_stats <- rbind(de_stats, mini_df)
            
            dev.off()
            { next }
        }
    }

    else {
        n_control_samps <- "NA"
        n_experimental_samps <- "NA"
    }

    if (length(unique(coldata[,contrast_var])) < 2) {

        er <- "only one contrast group left after filtering"
        print(er)

        mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(coldata),
                          base_formula = "NA",
                          n_base_DE_feats = "NA",
                          ruv_formula = "NA",
                          n_ruv_DE_feats = "NA",
                          error_message = er
                         )
        de_stats <- rbind(de_stats, mini_df)
        
        dev.off()
        { next }
    }

    #let's move this here
    covariates_list <- str_split(covariates, pattern = "_") %>% unlist()
    for (i in covariates_list) {
        if (length(unique(coldata[,i])) < 2) {
            covariates_list <- covariates_list[covariates_list != i]
        }
    }

    #let's save what coldata is actually being used
    write.table(coldata, paste0(outdir, cell.type, "_", contrast_id, "_meta_in_sc_251118.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    # quit()


    # Create dds object
    dds <- create_dds_obj(raw_mat, coldata)

    #plot initial volcano plot - Liza
    base_design <- paste0("~", paste0(covariates_list, collapse = "+"), "+", contrast_var)

                                    
    # base_design <- "~ diabetes_status_description+sex+age+bmi+ethnicity"
    additional_covs <- c("chemistry")

    # Run this chunk on a SLURM job in case k is relatively large. When k = 30, it takes ~ 45 mins to run the analysis
    dds <- create_dds_obj(raw_mat, coldata)

    skip <- FALSE
    tryCatch({
        if (is.factor(coldata[,contrast_var])) {
            print("var is factor")
            contrast_vec <- c(contrast_var, experimental_grp, control_grp)
    
            celltype_de_explore <- run_many_designs_deseq(
              dds,
              base_design,
              additional_covs,
              contrast = contrast_vec,
              shrink = FALSE
            )
    
            # Run this chunk on a SLURM job in case k is relatively large. When k = 30, it takes ~ 45 mins to run the analysis
            dds <- create_dds_obj(raw_mat, coldata)
        
            #plot initial volcano plot - Liza
            design <- base_design
            design(dds) <- as.formula(design)
            dds <- DESeq(dds)
            # extract results for the specified FDR threshold
            result <- DESeq2::results(dds, contrast = contrast_vec)
        }
        # contrast_vec = c("diabetes_status_description", "T1DM", "NonDiabetic") # fold change = numerator / denominator
        else {
            print("var is numeric")
            celltype_de_explore <- run_many_designs_deseq_continuous(
              dds,
              base_design,
              additional_covs,
              name = contrast_var,
              shrink = FALSE
            )
            # Run this chunk on a SLURM job in case k is relatively large. When k = 30, it takes ~ 45 mins to run the analysis
            dds <- create_dds_obj(raw_mat, coldata)
        
            #plot initial volcano plot - Liza
            design <- base_design
            design(dds) <- as.formula(design)
            dds <- DESeq(dds)
            # extract results for the specified FDR threshold
            result <- DESeq2::results(dds, name = contrast_var)
        }
    }, error = function(e) {
        skip <<- TRUE

        er <<- paste("Initial DESeq failed with message:", e$message)
        print(er)
    })

    if (skip) {
        mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(coldata),
                          base_formula = base_design,
                          n_base_DE_feats = "NA",
                          ruv_formula = "NA",
                          n_ruv_DE_feats = "NA",
                          error_message = er
                         )
        de_stats <- rbind(de_stats, mini_df)
        
        print("Initial DESeq failed on this cell type")
        dev.off()
        {next}
    }
    

    xs <- data.frame(result)
    print(head(xs))
    topxs <- tibble::rownames_to_column(xs[which(xs$padj < fdr), ], var = "geneid")
    n_degs_base <- nrow(xs[which(xs$padj < fdr),])

    to_save <- cbind(gene = rownames(xs), xs)
    to_save <- to_save[order(to_save$pvalue),]

    write.table(to_save, paste0(outdir, cell.type, "_", contrast_id, "_results_initial_DESeq.tsv"),
               sep = "\t", quote = FALSE, row.names = FALSE)

    print("making volcano plot")
    xs <- xs[!is.na(xs$padj),]
    plot <- ggplot(xs, aes(log2FoldChange, -log10(pvalue))) +
            geom_point(aes(col = ifelse(padj < fdr, "Signif.", "N.S")), size = .5) +
            scale_color_manual(values = c("gray", "firebrick")) +
            labs(col = "", title = design, subtitle = paste("number DEGs:", n_degs_base)) +
            theme_bw() +
            theme(plot.title = element_text(size = 8))
    plot <- plot + ggrepel::geom_text_repel(data = topxs, aes(x = log2FoldChange, y = -log10(pvalue), label = geneid), size = 3)
    # options(repr.plot.width = 6, repr.plot.height = 5, repr.plot.res = 300)
    print(plot)

    cat_vars <- c("sex", "diabetes.status.description", "ethnicity", "chemistry")
    cont_vars <- c("age", "bmi")
    n_subtract <- 0
    for (i in cat_vars) {
        i_len <- length(unique(coldata[,i]))
        n_subtract <- n_subtract + i_len
    }
    n_subtract <- n_subtract + length(cont_vars) + 1
    # nlatent <- min(nlatent.base, nrow(coldata)-n_subtract)
    nlatent <- min(nlatent.base, nrow(coldata)-length(covariates_list)-4)
    print(paste("nlatent =", nlatent))
    if (nlatent <= 0) {
        er <- "max nlatent is zero"

        mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(coldata),
                          base_formula = base_design,
                          n_base_DE_feats = n_degs_base,
                          ruv_formula = "NA",
                          n_ruv_DE_feats = "NA",
                          error_message = er
                         )
        de_stats <- rbind(de_stats, mini_df)
        
        print(er)
        dev.off()
        { next }
    }
    if (is.factor(coldata[,contrast_var])) {
        contrast_vec <- c(contrast_var, experimental_grp, control_grp)
        skip <- FALSE
        tryCatch({
            celltype_ruvseq <- run_ruvseq(
            dds,
            design = base_design,
            contrast = contrast_vec,
            k = min(nlatent, 10), #k at most is number of samples - 6 covariates in base model - 1 intercept
            p.val.thresh = 0.5,
            method = "ruvg"
        )
        }, error = function(e) {
            skip <<- TRUE

            er <<- paste("RUVseq failed with message:", e$message)
            print(er)
          
        })
    
        if (skip) {
            mini_df <- data.frame(celltype = cell.type,
                              contrast = contrast_id,
                              n_control_samps = n_control_samps,
                              n_experimental_samps = n_experimental_samps,
                              n_samples_total = nrow(coldata),
                              base_formula = base_design,
                              n_base_DE_feats = n_degs_base,
                              ruv_formula = "NA",
                              n_ruv_DE_feats = "NA",
                              error_message = er
                             )
            de_stats <- rbind(de_stats, mini_df)
            
            print("RUVseq failed on this cell type")
            dev.off()
            {next}
        }

    }
    else {
        skip <- FALSE
        tryCatch({
            celltype_ruvseq <- run_ruvseq_continuous(
            dds,
            design = base_design,
            name = contrast_var,
            k = min(nlatent, 10), #k at most is number of samples - 6 covariates in base model - 1 intercept
            p.val.thresh = 0.5,
            method = "ruvg"
        )
        }, error = function(e) {
            skip <<- TRUE

            er <<- paste("RUVseq failed with message:", e$message)
            print(er)

            
        })
    
        if (skip) {
            mini_df <- data.frame(celltype = cell.type,
                              contrast = contrast_id,
                              n_control_samps = n_control_samps,
                              n_experimental_samps = n_experimental_samps,
                              n_samples_total = nrow(coldata),
                              base_formula = base_design,
                              n_base_DE_feats = n_degs_base,
                              ruv_formula = "NA",
                              n_ruv_DE_feats = "NA",
                              error_message = er
                             )
            de_stats <- rbind(de_stats, mini_df)
            
            print("RUVseq failed on this cell type")
            dev.off()
            {next}
        }
    }

    
    coldata <- data.frame(colData(celltype_de_explore[[base_design]]$dds))
    print(head(coldata))
    cell.prop$Var1 <- gsub("-", ".", cell.prop$Var1)
    tmp_df <- distinct(cell.prop[, c("Var1", "Freq")])
    tmp_df <- cell.prop[, c("Var1", "Freq")]
    print(head(tmp_df))
    colnames(tmp_df) <- c("sample.id", "cell_counts")
    coldata <- inner_join(coldata, tmp_df, by = c("sample.id" = "sample.id"))
    rownames(coldata) <- coldata$sample.id
    print(head(coldata))


    k <- length(celltype_ruvseq)
    tmp <- data.frame(celltype_ruvseq[[k]]$W) %>% tibble::rownames_to_column("sample.id")
    tmp <- combine_by_sampleid(coldata, tmp)
    print("final tmp:")
    print(head(tmp))
    tmp <- dplyr::select(tmp, -any_of(c("sample.id", "rrid", "samples"))) %>%
        DataExplorer::dummify() %>%
        dplyr::select_if(is_almost_ok, ~.x)
    a <- psych::corr.test(tmp[, grep("W_", colnames(tmp))], tmp[, grep("W_|PC", colnames(tmp), invert = T)],
                         use = "na.or.complete", method = "spearman", adjust = "BH")
    # print(a)
    
    # tryCatch({
    #     corrplot::corrplot(a$r, p.mat = a$p.adj, tl.cex = .75, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.8,
    #                      insig = "label_sig", diag = FALSE)
    # }, error = function(e) {
    #     print("could not plot corrplot. number of latent vars may = 1")
    # })
    

    ## calculate which k to search
    tmp <- rowSums(a$p.adj < 0.05)
    print("here is number of p.adj < 0.05")
    print(tmp)
    if (all(tmp == 0) || length(tmp) == 0) {
        er <- "no known vars correlate with latent vars?"

        mini_df <- data.frame(celltype = cell.type,
                              contrast = contrast_id,
                              n_control_samps = n_control_samps,
                              n_experimental_samps = n_experimental_samps,
                              n_samples_total = nrow(coldata),
                              base_formula = base_design,
                              n_base_DE_feats = n_degs_base,
                              ruv_formula = "NA",
                              n_ruv_DE_feats = "NA",
                              error_message = er
                             )
        de_stats <- rbind(de_stats, mini_df)
        
        print(er)
        dev.off()
        {next}
    }
    tmp <- data.frame(tmp)
    tmp$k <- 1:nrow(tmp)
    colnames(tmp) <- c("n_known_vars", "k")
    k_stop <- tmp[tmp$n_known_vars == 0, "k"][1] # get the smallest k that does not significantly correlate with any known variables

    if (is.na(k_stop)) {
        k_stop <- tmp[nrow(tmp), "k"]
    }

    print(paste0("k_stop = ", k_stop))
    
    # options(repr.plot.width = 6, repr.plot.height = 4, repr.plot.res = 300)
    #try to plot upset plot
    # tryCatch({
    #     p <- ggplot(tmp, aes(k, n_known_vars)) + geom_line() +
    #       geom_point() + ggrepel::geom_text_repel(aes(label = glue("{n_known_vars}")), max.overlaps = Inf, box.padding = 0.1) +
    #       scale_x_continuous(breaks=seq(0, length(celltype_ruvseq), by=1)) +
    #       labs(x = 'Num. of latent variables', y = 'Total n known vars significantly corr.',
    #            title = glue("{cell.type}_nCells-{20}_reads-{10}_in-{1/4}-nSamples_fdr-{0.05}")) +
    #       theme_bw()
    #     # print(p)
    # }, error = function(e) {
    #     print("couldn't make corr line plot plot. keep going...")
    # })


    #find which known vars. corr. with latent vars
    tmp <- data.frame(a$p.adj < 0.05)[1:(k_stop-1),]
    print(tmp)
    x <- vector(mode = "list", length = length(1:(k_stop-1)))
    for (k in 1:(k_stop-1)) {
        j <- which(1:(k_stop-1) == k)
        x[[j]] <- colnames(tmp)[which(tmp[k,] == TRUE)]
    }
    names(x) <- paste0("k_", 1:(k_stop-1))

                                    
    # # options(repr.plot.width = 5, repr.plot.height = 6, repr.plot.res = 300)
    tryCatch({
        p <- upset(fromList(x), order.by = "freq", nsets = 100,
        text.scale = c(1, 1, 1, 1, 1, 2)) #c(intersection size title, intersection size tick labels,
                                            #set size title, set size tick labels, set names, numbers above bars).
        # print(p)

    }, error = function(e) {
        print("couldn't make upset plot. keep going...")
    })


    #choose k
    # Function to find the smallest k with a unique element
    find_smallest_k_with_unique <- function(lists) {
      k_with_unique <- c()
      all_elements <- unlist(lists)
      unique_elements <- unique(all_elements)
      element_count <- table(all_elements)
      
      for (i in seq_along(lists)) {
        k_elements <- lists[[i]]
        # Find unique elements in the current list
        unique_in_k <- k_elements[k_elements %in% unique_elements[element_count[k_elements] == 1]]
        if (length(unique_in_k) > 0) {
          k_with_unique <- c(k_with_unique, i)
        }
      }
        if (length(k_with_unique) > 0) {
            return(max(k_with_unique))
        } else {
          return(NULL) # Return NULL if no unique element is found
            }
    }
    
    if ( length(grep(contrast_var, x)) > 0 ) {
        k_chosen <- grep(contrast_var, x) - 1
    } else {
        k_chosen <- find_smallest_k_with_unique(x)
    }
                                    
    tryCatch({
        k_chosen <- max(k_chosen)
        print(paste("k chosen:", k_chosen))
    })

    if (is.null(k_chosen) || k_chosen == 0 || k_chosen == Inf || k_chosen == -Inf) {
        er <- "k_chosen is null or zero"

        mini_df <- data.frame(celltype = cell.type,
                              contrast = contrast_id,
                              n_control_samps = n_control_samps,
                              n_experimental_samps = n_experimental_samps,
                              n_samples_total = nrow(coldata),
                              base_formula = base_design,
                              n_base_DE_feats = n_degs_base,
                              ruv_formula = "NA",
                              n_ruv_DE_feats = "NA",
                              error_message = er
                             )
        de_stats <- rbind(de_stats, mini_df)
        
        print(er)
        dev.off()
        {next}
    }


    k <- max(k_chosen)

    #to save norm counts
    norm_cts <- celltype_ruvseq[[k]]$normCounts
    fname <- paste0(outdir, cell.type, "_", contrast_id, "_RUV_norm_cts.tsv")
    write.table(norm_cts, fname, sep = "\t", quote = FALSE, row.names = TRUE)

                                    
    de_object <- celltype_ruvseq[[k]]$de
    print(de_object$design)
    fdr <- 0.05
    xs <- data.frame(de_object$res)
    print(head(xs))
    topxs <- tibble::rownames_to_column(xs[which(xs$padj < fdr), ], var = "geneid")
    n_degs_ruv <- nrow(xs[which(xs$padj < fdr),])
    t <- unlist(strsplit(x = de_object$design, "+", fixed = T))

    to_save <- cbind(gene = rownames(xs), xs)
    to_save <- to_save[order(to_save$pvalue),]

    write.table(to_save, paste0(outdir, cell.type, "_", contrast_id, "_results_final_RUVseq.tsv"),
               sep = "\t", quote = FALSE, row.names = FALSE)

    print("making volcano plot")
    xs <- xs[!is.na(xs$padj),]
    plot <- ggplot(xs, aes(log2FoldChange, -log10(pvalue))) +
            geom_point(aes(col = ifelse(padj < fdr, "Signif.", "N.S")), size = .5) +
            scale_color_manual(values = c("gray", "firebrick")) +
            labs(col = "", title = de_object$design, subtitle = paste("number DEGs:", n_degs_ruv)) +
            theme_bw() +
            theme(plot.title = element_text(size = 8))
    plot <- plot + ggrepel::geom_text_repel(data = topxs, aes(x = log2FoldChange, y = -log10(pvalue), label = geneid), size = 3)
    # options(repr.plot.width = 6, repr.plot.height = 5, repr.plot.res = 300)
    print(plot)

    # dev.off()

                                    
    mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(coldata),
                          base_formula = base_design,
                          n_base_DE_feats = n_degs_base,
                          ruv_formula = de_object$design,
                          n_ruv_DE_feats = n_degs_ruv,
                          error_message = "NA"
                         )
    de_stats <- rbind(de_stats, mini_df)




    #to run fGSEA
    ### Run GSEA on each cell type and disease
    KEGG_react <- gmtPathways('/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/reactome_kegg.gmt.txt')

    mainDir <- outdir
    # dir.create(file.path(mainDir))

    rpl <- fread('/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/rpl_file_gsea.csv', fill=TRUE, header=TRUE)
    rps <- fread('/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/rps_file_gsea.csv', fill=TRUE, header=TRUE)
    mtr <- fread('/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/mts_file_gsea.csv', fill=TRUE, header=TRUE)

    ribo_proteins <- c(rpl$`Approved symbol`, rps$`Approved symbol`, mtr$`Approved symbol`)
    ribo_proteins <- ribo_proteins[which(ribo_proteins != 'Approved symbol')]
    
    
   # beta_t1d <- beta_t1d[which(!rownames(beta_t1d) %in% ribo_proteins),]
    res <- to_save
    res <- res[which(!rownames(res) %in% ribo_proteins),]
    res$rank = (-log10(as.numeric(res$pvalue)))*res$log2FoldChange
    res = data.frame("SYMBOL" = rownames(res),
                     "stat" = res$rank)
    res = res[!grepl(pattern = "NA", x = res$SYMBOL),]
    res = res[!grepl(pattern = "MT-", x = res$SYMBOL),]
    ranks <- deframe(res)
    ranks = sort(ranks, decreasing=TRUE)
    deseq_df <- data.frame(gene = res$SYMBOL, GSEArank = res$stat)
    head(ranks)

    
    KEGG_react_fgseaRes <- fgseaMultilevel(pathways=KEGG_react,
                                stats=ranks,
                                eps =0.0,
                                minSize  = 0, 
                                maxSize  = 1000)#,nPermSimple = 10000)
    message("Number of total enriched terms ", cell.type,": ", nrow(KEGG_react_fgseaRes))
    FDR_tresh = 0.10
    KEGG_react_fgseaRes.tresh = KEGG_react_fgseaRes[KEGG_react_fgseaRes$padj < FDR_tresh,]
    message("Number of significant terms: ", cell.type,": ",nrow(KEGG_react_fgseaRes.tresh))
    ## Add categories
    # res <- KEGG_react_fgseaRes.tresh

    res <- KEGG_react_fgseaRes[order(KEGG_react_fgseaRes$pval), ]
    res_sig <- res[res$padj < FDR_tresh,]
    
    fwrite(res, file=paste0(mainDir, cell.type, "_", contrast_id, "_fGSEA_res_all.tsv"), sep="\t", sep2=c("", " ", ""), quote = FALSE)
    fwrite(res_sig, file=paste0(mainDir, cell.type, "_", contrast_id, "_fGSEA_res_signif.tsv"), sep="\t", sep2=c("", " ", ""), quote = FALSE)

    # fgsea_res_ord <- res[order(res$NES, decreasing = TRUE),]
    # fgsea_res_top <- fgsea_res_ord %>% slice_head(n = 10)

    fgsea_res_ord <- res_sig[order(res_sig$NES, decreasing = TRUE),]
    fgsea_res_top <- fgsea_res_ord %>% slice_head(n = 10)
    fgsea_res_bottom <- fgsea_res_ord %>% slice_tail(n = 10)
    fgsea_tops <- rbind(fgsea_res_top, fgsea_res_bottom)

    nrow_sig <- nrow(res_sig)

    if (nrow_sig > 0) {
        top_n <- min(20, nrow_sig)
        fgsea_top_sig <- res_sig[1:top_n,]

        p <- ggplot(fgsea_top_sig, aes(reorder(pathway, NES), NES)) +
              geom_col(aes(fill=padj)) +
              coord_flip() +
              labs(x="Pathway", y="Normalized Enrichment Score",
                   title="Top significant pathways (p < 0.1) from GSEA") + 
              theme_minimal() +
              theme(text = element_text(size = 6))
        print(p)

    }

    
    dev.off()

        
}
                                    
print(de_stats)
write.csv(de_stats, paste0(outdir, "All_DE_stats_", contrast_id, ".csv"), row.names = FALSE)





