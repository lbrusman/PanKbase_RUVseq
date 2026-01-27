# Load packages
suppressMessages({
    library(dplyr)
    library(tidyr)
    library(data.table)
    library(ggplot2)
    library(R.utils)
    library(tibble)
    library(stringr)
    library(ggcorrplot)
    library(viridis)
    library(readr)
    library(DESeq2)
    library(RUVSeq)
    library(fgsea)
                     
})


source("0_clean_utils.R")


# Get arguments from command line
args = commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1]

keys <- attachLocally(args)
cat("Command-line arguments attached to global environment:\n");
str(mget(keys, envir=globalenv()))

set.seed(123)

# Let's get donor metadata straight from PanKbase Data Library
metadata_url <- "https://pankbase-data-v1.s3.us-west-2.amazonaws.com/analysis_resources/meta_data_table/human-donor_update_251210_v3.tar.gz"

download.file(metadata_url, destfile = "tmp.tar.gz")
untar("tmp.tar.gz",files="human-donor_update_251210_v3/pankbase_human_donor.txt")
meta_donor <- read.table("human-donor_update_251210_v3/pankbase_human_donor.txt", sep = "\t", header = TRUE)

# Biosample metadata from Data Library
biosample_url <- "https://pankbase-data-v1.s3.us-west-2.amazonaws.com/analysis_resources/meta_data_table/biosamples.2_2025-11-06.txt"

meta_biosample <- fread(biosample_url)


# We still need cell counts, so use that (local) file too... Should be uploaded to Data Library
metadata <- read.csv(metadata)


# Get cell types from metadata file
cell_types <- str_split(cell_types, "_") %>% unlist()
cell_types <- cell_types[order(cell_types)]
print(cell_types)

# Get out DE covariates
covariates_list <- str_split(covariates, pattern = "_") %>% unlist()

# Initial DESeq formula
deseq_formula <- paste0("~", paste0(covariates_list, collapse = "+"), "+", contrast_var)
print(deseq_formula)

# Get variables for correlation with latent vars
correlate_vars <- str_split(correlate_vars, pattern = "_") %>% unlist()


# Get subsets out
subset1 <- str_split(subset1, pattern = "_") %>% unlist()
subset2 <- str_split(subset2, pattern = "_") %>% unlist()

# Set FDR Threshold for DE
fdr <- 0.05


# General filters
ncells <- 20 # min number of cells in cell type for sample to be considered
minreads <- 10 # min reads in cell type for a gene to be kept
minprop <- 0.25 # min percent cells in cell type gene must be expressed in to be kept
nlatent.base <- 10

#get latent vars
latent_vars <- paste0("W_", 1:nlatent.base)

# Get base URL for pseudobulk matrices
base_url <- "https://pankbase-data-v1.s3.us-west-2.amazonaws.com/processed_data/single_cell_RNA/scRNA_matrices/"


de_stats <- data.frame() # Create df to hold info about DE in ALL cell types for this contrast
for (cell.type in cell_types) {
    print(cell.type)

    # Set up pdf to hold all plots
    pdf(paste0(outdir, cell.type, "_", contrast_id, "_all_plots.pdf"))

    # Gotta rename some cell types to match dataframe with celltype counts
    if (cell.type == "Immune") {
        cell.name <- "Immune..Macrophages."
    }

    else {
        cell.name <- gsub("[^[:alnum:]]", ".", cell.type)
    }

    print(cell.name)

    # Filter to only samples with >= the cutoff number of cells
    keep_meta <- metadata[metadata[,cell.name] >= ncells,]
    print(paste0("n samps with >=20 cells = ", nrow(keep_meta)))

    # Change names with special characters to match cellcount file
    if (cell.type == "MUC5B+ Ductal") {
        cell.url <- "MUC5b+Ductal"
        cell.url <- URLencode(cell.url, reserved = TRUE)
    }

    else {
        cell.url <- gsub(" ", "", cell.type)
        cell.url <- URLencode(cell.url, reserved = TRUE)
    }

    print(cell.url)

    # Read in celltype pseudobulk matrix from pankbase
    celltype_url <- paste0(base_url, cell.url, file_pattern)
    raw_mat <- fread(celltype_url)

    samples <- colnames(raw_mat)[colnames(raw_mat) != "V1"]
    
    pankbase_ids <- t(raw_mat[1:4,]) %>% as.data.frame()
    colnames(pankbase_ids) <- pankbase_ids[1,]
    pankbase_ids <- pankbase_ids[-1,]
    pankbase_ids <- pankbase_ids %>% tibble::rownames_to_column(var = sample_col)


    # Merge with donor metadata
    meta_donor$donor_accession <- meta_donor$Accession
    ids_with_donor <- left_join(pankbase_ids, meta_donor, by = "donor_accession")

    # Merge with biosample metadata
    meta_biosample$biosample_accession <- meta_biosample$Accession
    ids_with_donor_final <- left_join(ids_with_donor, meta_biosample, by = "biosample_accession")

    ids_with_donor_final <- ids_with_donor_final %>% filter(samples %in% keep_meta$samples)

    keep_meta <- keep_meta %>% merge(ids_with_donor_final, on = sample_col)

    # Now filter for only no_treatment samples
    keep_meta <- keep_meta %>% filter(treatment == "no_treatment")

    
    

    # Select only one sample when a donor has more than one
    a <- data.frame(table(keep_meta$donor_accession))
    coldata <- keep_meta[keep_meta$donor_accession %in% a[a$Freq == 1, "Var1"],]
    for (i in a[a$Freq > 1, "Var1"]) {
        tmp <- keep_meta[keep_meta$donor_accession == i,]
        tmp <- tmp[sample(1:nrow(tmp), 1),]
        coldata <- rbind(coldata, tmp)
    }

    keep_meta <- coldata

    # Get rid of metadata rows in pseudobulk matrix now
    raw_mat <- raw_mat[-c(1:4),]
    rownames(raw_mat) <- raw_mat$V1
    keep <- apply(raw_mat, 1, basic_filter) # Apply basic filter to genes
    raw_mat <- raw_mat[keep, ]
    names(raw_mat)[names(raw_mat) == "V1"] <- "gene"
    raw_mat <- as.data.frame(raw_mat)
    rownames(raw_mat) <- raw_mat$gene

    colnames(raw_mat) <- gsub("\\+", "-", colnames(raw_mat))
    print("here is the setdiff to match sample names")
    setdiff(colnames(raw_mat), keep_meta$samples)


    # If not enough samples left, exit early
    if (nrow(keep_meta) < 3) {
        er <- "less than 3 total samples left"
        mini_df <- data.frame(celltype = cell.type,
                              contrast = contrast_id,
                              n_control_samps = n_control_samps,
                              n_experimental_samps = n_experimental_samps,
                              n_samples_total = nrow(keep_meta),
                              base_formula = "NA",
                              n_base_DE_feats = "NA",
                              ruv_formula = "NA",
                              n_ruv_DE_feats = "NA",
                              error_message = er
                             )
        de_stats <- rbind(de_stats, mini_df)

        print(er)
            
        dev.off()
        { next }
    }

    # Now do some sample subsetting
    if (control_grp != "NA") {
        print(control_grp)
        print(experimental_grp)
        keep_meta <- keep_meta[which(keep_meta[[contrast_var]] %in% c(control_grp, experimental_grp)), ]

    }

    if (subset_on1 != "NA") {
        print(subset1)
        keep_meta <- keep_meta[which(keep_meta[[subset_on1]] %in% subset1), ]

    }

    if (subset_on2 != "NA") {
        print(subset2)
        keep_meta <- keep_meta[which(keep_meta[[subset_on2]] %in% subset2), ]

    }

    rownames(keep_meta) <- keep_meta$samples

    stash_genes <- rownames(raw_mat)
    raw_mat <- raw_mat[,keep_meta$samples]
    raw_mat <- sapply(raw_mat, as.numeric)
    rownames(raw_mat) <- stash_genes

    # Now turn numeric vars back to numeric
    keep_meta$Age..years. <- as.numeric(keep_meta$Age..years.)
    keep_meta$BMI <- as.numeric(keep_meta$BMI)


    # Subset samples to ones relevant for contrast
    skip <- FALSE
    if (is.factor(keep_meta[,contrast_var]) || is.character(keep_meta[,contrast_var])) {
        n_control_samps = nrow(keep_meta[keep_meta[,contrast_var] == control_grp, ])
        n_experimental_samps = nrow(keep_meta[keep_meta[,contrast_var] == experimental_grp, ])

        # Stop if not enough samples left
        if (n_control_samps < 3 | n_experimental_samps < 3) {
            er <- "less than 3 experimental or control samples left"
            print(er)

            mini_df <- data.frame(celltype = cell.type,
                              contrast = contrast_id,
                              n_control_samps = n_control_samps,
                              n_experimental_samps = n_experimental_samps,
                              n_samples_total = nrow(keep_meta),
                              base_formula = "NA",
                              n_base_DE_feats = "NA",
                              ruv_formula = "NA",
                              n_ruv_DE_feats = "NA",
                              error_message = er
                             )
            de_stats <- rbind(de_stats, mini_df)

            skip <- TRUE
            dev.off()
            { next }
        }

    }

    else {
        n_control_samps <- "NA"
        n_experimental_samps <- "NA"
    }

    # Make sure known covariates that are factors are indeed factors
    for (i in c(contrast_var, covariates_list)) {
        if (is.character(keep_meta[, i])) {
            keep_meta[,i] <- as.factor(keep_meta[,i])
        }
        
    }


    # Scale Age and BMI if desired
    if (scale == "TRUE") {
        for (i in c(contrast_var, covariates_list)) {
            if (is.numeric(keep_meta[, i])) {
                keep_meta[,i] <- scale(keep_meta[,i])
            } 
        }
    }


    # Can write csv to see what samples are being used for this contrast
    write.csv(keep_meta, paste0(outdir, cell.type, "_", contrast_id, "_meta_for_contrast.csv"))


    # Create dds object
    base_design <- paste0("~", paste0(covariates_list, collapse = "+"), "+", contrast_var)
    print(base_design)

    dds <- create_dds_obj(raw_mat, keep_meta)
    design(dds) <- as.formula(base_design)
    dds <- DESeq(dds)
    message("done with initial DESeq")
                                    

    # Do initial DESeq
    skip <- FALSE
    tryCatch({
        if (is.factor(keep_meta[,contrast_var])) {
            print("var is factor")
            contrast_vec <- c(contrast_var, experimental_grp, control_grp)
        
            # Get results
            result <- DESeq2::results(dds, contrast = contrast_vec)
        }
        else {
            print("var is numeric")

            # Get results
            result <- DESeq2::results(dds, name = contrast_var)
        }
    }, error = function(e) {
        skip <<- TRUE

        er <<- paste("Initial DESeq failed with message:", e$message)
    })

    if (skip) {
        mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(keep_meta),
                          base_formula = base_design,
                          n_base_DE_feats = "NA",
                          ruv_formula = "NA",
                          n_ruv_DE_feats = "NA",
                          error_message = er
                         )
        de_stats <- rbind(de_stats, mini_df)
        
        print("Initial DESeq failed on this cell type")
        dev.off()
        { next }
    }

    # Plot initial DESeq results
    xs <- data.frame(result)
    print(head(xs))
    topxs <- tibble::rownames_to_column(xs[which(xs$padj < fdr), ], var = "geneid")
    n_degs_base <- nrow(xs[which(xs$padj < fdr),])

    to_save <- cbind(gene = rownames(xs), xs)
    to_save <- to_save[order(to_save$pvalue),]

    # Save initial DESeq results
    write.table(to_save, paste0(outdir, cell.type, "_", contrast_id, "_results_initial_DESeq.tsv"),
               sep = "\t", quote = FALSE, row.names = FALSE)

    xs <- xs[!is.na(xs$padj),]
    plot <- ggplot(xs, aes(log2FoldChange, -log10(pvalue))) +
            geom_point(aes(col = ifelse(padj < fdr, "Signif.", "N.S")), size = .5) +
            scale_color_manual(values = c("gray", "firebrick")) +
            labs(col = "", title = base_design, subtitle = paste("number DEGs:", n_degs_base)) +
            theme_bw() +
            theme(plot.title = element_text(size = 8))
    plot <- plot + ggrepel::geom_text_repel(data = topxs, aes(x = log2FoldChange, y = -log10(pvalue), label = geneid), size = 3)
    print(plot)

    
    # Now let's run RUVseq
    # Calculate max number of latent variables based on number of samples remaining
    # Max k is number of samples - number covariates in base model - 1 intercept
    nlatent <- min(nlatent.base, nrow(keep_meta)-length(covariates_list)-1) # calculate 
    print(paste("nlatent =", nlatent))
    if (nlatent <= 0) {
        er <- "max nlatent is zero"

        mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(keep_meta),
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
    if (is.factor(keep_meta[,contrast_var])) {
        contrast_vec <- c(contrast_var, experimental_grp, control_grp)
        skip <- FALSE
        tryCatch({
            celltype_ruvseq <- run_ruvseq(
            dds,
            design = base_design,
            contrast = contrast_vec,
            k = nlatent,
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
                              n_samples_total = nrow(keep_meta),
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
            k = nlatent,
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
                              n_samples_total = nrow(keep_meta),
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


    k <- length(celltype_ruvseq)
    tmp <- data.frame(celltype_ruvseq[[k]]$W) %>% tibble::rownames_to_column(sample_col)
    tmp <- keep_meta %>% merge(tmp, on = sample_col)


    # Make corr mat with all known covariates and latent vars
    latent_names <- paste0("W_", 1:nlatent)
    t_mat <- matrix(nrow = length(correlate_vars), ncol = nlatent)
    rownames(t_mat) <- correlate_vars
    colnames(t_mat) <- latent_names

    p_mat <- matrix(nrow = length(correlate_vars), ncol = nlatent)
    rownames(p_mat) <- correlate_vars
    colnames(p_mat) <- latent_names

    # Run correlations
    for (v in correlate_vars) {
          print(v)
          for (x in latent_names) {
              this_form <- as.formula(paste0(x, " ~ ", v))
              stats_out <- lm(this_form, data = tmp, na.action=na.omit)
              stats_p <- summary(stats_out)$coefficients[2,4]
              stats_t <- summary(stats_out)$coefficients[2,3]
  
              t_mat[v, x] <- stats_t
              p_mat[v, x] <- stats_p
              
              }
        }

      # Plot correlation matrix
      # Get p-value matrix - bonferroni corrected if we want
      p.df <- as.data.frame(p_mat)*nrow(p_mat)*ncol(p_mat)

      # We can save correlations to check
      write.csv(p.df, paste0(outdir, cell.type, "_", contrast_id, "_test_corr_df.csv"))
      
      # Get asteriks matrix based on p-values
      p.labs <- p.df %>% mutate_all(labs.function)
      
      # Reshaping asterisks matrix to match ggcorrplot data output
      p.labs$Var1 <- as.factor(rownames(p.labs))
      p.labs <- reshape2::melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")
      
      # Initial ggcorrplot
      cor.plot <- ggcorrplot(t_mat, hc.order = FALSE,
                            lab = TRUE,
                            title = paste0(cell.type, " ", contrast_id, " correlation between variables"), 
                            lab_size = 2.25) + #original plots used 2.75
                scale_fill_gradient2(low = "#3A86FF", 
                                    mid = "white", 
                                    high = "#FD2244",
                                    limits = c(-15, 15),
                                    midpoint = 0,
                                    oob = scales::squish) +
                labs(fill = "t-statistic",
                    caption = "Bonferroni corrected p-values: * p < 0.05, ** p < 0.01, *** p < 0.001")
      p.labs$in.df <- ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                            paste0(cor.plot[["data"]]$Var1, cor.plot[["data"]]$Var2))),
                                "No", "Yes")
      
      p.labs <- dplyr::select(filter(p.labs, in.df == "Yes"), -in.df)
      
      # Add asterisks to ggcorrplot
      cor.plot.labs <- cor.plot +
        geom_text(aes(x = p.labs$Var1,
                      y = p.labs$Var2),
                  label = p.labs$lab,
                  nudge_y = 0.25,
                  size = 5) +
        guides(fill = guide_legend(title = "t-statistic"))
      
      print(cor.plot.labs)

      # Get out which latent vars associated with known co-variates
      corr_wide <- tibble::rownames_to_column(p.df, "correlate_vars")
      corr_wide <- setDT(corr_wide)

      corr_long <- reshape2::melt(corr_wide, id.vars = "correlate_vars", measure.vars = latent_names, variable.name = "latent_vars", value.name = "p_adj_bonf")
      corr_long <- corr_long %>% as.data.frame()
      pval_filt <- corr_long %>% filter(p_adj_bonf < 0.05)
    

      # Let's find best k by minimizing between-sample variance
      # Here Liza found the elbow of the RLE plot by finding the max second derivative

      anova_res_all <- data.frame()
      for (i in 1:length(celltype_ruvseq)) {
          print(i)
          this_res <- celltype_ruvseq[[i]]
          x1 = this_res$normCounts
          this_pheno <- this_res$W %>% tibble::rownames_to_column(sample_col)
          y1 <- log(x1+1)
          medi <- apply(y1, 1, median)
          rle_mat <- apply(y1, 2, function(x) x - medi)

          rle_long <- t(rle_mat) %>% as.data.frame()
          gene_names <- colnames(rle_long)
          rle_long[,sample_col] <- rownames(rle_long)

          rle_long <- rle_long %>% merge(this_pheno, by = sample_col)
          rle_long <- rle_long %>% pivot_longer(cols = gene_names)

          # ANOVA: RLE ~ sample
          aov_form <- as.formula(paste("value ~", sample_col))
          fit <- aov(aov_form, data = rle_long)
          anova_summary <- summary(fit)[[1]]
          f_stat <- anova_summary[sample_col, "F value"]
          res_var <- anova_summary["Residuals", "Mean Sq"]

          # Store ANOVA results
          Anove_stat_tmp <- data.frame(
                                  celltype = cell.type,
                                  contrast = contrast_id,
                                  donors = nrow(keep_meta),
                                  k = paste0("k_", i),
                                  k_num = i,
                                  # formula = i,
                                  f_statistic = f_stat,
                                  residual_variance = res_var
          )
        anova_res_all <- rbind(anova_res_all, Anove_stat_tmp)
        
      }

      # # We can save this to check if we want
      # write.csv(anova_res_all, paste0(outdir, cell.type, "_", contrast_id, "_anova_res_all.csv"))

                         
      # To find best k
      diff_df <- anova_res_all[order(anova_res_all$k_num),]

      # Find by maximizing second derivative
      f_stat <- diff_df$f_statistic
      d1 <- diff(f_stat)
      d1 <- c(NA, d1)
      d2 <- diff(d1)
      d2 <- c(d2, NA)

      diff_df$d1 <- d1
      diff_df$d2 <- d2

      best_k <- diff_df$k_num[which.max(diff_df$d2)]

      # # We can save this to check if we want
      # write.csv(diff_df, paste0(outdir, cell.type, "_", contrast_id, "_diff_df.csv"))

      message("  - Found new BestK -> ", best_k)


      # Plot RLE plot
      k_levels <- paste0("k_", 1:nlatent)
      anova_res_all$k <- factor(anova_res_all$k, levels = k_levels)
      gg <- ggplot(anova_res_all, aes(x = k, y = f_statistic, group = 1)) +
            geom_line() +
    
            geom_point(
              aes(size = residual_variance, fill = residual_variance),
              shape = 21,
              color = "black"
            ) +
    
            scale_fill_viridis(
              name = "Residual variance",
              option = "D"
            ) +
            scale_size_continuous(
              name = "Residual variance"
            ) +
            theme_classic() +
            labs(
              x     = "Number of RUV factors (k)",
              y     = "\n F-statistic\n(ANOVA: RLE ~ sample)",
              title = paste0(cell.type, " ", contrast_id, "\nBest-k -> ", best_k)
            ) +
            theme(
              axis.title    = element_text(size = 14),
              axis.text     = element_text(size = 12),
              axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1)
            )

      print(gg)

      # Find disease associated latent vars (Ws)
      corr_lats <- pval_filt %>% filter(correlate_vars == contrast_var)
      disease_associated_Ws <- unique(corr_lats$latent_vars)

      message("  - association found with disease status: \n    ", paste0(disease_associated_Ws,
                                                collapse = "|"))

      # Final formula!
      all_ws <- 1:best_k            
      all_ws <- paste0("W_", all_ws)

      W_terms <- all_ws[!all_ws %in% disease_associated_Ws]

      best_formula <- paste0("~", paste0(c(covariates_list, W_terms), collapse = "+"), "+", contrast_var)
      message("  - Best formula ->  ", best_formula)
              
      # Get best coldata
      best_ruv <- anova_res_all$formula[which(anova_res_all$k_num == best_k)]
      best_coldata <- celltype_ruvseq[[best_k]]$W 
      best_coldata$samples <- rownames(best_coldata)
      best_coldata <- best_coldata %>% merge(keep_meta, on = sample_col)
      rownames(best_coldata) <- best_coldata$samples

      # write.csv(best_coldata, paste0(outdir, cell.type, "_", contrast_id, "_best_ruv_coldata.csv"))

      # Reorder raw mat colnames
      raw_mat <- raw_mat[,rownames(best_coldata)]


      # Create final dds object
      dds <- create_dds_obj(raw_mat, best_coldata)
      

      if (is.factor(keep_meta[,contrast_var])) {
            print("var is factor")
            contrast_vec <- c(contrast_var, experimental_grp, control_grp)
            dds <- run_differential(dds, design = best_formula, contrast = contrast_vec)
        
            # extract results for the specified FDR threshold
            result <- dds$res
        }
      else {
            print("var is numeric")

            dds <- run_differential_continuous(dds, design = best_formula, name = contrast_var)

            result <- dds$res
      }
                           
      # Plot final RUVSeq results
      xs <- data.frame(result)
      topxs <- tibble::rownames_to_column(xs[which(xs$padj < fdr), ], var = "geneid")
      n_degs_ruv <- nrow(xs[which(xs$padj < fdr),])
    
      to_save <- cbind(gene = rownames(xs), xs)
      to_save <- to_save[order(to_save$pvalue),]

        
      write.table(to_save, paste0(outdir, cell.type, "_", contrast_id, "_results_final_RUVSeq.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
    
      xs <- xs[!is.na(xs$padj),]
      plot <- ggplot(xs, aes(log2FoldChange, -log10(pvalue))) +
                geom_point(aes(col = ifelse(padj < fdr, "Signif.", "N.S")), size = .5) +
                scale_color_manual(values = c("gray", "firebrick")) +
                labs(col = "", title = best_formula, subtitle = paste("number DEGs:", n_degs_ruv)) +
                theme_bw() +
                theme(plot.title = element_text(size = 8))
      plot <- plot + ggrepel::geom_text_repel(data = topxs, aes(x = log2FoldChange, y = -log10(pvalue), label = geneid), size = 3)
      print(plot)

      # Add stats/info about this contrast to larger csv to save at the end (for all cell types)
      mini_df <- data.frame(celltype = cell.type,
                          contrast = contrast_id,
                          n_control_samps = n_control_samps,
                          n_experimental_samps = n_experimental_samps,
                          n_samples_total = nrow(keep_meta),
                          base_formula = base_design,
                          n_base_DE_feats = n_degs_base,
                          ruv_formula = best_formula,
                          n_ruv_DE_feats = n_degs_ruv,
                          error_message = "NA"
                         )
      de_stats <- rbind(de_stats, mini_df)


      # Run fGSEA
      # To run fGSEA in Google colab, you will have to get these files another way!
      KEGG_react <- gmtPathways('/tscc/projects/ps-gaultonlab/lebrusman/from_TSCC_251212/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/reactome_kegg.gmt.txt')

      rpl <- fread('/tscc/projects/ps-gaultonlab/lebrusman/from_TSCC_251212/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/rpl_file_gsea.csv', fill=TRUE, header=TRUE)
      rps <- fread('/tscc/projects/ps-gaultonlab/lebrusman/from_TSCC_251212/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/rps_file_gsea.csv', fill=TRUE, header=TRUE)
      mtr <- fread('/tscc/projects/ps-gaultonlab/lebrusman/from_TSCC_251212/Gaulton_lab/code/RUVseq_pankbase/pipeline_just_for_pankbase/GSEA_files/mts_file_gsea.csv', fill=TRUE, header=TRUE)

      ribo_proteins <- c(rpl$`Approved symbol`, rps$`Approved symbol`, mtr$`Approved symbol`)
      ribo_proteins <- ribo_proteins[which(ribo_proteins != 'Approved symbol')]
    
    
      res <- to_save
      res <- res[which(!rownames(res) %in% ribo_proteins),]
      res$rank <- (-log10(as.numeric(res$pvalue)))*res$log2FoldChange
      res <- data.frame("SYMBOL" = rownames(res),
                        "stat" = res$rank)
      res <- res[!grepl(pattern = "NA", x = res$SYMBOL),]
      res <- res[!grepl(pattern = "MT-", x = res$SYMBOL),]
      ranks <- deframe(res)
      ranks <- sort(ranks, decreasing=TRUE)
      deseq_df <- data.frame(gene = res$SYMBOL, GSEArank = res$stat)

    
      KEGG_react_fgseaRes <- fgseaMultilevel(pathways = KEGG_react,
                                            stats = ranks,
                                            eps = 0.0,
                                            minSize  = 0, 
                                            maxSize  = 1000)
      message("Number of total enriched terms ", cell.type,": ", nrow(KEGG_react_fgseaRes))
                           
      # Set FDR threshold for fGSEA
      FDR_tresh <- 0.10
      KEGG_react_fgseaRes.tresh <- KEGG_react_fgseaRes[KEGG_react_fgseaRes$padj < FDR_tresh,]
      message("Number of significant terms: ", cell.type,": ",nrow(KEGG_react_fgseaRes.tresh))

      res <- KEGG_react_fgseaRes[order(KEGG_react_fgseaRes$pval), ]
      res_sig <- res[res$padj < FDR_tresh,]

      # Save fGSEA results
      fwrite(res, file=paste0(outdir, cell.type, "_", contrast_id, "_fGSEA_res_all.tsv"), sep="\t", sep2=c("", " ", ""), quote = FALSE)
      fwrite(res_sig, file=paste0(outdir, cell.type, "_", contrast_id, "_fGSEA_res_signif.tsv"), sep="\t", sep2=c("", " ", ""), quote = FALSE)

      # Plot top fGSEA results
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


    # Close pdf
    dev.off()
}



# Write csv with all DE stats/info
write.csv(de_stats, paste0(outdir, "All_DE_stats_", contrast_id, ".csv"), row.names = FALSE)





