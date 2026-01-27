# Helper functions

suppressPackageStartupMessages({
  require(magrittr)
  require(ggplot2)
  require(org.Hs.eg.db)
  require(glue)
})

# filter genes that have some min_reads raw counts in at least min_prop ratio of samples
basic_filter <- function (row, min_reads = minreads, min_prop = minprop) {
  mean(row >= min_reads) >= min_prop
}


#' DESeq2's Variance Stabilizing Transformation
normalize_deseq <- function(counts, ...) {
  # Outputs a "log-like" transformed value
  DESeq2::varianceStabilizingTransformation(counts, ...)
}


#' Create DDS object.
#' This object rarely changes, so do not compute it everytime.
create_dds_obj <- function(counts_mat, metadata) {
  sample.ids <- colnames(counts_mat)
  message(glue("info: counts matrix: {paste(dim(counts_mat), collapse = ', ')}"))
  message(glue("info: covariate matrix: {paste(dim(metadata), collapse = ', ')}"))

  col_data <- metadata


  if (any(rownames(col_data) != colnames(counts_mat))) {
    print("names don't match. check or error!")
    print(setdiff(rownames(col_data), colnames(counts_mat)))
    return()
  }

  message("info: creating dds object..")
  DESeq2::DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~1)
}

#' Run the differential expression logic for the given design formula
run_differential <- function(dds, design, contrast, interaction = NULL, fdr = 0.05, shrink = FALSE) {
  if (!is.null(interaction))
    design <- paste0(design, "+", interaction)

  writeLines(glue("info: running design {design}"))
  design(dds) <- as.formula(design)  # replace / update design

  dds <- DESeq2::DESeq(dds)
  vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))

  # extract results for the specified FDR threshold
  result <- DESeq2::results(dds, contrast = contrast, alpha = fdr)

  # if shrink is specifies, then report shrunken log2FC using apeglm method
  if (shrink) {
    writeLines("info: estimating shrunk log2FC")
    coef <- gsub(' ', '.', glue("{contrast[1]}_{contrast[2]}_vs_{contrast[3]}"))
    result <- DESeq2::lfcShrink(dds, coef = coef, type = "apeglm")
  }

  return(list(dds = dds, vst = vst, res = result, design = design, shrink = shrink))
}

# Added by Liza
run_differential_continuous <- function(dds, design, name, fdr = 0.05, shrink = FALSE) {
    design(dds) <- as.formula(design)
    dds <- DESeq2::DESeq(dds)
    vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))

    result <- DESeq2::results(dds, alpha = fdr, name = name)

    return(list(dds = dds, vst = vst, res = result, design = design))

}


run_ruvseq <- function(dds, design, contrast, k = 4, p.val.thresh = 0.5, method = "ruvg", ...) {
  # get counts from DESeq2 -- don't normalize yet
  counts <- DESeq2::counts(dds, normalized = FALSE)

  set <- EDASeq::newSeqExpressionSet(counts)
  set <- EDASeq::betweenLaneNormalization(set, which="upper")

  if (method == "ruvg") {
    # Run first-pass DESeq scan to get empirical control genes
    suppressMessages({
      de_res <- run_differential(dds, design = design, contrast = contrast, ...)
    })

    # Use non-significant genes as as empirical control genes
    not.sig <- rownames(de_res$res)[which(de_res$res$pvalue > p.val.thresh)]
    empirical <- rownames(set)[rownames(set) %in% not.sig]
  } else if (method == "ruvs") {
    empirical <- rownames(set)
  }

  print(glue("info: using {length(empirical)} control genes to estimate variation.."))

  result <- list()
  for (num_latent in 1:k) {
    # Run specified method of unwanted variation control
    if (method == "ruvg") {
      print(glue("info: estimating k={num_latent} factors using RUVg.."))
      new_set <- RUVSeq::RUVg(set, empirical, k = num_latent) # Run RUVSeq
    } else if (method == "ruvs") {
      print(glue("info: estimating k={num_latent} factors using RUVs.."))
      groups <- RUVSeq::makeGroups(SummarizedExperiment::colData(dds)[, contrast[1]])
      new_set <- RUVSeq::RUVs(set, empirical, k = num_latent, groups) # Run RUVSeq
    }

    ddsruv <- dds
    for (w in 1:num_latent) {
      col_idx <- glue("W_{w}")
      ddsruv[[col_idx]] <- new_set[[col_idx]]
    }

    new_design <- glue("{design} + {paste('W', 1:num_latent, collapse = '+', sep = '_')}")

    suppressMessages({
      de_res <- run_differential(ddsruv, design = new_design, contrast = contrast, ...)
    })

    result[[new_design]] <- list(de = de_res, W = pData(new_set), normCounts = new_set@assayData$normalizedCounts)
  }

  return(result)
}


run_ruvseq_continuous <- function(dds, design, name, k = 4, p.val.thresh = 0.5, method = "ruvg", ...) {
  # get counts from DESeq2 -- don't normalize yet
  counts <- DESeq2::counts(dds, normalized = FALSE)

  set <- EDASeq::newSeqExpressionSet(counts)
  set <- EDASeq::betweenLaneNormalization(set, which="upper")

  if (method == "ruvg") {
    # Run first-pass DESeq scan to get empirical control genes
    suppressMessages({
      de_res <- run_differential_continuous(dds, design = design, name = name, ...)
    })

    # Use non-significant genes as as empirical control genes
    not.sig <- rownames(de_res$res)[which(de_res$res$pvalue > p.val.thresh)]
    empirical <- rownames(set)[rownames(set) %in% not.sig]
  } else if (method == "ruvs") {
    empirical <- rownames(set)
  }

  print(glue("info: using {length(empirical)} control genes to estimate variation.."))

  result <- list()
  for (num_latent in 1:k) {
    # Run specified method of unwanted variation control
    if (method == "ruvg") {
      print(glue("info: estimating k={num_latent} factors using RUVg.."))
      new_set <- RUVSeq::RUVg(set, empirical, k = num_latent) # Run RUVSeq
    } else if (method == "ruvs") {
      print(glue("info: estimating k={num_latent} factors using RUVs.."))
      groups <- RUVSeq::makeGroups(SummarizedExperiment::colData(dds))
      new_set <- RUVSeq::RUVs(set, empirical, k = num_latent, groups) # Run RUVSeq
    }

    ddsruv <- dds
    for (w in 1:num_latent) {
      col_idx <- glue("W_{w}")
      ddsruv[[col_idx]] <- new_set[[col_idx]]
    }
    init_vars <- str_split_i(design, "~", 2)
    new_design <- glue("{design} + {paste('W', 1:num_latent, collapse = '+', sep = '_')}")
      print(new_design)

    suppressMessages({
      de_res <- run_differential_continuous(ddsruv, design = new_design, name = name, ...)
    })

    result[[new_design]] <- list(de = de_res, W = pData(new_set), normCounts = new_set@assayData$normalizedCounts)
  }

  return(result)
}


# Function to get asterisks for correlation plot
labs.function <- function(x){
    case_when(x >= 0.05 ~ "",
              x < 0.05 & x >= 0.01 ~ "*",
              x < 0.01 & x >= 0.001 ~ "**",
              x < 0.001 ~ "***")
}