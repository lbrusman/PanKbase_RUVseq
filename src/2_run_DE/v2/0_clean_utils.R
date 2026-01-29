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

new_filter <- function(row, min_reads = minreads, min_prop = minprop) {
    sum(row >= min_reads) >= min_prop
}

fancy_filter <- function(raw_mat, keep_meta, sample_col, contrast_var, control_grp, experimental_grp, ...) {
    if (is.factor(keep_meta[,contrast_var]) || is.character(keep_meta[,contrast_var])) {
        ctrl_samps <- keep_meta[keep_meta[[contrast_var]]==control_grp, sample_col]
        ctrl_counts <- raw_mat[,ctrl_samps]

        exp_samps <- keep_meta[keep_meta[[contrast_var]]==experimental_grp, sample_col]
        exp_counts <- raw_mat[,exp_samps]

        #get min donors for expression
        ctrl_cutoff <- floor(length(ctrl_samps)/2)
        exp_cutoff <- floor(length(exp_samps)/2)

        all_genes <- rownames(raw_mat)

        ctrl_genes_to_keep <- c()
        # loop through each row in the data frame
          for (i in 1:nrow(ctrl_counts)) {
          
            # check if there are at least n_ND values greater than 5 in the current row
            if (sum(ctrl_counts[i, ] >= 5) >= ctrl_cutoff) {
                
                # add the row name to the result list
                ctrl_genes_to_keep <- c(ctrl_genes_to_keep, all_genes[i])
            }
          }
          
          exp_genes_to_keep <- c()
          # loop through each row in the data frame
          for (i in 1:nrow(exp_counts)) {
          
            # check if there are at least n_ND values greater than 5 in the current row
            if (sum(exp_counts[i, ] >= 5) >= exp_cutoff) {
                
                # add the row name to the result list
                exp_genes_to_keep <- c(exp_genes_to_keep, all_genes[i])
            }
          }
          genes_to_keep <- unique(c(ctrl_genes_to_keep, exp_genes_to_keep))
    }

    else {
        all_samps_cutoff <- floor(nrow(keep_meta)/2)

        all_genes <- rownames(raw_mat)

        genes_to_keep <- c()
        for (i in 1:nrow(raw_mat)) {
          if (sum(raw_mat[i, ] >= 5) >= all_samps_cutoff) {
            genes_to_keep <- c(genes_to_keep, all_genes[i])
          }
        }
      }

    genes_to_keep <- unique(genes_to_keep)

    return(genes_to_keep)

}



run_ruvseq <- function(raw_mat, result, k = 10, p.val.thresh = 0.5, ...) {
    not.sig <- rownames(result)[result$pvalue > p.val.thresh]

    set <- EDASeq::newSeqExpressionSet(as.matrix(raw_mat))
    set <- EDASeq::betweenLaneNormalization(set, which="upper")

    empirical <- rownames(set)[rownames(set) %in% not.sig]
    print("length empirical:")
    print(length(empirical))

      all_result <- list()
      
      for (num_latent in 1:k) {
          new_set <- RUVSeq::RUVg(set, empirical, k = num_latent) # Run RUVSeq
          all_result[[num_latent]] <- list(W = pData(new_set), normCounts = new_set@assayData$normalizedCounts) #de = de_res, 
        }
    return(all_result)
    
}



# Function to get asterisks for correlation plot
labs.function <- function(x){
    case_when(x >= 0.05 ~ "",
              x < 0.05 & x >= 0.01 ~ "*",
              x < 0.01 & x >= 0.001 ~ "**",
              x < 0.001 ~ "***")
}