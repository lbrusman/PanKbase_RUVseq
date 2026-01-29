suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(parallel))
suppressMessages(library(readr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressMessages(library(patchwork))
suppressMessages(library(RUVSeq))


outdir <- "/tscc/nfs/home/lebrusman/Gaulton_lab/code/RUVseq_pankbase/outputs/250903_outs/pseudobulk_counts/"

data <- readRDS("/tscc/nfs/home/lebrusman/Gaulton_lab/data/ADA_object/250424_ADA_object_metadata_v3_3.rds")
print(data)

data@meta.data$coarse_annot <- gsub(" ", "", data@meta.data$coarse_annot)
Idents(data) <- data@meta.data$coarse_annot
unique_cell_types <- unique(data$coarse_annot)
samples <- unique(data@meta.data$samples)

for (s in c("SRR27326986", "SRR27326987", "SRR27326992", "SRR27326993",
            "SRR27326994", "SRR27326995", "SRR27326996", "SRR27326997")) {
    data@meta.data[data@meta.data$samples == s, "samples"] <- paste0(data@meta.data[data@meta.data$samples == s, "samples"],
                                                                    "__", data@meta.data[data@meta.data$samples == s, "treatments"])
}


# do pseudobulking

## Get barcodes from each sample
sample_bcs <- list()
samples <- unique(data@meta.data$samples)
for (sample in samples){
    sample_bcs[[sample]] <- row.names(data[[]][data[[]]$samples == sample,])
}

# code adapted from Elgamal et al. https://diabetesjournals.org/diabetes/article/72/11/1719/153488/An-Integrated-Map-of-Cell-Type-Specific-Gene
data_matrices <- data

get_per_sample_gex_SUMS <- function(cell.type, mtx.fp){
    print(cell.type)
    #Pull out rows of gex.counts where barcode Ident matches cell.type
    bcs <- names(Idents(data_matrices)[Idents(data_matrices) == cell.type])
    counts <- gex.counts[,colnames(gex.counts) %in% bcs]
    print(dim(counts))

    #Initialize the sample gex matrix
    counts.df <- as.data.frame(rep(0,length(row.names(gex.counts))))
    row.names(counts.df) <- row.names(gex.counts)
    colnames(counts.df) <- c('temp')

    #Loop through samples and calculate sum of gex values
    for (sample in samples){
        sample_cols <- colnames(counts) %in% sample_bcs[[sample]]
        counts.cut <- counts[,sample_cols]
        
        #If only one barcode, this becomes a vector which is an issue
        if (typeof(counts.cut) == 'double'){
            sum.counts <- counts.cut
        #If there are no barcodes, this will return NA (just return 0 for everything)
        } else if(length(colnames(counts.cut)) == 0){
            sum.counts <- rep(0,length(row.names(counts)))
        } else {
            sum.counts <- rowSums(counts.cut)
        }
        counts.df <- cbind(counts.df,as.data.frame(sum.counts))
     }
    fin.counts.df <- counts.df[,-c(1)]
    colnames(fin.counts.df) <- samples
    head(fin.counts.df)

    #Export cell type specific gene by sample matrices
    write.table(fin.counts.df, mtx.fp, sep='\t', quote=FALSE)
}


DefaultAssay(data) <- 'RNA'
gex.counts <- GetAssayData(data,slot='counts') #get gene expression counts

#Run function to make matrices
for (cell.type in unique_cell_types){
    fp <- paste0(outdir, cell.type, "_sample_gex_total_counts.txt")
    get_per_sample_gex_SUMS(cell.type, fp)
}

