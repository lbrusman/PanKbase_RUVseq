suppressMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(Seurat)
    library(data.table)
})


# Read in rds from local
data <- readRDS("/tscc/projects/ps-gaultonlab/lebrusman/from_TSCC_251212/Gaulton_lab/data/ADA_object/250424_ADA_object_metadata_v3_3.rds")

# Fix some issues in the single cell object
for (s in c("SRR27326986", "SRR27326987", "SRR27326992", "SRR27326993",
            "SRR27326994", "SRR27326995", "SRR27326996", "SRR27326997")) {
    data@meta.data[data@meta.data$samples == s, "samples"] <- paste0(data@meta.data[data@meta.data$samples == s, "samples"],
                                                                    "__", data@meta.data[data@meta.data$samples == s, "treatments"])
}

# Get metadata for chemistry info
meta <- data@meta.data
chem_meta <- meta %>% distinct(samples, chemistry)

dir.create(file.path("metadata_files"))
write.csv(chem_meta, "metadata_files/chemistry_metadata.csv", row.names = FALSE)


# Get cell count metadata
Idents(data) <- "cell_type"
counts <- table(Idents(data), data$samples)
counts <- as.data.frame(counts)

# Pivot to wide df
counts_wide <- pivot_wider(counts, names_from = Var1, values_from = Freq)
names(counts_wide)[names(counts_wide) == "Var2"] <- "samples"
counts_wide$samples <- as.character(counts_wide$samples)

write.csv(counts_wide, "metadata_files/sc_object_cellcounts.csv", row.names = FALSE)



