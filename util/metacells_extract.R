# require R/4.0.max


# Main initial code from Enrique Blanco Carmona  
# Check for more details here: https://github.com/gcapture/glioma-het


library(Seurat)
library(dplyr)
options(stringsAsFactors=FALSE)

args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    sId = args[1]
    inFile = args[2]
}

resDir = dirname(inFile)

print(paste("Processing",inFile))

# METADATA VARIABLES
sample <- readRDS(inFile)
metacell_content <- 5 # Number of cells to be merged into metacells.

TARG="cnvBlock"

if (length(args) > 2) {
    print("Loading annotation...")
    annDf <- read.table(args[3])
    summary(rownames(annDf) == rownames(sample@meta.data))
    sample$seurat_clusters <- as.factor(annDf[, TARG])      

}


# METACELL GENERATION.


# Will store all the metacells. The test column will be removed at the end.
whole_metacells <- data.frame(test = rownames(sample), row.names = rownames(sample))


# Will store the complete annotation for the metacells.
whole_annotation <- data.frame(cluster_names = "test", row.names = "test")
meta_counter <- 0 # To keep a count of the metacells that are created.

# Generate a new metadata column storing the mapping cell-metacell.
sample[["metacell_mapping"]] <- "not_mapped"

# Here my sample has already been labelled, so the cluster labels are actually in the 'cluster_names' column. Use 'seurat_clusters' instead if not labelled.

# Checking list:
#Seurat::Idents(sample) <- sample$cluster_names

# for restored object
for (cluster_id in levels(sample$seurat_clusters)){
# custom version
#for (cluster_id in levels(sample$)){
    print(sprintf("Computing metacells for cluster %s.", cluster_id))
    # Will store the metacells per cluster.
    metacells <- data.frame(test = rownames(sample), row.names = rownames(sample))
    #chunksample <- sample[, sample$cluster_names == cluster_id]
    # Subset the sample by each cluster ID.
    chunksample <- sample[, sample$seurat_clusters  == cluster_id]
     # Subset the sample by each cluster ID.
    #chunksample <- sample[, sample$seurat_clusters  == cluster_id]
    # Get the count data as a data frame and transpose it so columns are GENES and rows are CELLS.
    countdata <- t(as.data.frame(Seurat::GetAssayData(chunksample, slot = "counts")))
    # Get the possible amount of metacells.
    times <- trunc(dim(countdata)[1] / metacell_content)
    for (i in seq(1,times)){
        meta_counter <- meta_counter + 1
        # Generate slice points for each metacell. i.e: 1-5, 6-10, 11-15...
        start <- ((i -1) * metacell_content + 1)
        end <- i * metacell_content
        # Compute the slice as a data frame containing the sum of the subsetted cells. dims = 1 row (metacell), X columns (genes)
        slice <- as.data.frame(colSums(countdata[start:end, ]))
        # Get the name of the cells merged.
        cell_names <- rownames(countdata[start:end, ])
        # Add the metacell.
        col_name <- sprintf("metacell_%s", meta_counter)
        metacells[[col_name]] <- slice[,1]
        # Add the metacell mapping to the cells in the sample.
        sample$metacell_mapping[colnames(sample) %in% cell_names] <- col_name
    }
    # Delete the test column as we already have more than 1 column in our data frame.
    metacells[["test"]] <- NULL
    # Will contain the annotation of the generated metacells. Columns: cluster identities. Rows: each metacell.
    annotation <- data.frame(cluster_names = colnames(metacells), row.names = colnames(metacells))
    # Replace the dummy cluster_names column's values for the actual label for the cluster.
    annotation$cluster_names <- cluster_id
    # Add the annotation data and the metacell data to the "whole" dataframe. \
    # In the end: Number of Columns for metacell object = Number of rows for annotation object.
    whole_metacells <- cbind(whole_metacells, metacells)
    whole_annotation <- rbind(whole_annotation, annotation)
}

# Turn the names into characters for the sake of avoiding errors when subsetting.
if (length(args) == 2) {
   whole_annotation$cluster_names <- paste0("cl",as.character(whole_annotation$cluster_names))
}

# Delete the test row from the global annotation data.
whole_annotation <- whole_annotation[!rownames(whole_annotation) %in% c("test"), , drop = FALSE]

# Delete the test column from the global metacell data.
whole_metacells$test <- NULL

# Path and name of the annotation file that will be used in the inferCNV call.
#cnv_analysis_folder <- "" # Path to the folder that will store the annotation file.
#dir.create(cnv_analysis_folder, recursive = TRUE)

annotation_file <- sprintf("%s/%s_annotation_metacells.tsv", resDir, sId)

# Save the annotation object.
utils::write.table(whole_annotation,
                   file = annotation_file,
                   sep = "\t",
                   row.names = TRUE,
                   col.names = FALSE,
                   quote = FALSE)

# Return the metacell object as a matrix (required for running inferCNV).
whole_metacells <- as.matrix(whole_metacells)

# It would be wise to save the Seurat object with the metacell mapping.

ann_file2 <- sprintf("%s/%s_metacells_mapping.tsv", resDir , sId)
write.table(sample@meta.data, ann_file2, sep="\t", quote=F)

mtxFile <- paste0(resDir,"/", sId, "_raw_matrix_metacells.txt.gz")
gz1 <- gzfile(mtxFile, "w")
write.table(whole_metacells, gz1, quote=F,sep="\t")
close(gz1)

