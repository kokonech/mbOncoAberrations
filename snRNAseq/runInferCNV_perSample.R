setwd("/b06x-isilon/b06x-m/mbCSF/results/humanTumor/cnvProfiling/globalInferCNV")

library(infercnv)
library(config)
library(stringr)

args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    inFile = args[1]
}

print(paste("Read config:", inFile))
cfg <- config::get(file=inFile)

print(paste("Processing",cfg$resName))
print(paste("Input:",cfg$countsFile))
print(paste("Annotation:",cfg$annFile))
print(paste("Normal clusters:",cfg$refGroup))
print(paste("Cut off:",cfg$cutOff))

numClusters <- 1
groupUsage <- TRUE
if (!is.null(cfg$numClusters)) {
    numClusters <- cfg$numClusters
    groupUsage <- FALSE
}
print(paste("Num tumor clusters:",numClusters)) 
  
refGroups <- str_split(cfg$refGroups,",")[[1]]
geneOrderRef = "/omics/odcf/analysis/OE0290_projects/Ependymoma/annotations/cell10x/refdata-cellranger-hg19-1.2.0/genes/genes.inferCNV_ref.sorted.txt"

if (!is.null(cfg$customRef)) {
    print(paste("Assign custom reference:",cfg$customRef))
    geneOrderRef = cfg$customRef
}

# specific:  selected cluster as reference
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=cfg$countsFile,
                                    annotations_file=cfg$annFile,
                                    delim="\t",
                                    gene_order_file=geneOrderRef,
                                    ref_group_names=refGroups,
                                    chr_exclude=c("Y","MT") # default
                                    #chr_exclude=c("MT")
)

#options(scipen = 100) # for subclusters


infercnv_obj = infercnv::run(infercnv_obj,
                            # cutoff: 1 for SmartSeq, 0.1 for 10x Genomics
                             cutoff=cfg$cutOff,
                             out_dir=cfg$resName,
                             cluster_by_groups=groupUsage,
                             k_obs_groups =  numClusters, 
                             smooth_method="runmeans",
                             analysis_mode="subclusters",
                             output_format = "pdf", # issue with atac
                             denoise=T,
                             plot_steps=F,
                             HMM=F
                    )

print("Done!")
