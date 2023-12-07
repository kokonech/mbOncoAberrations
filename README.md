This repository contains R/Python source code for the data analysis in the context of a study focused on investigation of oncogenic aberrations in medulloblastoma Group 3/4 brain tumors.
The folders are constructed based on the analysis of specific data types.

## Directory contents ##

### snRNAseq ###

Analysis of single nuclei gene expression data.

General processing (normalization, clustering, UMAP, output):  *processSnRNAseq.R*

Merged samples processing (UMAPs, markers visualization): _mergeSnRNAseq.R_

Copy number profiling: _runInferCNV_perSample.R_

Gene set varaince enrichment analysis: _runGsvaPerCell.R_

R session details: _sessionInfo.RNAseq.txt_

### snATACseq ###

Analysis of single nuclei ATAC chromatin data.

General processing (filtering, normalization, preparation for InferCNV): _processAtacTumorSampleSV4.R_

Merged samples processing (UMAPs): _mergeAtac.R_

R session details: _sessionInfo.ATAC.txt_






