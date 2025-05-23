This repository contains R/Python source code for the data analysis in the context of a study in Nature focused on investigation of oncogenic aberrations in medulloblastoma Group 3/4 brain tumors:  ["Oncogene aberrations drive medulloblastoma progression, not initiation"](https://www.nature.com/articles/s41586-025-08973-5) 

The folders are constructed based on the analysis of specific data types. Check the [wiki](https://github.com/kokonech/mbOncoAberrations/wiki) for some tutorials on the certain anlaysis steps. If case of any questions, please create an issue in github repo or contact the authors directly.

## Directory contents ##

### snRNAseq ###

Analysis of single nuclei gene expression data counts matrix computed by [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) and futher filtered from quality control inspections. More details about these procedures are provided in the [corresponding study](https://www.biorxiv.org/content/10.1101/2024.02.09.579680v1.full).

General processing (normalization, clustering, UMAP, output) per sample:  [processSnRNAseq.R](snRNAseq/processSnRNAseq.R)

Merged samples processing (UMAPs, markers visualization): [mergeSnRNAseq.R](snRNAseq/mergeSnRNAseq.R)

Copy number profiling per sample: [runInferCNV_perSample.R](snRNAseq/runInferCNV_perSample.R)

Gene set varaince enrichment analysis per sample: [runGsvaPerCell.R](snRNAseq/runInferCNV_perSample.R)

R session details: [sessionInfo.RNAseq.txt](snRNAseq/sessionInfo.RNAseq.txt)

### snATACseq ###

Analysis of single nuclei ATAC chromatin data generated by [CellRanger-Arc](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc?src=social&lss=wechat&cnm=soc-we-ra_g-program&cid=7011P000000oOU6) toolkit.

General processing (filtering, normalization, preparation for InferCNV) per sample: [processAtacTumorSampleSV4.R](snATACseq/processAtacTumorSampleSV4.R)

Merged samples processing (UMAPs): [mergeAtac.R](snATACseq/mergeAtac.R)

Analysis of cis-regulatory elements (subclone-specific, gene assoications,etc): [list of scripts](snATACseq/creAnalysis)

Analysis of mutations (SComatic script, filter somatic based on bulk WGS,  plot results): [list of scripts](snATACseq/mutAnalysis)

R session details: [sessionInfo.ATAC.txt](snATACseq/sessionInfo.ATAC.txt)

### scSpatial ###

Analysis of single cell spatial RNA data obtained from the protocol Molecular Cartography of [Resolve Bioscience](https://resolvebiosciences.com/). 

General per sample processing (filtering, normalization, visualziation): [resolveSeurat.R](scSpatial/resolveSeurat.R)

Merged samples processing (UMAPs): [resolveMerged.R](scSpatial/resolveMerged.R)

Transfer of the signals obtained from snRNA-seq data into spatial: [transferAnchorsPerSample.R](scSpatial/transferAnchorsPerSample.R)

Additional processing with Giotto R package (neigbourhood, etc): [giotoAnalysis.R](scSpatial/giotoAnalysis.R)

R session details: [sessionInfo.spatial.txt](scSpatial/sessionInfo.spatial.txt) 

### ResolveOME data ###

DNA processing (prepare alignment/QC; server jobs; draw CNV results): [list of scripts](resolveBio/dna)

RNA processing (prepare alignment; server jobs for alignment and compute counts; process counts, transfer from 10X snRNAseq): [list of scripts](resolveBio/rna)


### WGS ###

This folder contains scripts to reproduce the WGS analysis. To start with, download the data from Mendeley (doi: 10.17632/g4r22w9pp8.1), adjust the paths and directories in [Settings.R](WGS/Settings.R) and provide Extended Data Table 5 in the folder "Meta_data". Install the missing libraries as listed in [Settings.R](WGS/Settings.R). In particular, install the R package [NBevolution](https://github.com/VerenaK90/Neuroblastoma_evolution/tree/main/NBevolution_0.0.0.9000.tar.gz). Then source `Settings.R`. Take care that all input and output directories exist.

#### Cohort overview ####

The script [Cohort_characteristics.R](WGS/Cohort_characteristics.R) plots the distribution of subtypes that went into the analysis as shown in **Fig. 3a** and **Extended Fig. 3a**. 

#### Mutation timing ####

SNV densities were quantified with the script [Mutation_density_quantification.R](WGS/SNVdensities/Mutation_density_quantification.R). Run [Adjust_purity.R](WGS/Adjust_purity.R) first to collect all ploidies and purities along the cohort and adjust the purity in cases where it could not be estimated by ACEseq. The script [Mutation_density_quantification.R](WGS/SNVdensities/Mutation_density_quantification.R) generates the plots shown in **Extended Data Fig. 6b-f**, **Fig. 3b and d**, and the statistics shown in **Fig. 3i,j**. For comparison, the script [MutationTimeR.R](WGS/SNVdensities/MutationTimeR.R) computes mutation densities using MutationTimeR (Gerstung et al., Nature, 2020); this script generates the output shown in **Extended Data Fig. 6g** and the input data for the shinyApp.

#### Population-genetics modeling ####

This was performed in 2 steps. First, we fit a population-genetics model of tumor initiation to the SNV densities at MRCA in G3/4 tumors. To re-run the analysis, first generate the summary data by sourcing the script [Input_data_MB.R](WGS/PopGen/Input_data_MB.R). In a second step, the model is fit to the data. Here, you need to install [pyABC](https://pyabc.readthedocs.io/en/latest/). The model fit is done by executing the script [Expansion_decay_continuous_evol.py](WGS/PopGen/Expansion_decay_continuous_evol.py), which sources the script [Expansion_decay_continuous_evol.R](WGS/PopGen/Expansion_decay_continuous_evol.R), which, in turn, sources the script [ABC.R](WGS/PopGen/ABC.R). Parameter estimation is performed according to the following pseudo-code:

**FOR EACH** optimization step: 
- Sample a parameter set from the prior probabilities given in Table 1 below 
- Sample for each tumor a time point of the MRCA, $t_2$, using the function `Sample.timepoint()` from NBevolution (the function accounts for an overall incidence of $10^{–5}$). 
- Sample for each tumor a time point of the ECA, $t_1$ using the function `Sample.timepoint.eca()` from NBevolution. 
- Sample for each tumor a neutral mutation count at ECA and MRCA using the function `Sample.mutations()` from NBevolution. Determine the mutation density by dividing with the haploid genome length of 3.3× 10^9, yielding $`\tilde{m}_\mathrm{MRCA,sim}$ and $\tilde{m}_\mathrm{ECA,sim}`$.
- Determine the simulated incidence of the ECA at the experimentally determined mutation loads, $`I_{\mathrm{ECA,sim},i} = \sum \tilde{m}_\mathrm{ECA,sim} < i; i \in \tilde{m}_\mathrm{ECA,exp}`$, and the simulated incidence of the MRCA at the experimentally determined mutation loads,  $`I_{\mathrm{MRCA,sim},i} = \sum \tilde{m}_\mathrm{MRCA,sim} < i; i \in \tilde{m}_\mathrm{MRCA,exp}`$ (where "sim" and "exp" denote simulated and experimentally determined mutation densities).
- Determine the simulated incidence of the MRCA at age 50 years,  $`I_\mathrm{MRCA,sim,ten years}`$ using the function `P.2nd.driver()` from NBevolution. 
- Compute the cost function
 $`d=\sum_{i \in m_\mathrm{MRCA,exp}} w_i  \frac{(I_{\mathrm{MRCA,sim},i} -I_{\mathrm{MRCA,exp},i})^2}{\Delta I_{\mathrm{MRCA,exp},i}^2}+\frac{(I_{\mathrm{ECA,sim},i} -I_{\mathrm{ECA,exp},i})^2}{\Delta I_{\mathrm{ECA,exp},i}^2 }+\frac{(I_\mathrm{MRCA,sim,ten years} -10^{-5})}{(10^{-4} )^2}`$,
where the weights were chosen as
```math
w_i=\begin{cases}
10, \quad \mathrm{if} \quad \tilde{m}_\mathrm{MRCA} \le 0.2/10^6 \\
1, \quad \mathrm{if} \quad \tilde{m}_\mathrm{MRCA} > 0.2/10^6,
\end{cases}
```
&nbsp;&nbsp;&nbsp;&nbsp; to enforce good fits for the initial expansion phase. The experimental incidences were computed as 
```math
I_{\mathrm{MRCA,exp},i} =\sum_{\tilde{m}_\mathrm{MRCA,exp}} < i ; i \in \tilde{m}_\mathrm{MRCA,exp},
```
&nbsp;&nbsp;&nbsp;&nbsp; and uncertainties were estimated to 
```math
\Delta I_{\mathrm{MRCA,exp},i} =\sum{\tilde{m}_{\mathrm{MRCA,exp},l} < i} - \sum{\tilde{m}_{\mathrm{MRCA,exp},u} < i}  ; i \in \tilde{m}_\mathrm{MRCA,exp},
```
&nbsp;&nbsp;&nbsp;&nbsp; where $`\tilde{m}_{\mathrm{MRCA,exp},l}`$ and $`\tilde{m}_{\mathrm{MRCA,exp},u}`$ denote the lower and upper bounds of the 95% confidence interval of $`\tilde{m}_{\mathrm{MRCA,exp}}`$, respectively. Incidences and uncertainties of the ECA, $`I_{\mathrm{ECA,exp},i}`$ were computed in analogy. 

**DONE**

**Table 1. Prior probabilities to model medulloblastoma initiation**

| Parameter  | Description | Unit | Distribution | Min | Max |
| ------------- | ------------------------------------------- | -- | ----- | -- | --- |
| $`\delta_1`$  | Relative loss rate during tissue expansion  | 1 | uniform | 0 | 0.9 |
| $`N`$  | Number of cells in tissue at peak  | 1 | log-uniform | $`10^3`$ | $`10^9`$ |
| $`\delta_2`$ | Relative loss rate during tissue contraction | 1 | uniform | 1.0001 | 1.5 |
| $`r`$ | Reduction in cell loss due to first driver | 1 | uniform | 0 | 0.999 |
| $`\nu_2`$ | Survival probability attributed to 2nd driver | 1 | uniform | 0.01 | 0.49 |
| $`\mu_{D1}`$ | Mutation rate 1st driver | 1 | log-uniform | $`10^{-9}`$ | $`10^{-4}`$ |
| $`\mu_{D2}`$ | Mutation rate 2nd driver | 1 | log-uniform | $`10^{-9}`$ | $`10^{-4}`$ |
| $`\mu`$ | Neutral mutation rate | 1/division | uniform | 1 | 15 |



In the second step, we fit a model of tumor growth to the subclonal tail of 39 Group3/4 medulloblastomas with sufficient data quality. To reproduce this part, run [Neutral_fit_pre_clonal_and_clonal.py](WGS/PopGen/Neutral_fit_pre_clonal_and_clonal.py) in conjunction with [Mutation_accumulation_fit_pre_clonal_and_clonal.R](WGS/PopGen/Mutation_accumulation_fit_pre_clonal_and_clonal.R) and [Mutation_accumulation_exponential_growth.R](WGS/PopGen/Mutation_accumulation_exponential_growth.R). The fitting procedure encompasses the following steps:

1.	Compute summary statistics from measured data:

**FOR** each copy number state $1\le k \le4$:<br />
&nbsp;&nbsp;&nbsp;&nbsp; **IF** $`g_k<10^8 \mathrm{bp}`$ , where $g_k$ is the number of bases with copy number $k$ in the genome:<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **NEXT**<br />
&nbsp;&nbsp;&nbsp;&nbsp; **ELSE:** <br />
- Compute the total number of clonal mutations as the sum of the different clonal VAF peaks constituted by amplified and non-amplified clonal mutations. Denoting with $\hat{C_k}$ the average coverage from all mutations falling on segments with copy number $k$, we classified mutations as amplified clonal if $`\frac{Q_{l-1}^{0.95}}{\hat{C_k}} < \mathrm{VAF}_k \le \frac{Q_l^{0.95}}{\hat{C_k}}`$,  where $`Q_l^{0.95}`$ is the 95% quantile of a binomial distribution with success probability \frac{\rho l}{\zeta}$, where $\rho$ is the purity, $\zeta$ is the average copy number $l$ is the B-allele count.
- Merge these mutations with those of the non-amplified clonal peak by multiplying their frequencies with $l/k$ and adding them $l$ times. 	
- Compute the cumulative mutation counts of the measured data, $`F_{k,\mathrm{exp}}(f) = \sum \mathrm{VAF}_k>f `$, where $f$ runs from 0.05 to 1 in steps of size 0.05 and  extrapolate them to the whole genome by multiplication with $\frac{\sum_k g_k}{g_k}$.<br />
**DONE**

2. Fit the model to the data: 

**FOR EACH** optimization step:<br />
- Sample values for $`\mu, \frac{\delta_\mathrm{T}}{\lambda_\mathrm{T}}, n_\mathrm{clonal}`$ and $\Delta_\rho$ from the prior distributions given in Table 2 below, where $\Delta_\rho$ is a correction factor to the purity estimate by ACEseq.<br />
&nbsp;&nbsp;&nbsp;&nbsp;**FOR EACH** copy number state $1\le k \le 4$:__
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**IF** $g_k<10^8$ bp, where $g_k$ is the length of the genome at copy number $k$:<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**NEXT**<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**ELSE:** <br />
			- Determine $n_{f,k}$ from equation (8), assuming a tumor size of $10^9$ cells at diagnosis, and evaluating equation (8) in bins of size 0.05 at the lower limit:
```math
n_{f,k}\approx
\begin{cases}
\sum_{i=f10^9}^{(f+0.05)10^9} S_k(i,\mu) + k\times n_\mathrm{clonal}, \quad \mathrm{if} \quad f = 0.95,\\
\sum_{i=f10^9}^{(f+0.05)10^9} S_k(i,\mu), \quad \mathrm{else},
\end{cases}
```
&nbsp;&nbsp;&nbsp;&nbsp;where $`n_\mathrm{clonal}`$  is the number of clonal variants per haploid genome already present in the tumor’s MRCA.
			- Sample for each mutation a sequencing coverage $C_k$  according to $Pois(\hat{C_k})$, where $\hat{C_k}$ is the average coverage at copy number $k$ in the data.
			- Sample for mutation a VAF according to a Binomial distribution with $C_k$ draws and success probability $\frac{f \min{\rho+\Delta_\rho,1}}{\zeta}$, where $\rho$ is the tumor cell content estimated by ACEseq, $\Delta_\rho$ is a correction factor for the purity estimate, and $\zeta$ is the average copy number at a locus with tumor copy number $k$ in the impure sample (equation (2)).
			- Compute the cumulative mutation counts, $`F_{k\mathrm{sim}}(f)= \sum \mathrm{VAF}_k>f `$, where $f$ runs from 0.05 to 1 in steps of size 0.05.
			- Compute the cost function:
```math
d=\sum_k \sum_f(F_{k,\mathrm{sim}} - F_{k,\mathrm{exp}})^2 \frac{g_k}{\sum_{k'}g_{k'}}
```
&nbsp;&nbsp;&nbsp;&nbsp;**DONE** <br />
**DONE**

The model fits were integrated with the script [Compute_evolutionary_parameters_from_growth_model.R](WGS/PopGen/Compute_evolutionary_parameters_from_growth_model.R) and visualized with the script [Plot_dynamic_model.R](WGS/PopGen/Plot_dynamic_model.R), generating **Fig. 3e,f,g** and **Extended Data Fig. 6l**. Individual model fits to the subclonal tail can be visualized with the script [Plot_neutral_fit.R](WGS/PopGen/Plot_neutral_fit.R).

**Table 2. Prior probabilities to model medulloblastoma growth**

| Parameter  | Description | Unit | Distribution | Min | Max |
| ------------- | ------------------------------------------- | -- | ----- | -- | --- |
| $`\mu`$  | Number of SSNVs per division  | 1 | uniform | 0.1 | 20 |
| $`\frac{\delta_\mathrm{T}}{\lambda_\mathrm{T}}`$  | Relative loss during tumor growth  | 1 | uniform | $0$ | $0.99$ |
| $`\Delta_\rho`$ | Adjustment of the estimated purity ($`\pm`$) | 1 | uniform | -0.05 | 0.05 |
| $`n_{clonal}`$ | Number of clonal variants | 1 | uniform | 0 | 10000 |

#### Oncoprint ####

Significantly enriched large-scale CNVs were identified with the script [Gains_and_losses.R](\WGS/Oncoprint/Gains_and_losses.R) (which generates **Extended Data. Fig. 7a**). The regions were then integrated with the timing information in the script [CNV_timing_overview.R](WGS/Oncoprint/CNV_timing_overview.R, which also generates **Extended Data Fig. 7b**). Finally, the script [Driver_mutations.R](WGS/Oncoprint/Driver_mutations.R) generates the oncoprint in **Fig. 3h**.








