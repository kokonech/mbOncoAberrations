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


### WGS ###

This folder contains scripts to reproduce the WGS analysis. To start with, adjust the paths and directories in [Settings.R](WGS/Settings.R). Install the missing libraries as listed in [Settings.R](WGS/Settings.R). In particular, install the R package [NBevolution](https://github.com/VerenaK90/Neuroblastoma_evolution/tree/main/NBevolution_0.0.0.9000.tar.gz). Then source `Settings.R`. Take care that all input and output directories exist.

#### Cohort overview ####

The script [Cohort_characteristics.R](WGS/Cohort_characteristics.R) plots the distribution of subtypes that went into the analysis as shown in **Fig. 3a** and **Extende Fig. 3a**. 

#### Mutation timing ####

SNV densities were quantified with the script [Mutation_density_quantification.R](WGS/SNVdensities/Mutation_density_quantification.R). Run (Adjust_purity.R)[WGS/Adjust_purity.R) first to collect all ploidies and purities along the cohort and adjust the purity in cases where it could not be estimated by ACEseq. The script [Mutation_density_quantification.R](WGS/SNVdensities/Mutation_density_quantification.R) generates the plots shown in **Extended Data Fig. 3b-e** and **Fig. 3b and d**.

#### Population-genetics modeling ####

This was performed in 2 steps. First, we fit a population-genetics model of tumor initiation two the SNV densities at MRCA in G3/4 tumors. To re-run the analysis, first generate the summary data by sourcing the script (Input_data_MB.R)[WGS/PopGen/Input_data_MB.R]. In a second step, the model is fit to the data. Here, you need to install [pyABC](https://pyabc.readthedocs.io/en/latest/). The model fit is done by executing the script [Expansion_decay_continuous_evol.py](WGS/PopGen/Expansion_decay_continuous_evol.py), which sources the script [Expansion_decay_continuous_evol.R](WGS/PopGen/Expansion_decay_continuous_evol.R), which, in turn, sources the script [ABC.R](WGS/PopGen/ABC.R). Parameter estimation is performed according to the following pseudo-code:

**FOR EACH** optimization step: 
- Sample a parameter set from the prior probabilities given in Extended Data Table 6 
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
 to enforce good fits for the initial expansion phase. The experimental incidences were computed as 
```math
I_{\mathrm{MRCA,exp},i} =\sum_{\tilde{m}_\mathrm{MRCA,exp}} < i ; i \in \tilde{m}_\mathrm{MRCA,exp},
```
 and uncertainties were estimated to 
```math
\Delta I_{\mathrm{MRCA,exp},i} =\sum{\tilde{m}_{\mathrm{MRCA,exp},l} < i} - \sum{\tilde{m}_{\mathrm{MRCA,exp},u} < i}  ; i \in \tilde{m}_\mathrm{MRCA,exp},
```
 where $`\tilde{m}_{\mathrm{MRCA,exp},l}`$ and $`\tilde{m}_{\mathrm{MRCA,exp},u}`$ denote the lower and upper bounds of the 95% confidence interval of $`\tilde{m}_{\mathrm{MRCA,exp}}`$, respectively. Incidences and uncertainties of the ECA, $`I_{\mathrm{ECA,exp},i}`$ were computed in analogy. 

**DONE**

In the second step, we fit a model of tumor growth to the subclonal tail of 39 Group3/4 medulloblastomas with sufficient data quality. To reproduce this part, run [Neutral_fit_pre_clonal_and_clonal.py](WGS/PopGen/Neutral_fit_pre_clonal_and_clonal.py) in conjunction with [Mutation_accumulation_fit_pre_clonal_and_clonal.R](WGS/PopGen/Mutation_accumulation_fit_pre_clonal_and_clonal.R) and [Mutation_accumulation_exponential_growth.R](WGS/PopGen/Mutation_accumulation_exponential_growth.R). The fitting procedure encompasses the following steps:

1.	Compute summary statistics from measured data:

**FOR** each copy number state $1\le k \le4$:<br />
 **IF** $`g_k<10^8 \mathrm{bp}`$ , where $g_k$ is the number of bases with copy number $k$ in the genome:<br />
 **NEXT**<br />
 **ELSE:** <br />
- Compute the total number of clonal mutations as the sum of the different clonal VAF peaks constituted by amplified and non-amplified clonal mutations. Denoting with $\hat{C_k}$ the average coverage from all mutations falling on segments with copy number $k$, we classified mutations as amplified clonal if $`\frac{Q_{l-1}^{0.95}}{\hat{C_k}} < \mathrm{VAF}_k \le \frac{Q_l^{0.95}}{\hat{C_k}}`$,  where $`Q_l^{0.95}`$ is the 95% quantile of a binomial distribution with success probability \frac{\rho l}{\zeta}$, where $\rho$ is the purity, $\zeta$ is the average copy number $l$ is the B-allele count.
- Merge these mutations with those of the non-amplified clonal peak by multiplying their frequencies with $l/k$ and adding them $l$ times. 	
- Compute the cumulative mutation counts of the measured data, $`F_{k,\mathrm{exp}}(f) = \sum \mathrm{VAF}_k>f `$, where $f$ runs from 0.05 to 1 in steps of size 0.05 and  extrapolate them to the whole genome by multiplication with $\frac{\sum_k g_k}{g_k}$.<br />
**DONE**

2. Fit the model to the data: 

**FOR EACH** optimization step:<br />
-Sample values for $`\mu, \frac{\delta_\mathrm{T}}{\lambda_\mathrm{T}}, n_\mathrm{clonal}`$ and $\Delta_\rho$ from the prior distributions given in Extended Data Table 7, where $\Delta_\rho$ is a correction factor to the purity estimate by ACEseq.<br />
	**FOR EACH** copy number state $1\le k \le 4$:__
		**IF** $g_k<10^8$ bp, where $g_k$ is the length of the genome at copy number $k$:<br />
			**NEXT**<br />
		**ELSE:** <br />
			- Determine $n_(f,k)$ from equation (8), assuming a tumor size of $10^9$ cells at diagnosis, and evaluating equation (8) in bins of size 0.05 at the lower limit:
```math
n_{f,k}\approx
\begin{cases}
\sum_{i=f10^9}^{(f+0.05)10^9} S_k(i,\mu) + k\times n_\mathrm{clonal}, \quad \mathrm{if} \quad f = 0.95,\\
\sum_{i=f10^9}^{(f+0.05)10^9} S_k(i,\mu), \quad \mathrm{else},
\end{cases}
```
where $`n_\mathrm{clonal}`$  is the number of clonal variants per haploid genome already present in the tumor’s MRCA.
			- Sample for each mutation a sequencing coverage $C_k$  according to $Pois(\hat{C_k})$, where $\hat{C_k}$ is the average coverage at copy number $k$ in the data.
			- Sample for mutation a VAF according to a Binomial distribution with $C_k$ draws and success probability $\frac{f \min{\rho+\Delta_\rho,1}}{\zeta}$, where $\rho$ is the tumor cell content estimated by ACEseq, $\Delta_\rho$ is a correction factor for the purity estimate, and $\zeta$ is the average copy number at a locus with tumor copy number $k$ in the impure sample (equation (2)).
			- Compute the cumulative mutation counts, $`F_{k\mathrm{sim}}(f)= \sum \mathrm{VAF}_k>f `$, where $f$ runs from 0.05 to 1 in steps of size 0.05.
			- Compute the cost function:
```math
d=\sum_k \sum_f(F_{k,\mathrm{sim}} - F_{k,\mathrm{exp}})^2 \frac{g_k}{\sum_{k'}g_{k'}}
```
**DONE** <br />
**DONE**

The model fits were integrated with the script [Compute_evolutionary_parameters_from_growth_model.R](WGS/PopGen/Compute_evolutionary_parameters_from_growth_model.R) and visualized with the script [Plot_dynamic_model.R](WGS/PopGen/Plot_dynamic_model.R), generating **Fig. 3e,f,g** and **Extended Data Fig. 3g**. Individual model fits to the subclonal tail can be visualized with the script [Plot_neutral_fit.R](WGS/PopGen/Plot_neutral_fit.R).


#### Oncoprint ####

Significantly enriched large-scale CNVs were identified with the script [Gains_and_losses.R](\WGS/Oncoprint/Gains_and_losses.R) (which generates **Extended Data. Fig. 3h**). The regions were then integrated with the timing information in the script [CNV_timing_overview.R](WGS/Oncoprint/CNV_timing_overview.R). Finally, the script [Driver_mutations.R][WGS/Oncoprint/Driver_mutations.R] generates the oncoprint in **Fig. 3h** and the statistics shown in **Fig. 3i-l**.








