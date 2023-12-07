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

SNV densities were quantified with the script [Mutation_density_quantification.R](WGS/SNVdensities/Mutation_density_quantification.R). Run (Adjust_purity.R)[WGS/Adjust_purity.R) first to collect all ploidies and purities along the cohort and adjust the purity in cases where it could not be estimated by ACEseq. The script [Mutation_density_quantification.R](WGS/SNVdensities/Mutation_density_quantification.R) generates the plots shown in **Extended Data Fig. 3b-e** and **Fig. 3b and d***.

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
and uncertainties were estimated to $`\Delta I_{\mathrm{MRCA,exp},i} =\sum_{\tilde{m} _{\mathrm{MRCA,exp},l}<i -\sum_{\tilde{m}_{\mathrm{MRCA, exp},u}}<i;i \in \tilde{m}_\mathrm{MRCA,exp}`$, where `$\tilde{m}_{\mathrm{MRCA,exp},l}`$ and `$\tilde{m}_{\mathrm{MRCA,exp},u}`$ denote the lower and upper bounds of the 95% confidence interval of `$\tilde{m}_{\mathrm{MRCA,exp}}`$, respectively. Incidences and uncertainties of the ECA, $`I_{\mathrm{ECA,exp},i}`$ were computed in analogy. 
**DONE**







