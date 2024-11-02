Elucidating Pathways Underlying Sporadic and Lynch Syndrome Forms of Colorectal Cancer
This repository contains code to analyze clinical and genetic data related to colorectal cancer (CRC), aiming to identify significant genes, conduct statistical tests, and triage data from external sources like cBioPortal.
Contents
Clinical analysis.R: Script to analyze clinical data, which includes patient information such as their age, CRC localization, and TNM stage.
Finding_genes.R: Script to identify and analyze genes associated with CRC, potentially identifying significant genetic mutations.
Fishers_test_julia.R: Performs Fisher's exact test to find statistical significant associations
cbioportal_traige_CRC.R: Script to handle CRC data from cBioPortal, filtering or analyzing it for relevant features.
Prerequisites
R: The primary language for analysis scripts.
Required Libraries: Specific libraries used within each R script (e.g., dplyr, ggplot2) should be installed.

Clinical Analysis:
Run Clinical analysis.R to analyze clinical aspects of the CRC dataset. Results will include tables and visualizations related to patient characteristics.
Finding Genes:
Use Finding_genes.R to identify genes with significant associations in CRC. This script will output lists of genes, with counts, frequencies, and p-values.
Statistical Testing with Fisher's Test:
Run Fishers_test_julia.R to perform Fisher's exact test on specified datasets. 
cBioPortal Data Triage:
Execute cbioportal_traige_CRC.R to filter and analyze CRC data from cBioPortal. The output may include a cleaned dataset or key features for further analysis.
