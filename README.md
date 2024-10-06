# ESVT
This is the code package for the manuscript “A Tale of Two Environments: Divisive Normalization and the (In)Flexibility of Choice” by Kurtz-David et al. (2024).
The manuscript presents the results from a behavioral study that consisted of two stages, spanning over three days. 
In this repository, we share the raw data and all the code required for the analyses detailed in the manuscript. 

## System requirements
The code for the code package to create the datasets for STAGE II of the study was run and tested on MATLAB 2023b on Windows 11. 
The code for all the figures was run and tested on MATLAB 2023b on Windows 11.
The code for the regression tables (Table 1 and Supp. Table 3) was run and tested on STATA 18SE with Windows 11.
The code for the pooled model estimates (Supp. Table 5) was run and tested on STATA 18SE with Windows 11.
For the subject-level estimation, please see a separate READ.MD file saved on export/ SubjectLevelEstimations.
 
## Installation guide
Install MATLAB 2023b edition.
Install STATE 18SE.

## Data
All the raw data can be found in export/RawData.
This is the output from oTree.

## Demo
All the code in this capsule was tested. 

## Estimate rho based on STAGE I choices
To estimate rho for sample data – 
Run export/Estimate_STAGEI_rho/estimate_rho.do (sample data from session 16).
Runtime should take less than a minute.

## Subject-Specific Datasets for STAGE II
To generate datasets for STAGE II of the study based on the results from STAGE I (subject-specific Pareto and Uniform distributions) - 
Run export/DatasetsGeneratorCodePackage/main.m. (sample data).
Expected output: csv files with the datasets, one file per subject.
Runtime should take up to five minutes per subject.
 
The pre-processed data from STAGE II can be found in export/figures/ESVT_data_allSubjects.xlsx and in export/regressions/ESVT_data_allSubjects.dta.

## Figures and Tables
See detailed instructions below.

### Runtime: 
Figure 1: 1h20m. 
Subject-level estimation: ~5 minutes per subject (see separate guide in export/ SubjectLevelEstimations). 
All other analyses reported below take less that 2 minutes per function.

### Figure 1
Run export/figures/fig1.m.

### Figure 2
Fig 2A-B: No code was created for these figures.
Fig 2C: run export/figures/fig2C.m.
Fig 2E: run export/figures/fig2E.m.

### Figures 3-4
1) run the subject-level estimation (see a separate guide in export/ SubjectLevelEstimations). 
2) Run export/figures/figs3_4.m.

### Figure 5
No code was created for this figure.

### Table 1
Run regressions/table_1_supp_table_3.do.

### Table 2
Run the subject-level estimation (see a separate guide in export/ SubjectLevelEstimations). 

### Supplementary Fig. 1
Run export/figures/suppFig1.m

### Supplementary Fig. 2
Run export/figures/suppFig2.m

### Supplementary Fig. 3
Run export/figures/figs3_4.m.

### Supplementary Table 1
No code was created for this table.

### Supplementary Table 2
Run export/Estimate_STAGEI_rho/estimate_rho.do (sample data from session 16).

### Supplementary Table 3
Run regressions/table_1_supp_table_3.do.

### Supplementary Table 4
Run the subject-level estimation (see a separate guide in export/ SubjectLevelEstimations). 

### Supplementary Table 5
Run regressions/supp_table_5.do.
