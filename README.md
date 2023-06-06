# Analysis
Sample survey dataset and code for conducting data analysis in both Stata and R

This repository provides example Stata and R code to illustrate the analysis of a dataset comprising simulated results for an anti-TB drug resistance survey. 
The R code uses the simulated data contained in “drsSampledata.csv”. The Stata code uses the same data, reformatted for use in Stata software 
(“Stata_sampleDataFinal.csv”). The simulated data assumes a stratified cluster sampling survey design, involving three strata and variable cluster size. 

The R and Stata codes are similar but not identical. These codes emphasise different aspects of the analysis building upon the strengths of each software. 
For example, the Stata code details step by step all the script lines needed to produce the survey consort chart given in  
“Powerpoint_Patient_Enrollment_Flowchart_vs1.ppt”. In contrast, the R code focuses predominantly on providing example scripts for common table visualisations 
applicable to these survey types (see “Results R code VS1.ppt”).

To begin the analysis of the simulated dataset in R, download “R_Code_VS1.R” and follow the instructions. 
The results of this analysis are summarised in “Results R code VS1.ppt”. A stand-alone script to produce a plot visualising 
the size of clusters - in case of surveys involving fixed cluster size - is also provided in “ClusterSizeVisualisation_in_R_VS1.R”.

To begin the analysis of the simulated dataset in Stata, download “Stata_dofile.do” and follow the instructions. 
The results of this analysis are summarised in “Powerpoint_Patient_Enrollment_Flowchart_vs1.ppt”, “Stata_analysis_results_v1.xlsx” 
and “Stata_age_sex_distribution_v1.png”. 

