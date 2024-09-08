# OliparabDataAnalysis

This repository contains R scripts and data related to the study of olaparib treatment on prostate cancer cell lines. The study focuses on the effects of immunotherapy utilizing NK cells and ADCC mediating agents in conjunction with Poly (ADP-ribose) polymerase (PARP) inhibition.

## Study Overview
***Title: Immunotherapy Utilizing NK cells and ADCC Mediating Agents with Poly (ADP-ribose) polymerase (PARP) Inhibition***
## Overall Design:
Two prostate cancer cell lines, 22Rv1 and DU145, were tested. Each cell line was treated with 20µM olaparib or a vehicle control. Samples were harvested at 6, 12, and 18 hours. Gene expression changes were compared between olaparib-treated cells and corresponding vehicle control-treated cells at each time point.

## Sample Details:
The study includes a total of 12 samples, categorized by treatment type and duration:
22Rv1 Cell Line:
No Treatment: 3 samples (6h, 12h, 18h)
Olaparib Treatment: 3 samples (6h, 12h, 18h)
DU145 Cell Line:
No Treatment: 3 samples (6h, 12h, 18h)
Olaparib Treatment: 3 samples (6h, 12h, 18h)

## Data Processing
Normalization and Transformation:
Median Ratio Normalization:
The raw counts were normalized using the DESeq2 package, which applies median ratio normalization to account for differences in sequencing depth.
Variance Stabilizing Transformation (VST):
The variance stabilizing transformation from the DESeq2 package was applied to the normalized counts. VST aims to stabilize the variance across mean expression levels, accounting for the dependence of variance on the mean.

## Statistical Analysis:
Pearson Correlation Coefficient:
The Pearson correlation coefficient was used to compute the correlation coefficients between normal and breast tumor samples. This method measures the linear relationship between two variables, producing a value between -1 and 1

## Calculating the Correlation Matrix:
The cor() function in R was used to compute the correlation matrix of the remaining columns, calculating the Pearson correlation coefficient by default.
Squaring the Correlation Matrix:
The resulting correlation matrix was squared to obtain the coefficient of determination (R²), indicating the proportion of variance in one variable that is predictable from the other variable.

## Reference
For more details on the dataset, please refer to the GEO database: GSE121682.
