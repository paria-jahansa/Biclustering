Biclustering: Nonparametric Bayesian Biclustering in R
Overview
This repository contains R code for data simulation and the implementation of a nonparametric Bayesian biclustering method, as described in:

Lee, J., Müller, P., Zhu, Y., & Ji, Y. (2013).
A Nonparametric Bayesian Model for Local Clustering with Application to Proteomics.
Journal of the American Statistical Association, 503, 775-788.

This method is particularly useful for applications in proteomics and high-dimensional clustering problems.

Features:
1. Simulates biclustering datasets
2. Implements a Bayesian nonparametric approach for local clustering
3. Based on the Dirichlet Process Mixture Model (DPMM)

Usage
Run the main R script to perform biclustering:
source("main_script.R")
Make sure to adjust parameters in the script as needed.

Dependencies
The code requires the following R packages:

MCMCpack
coda
ggplot2
Install them using:
install.packages(c("MCMCpack", "coda", "ggplot2"))
Citation
If you use this code in your research, please cite:

 Lee, J., Müller, P., Zhu, Y., & Ji, Y. (2013).
A Nonparametric Bayesian Model for Local Clustering with Application to Proteomics.
Journal of the American Statistical Association, 503, 775-788.

Additionally, please cite this GitHub repository:

Jahansa, P. (2025). Biclustering: Nonparametric Bayesian Biclustering in R.
GitHub repository: https://github.com/paria-jahansa/Biclustering

Contributions
Contributions are welcome! Feel free to open issues or submit pull requests.
