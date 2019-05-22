# bayesianSFPCA
#### R code for the manuscript of Bayesian SFPCA 

***
### required R packages: 
parallel(v.3.4.3), rstan(v.2.18.2), loo(v.2.0.0.9000), Matrix(v.1.2.14), ggplot2(v.3.1.0), splines(v.3.4.3), bayesplot(v.1.6.0)

#### Reproduce results for simulation I 
simulate single and double region pertrubed data (datsim_splinectomeR.R)
* results on simulated single region perturbed data (dfm_singleRegion)
** observed data (sfpca_dfm.R; Observed_spaghetti_group.pdf as Figure 1a; Observed_patient3.pdf as Figure 2a;
Observed_patient5.pdf as Figure 2b
** results from using splinectomeR (figure_spline.R; splinectome_dfm_curves.pdf as Figure 1b)
** results from applying Bayesian SFPCA (sfpca_dfm.R)
FPCs_mean_PC1.pdf as Figure 2c; FPCs_mean_PC2.pdf as Figure 2d; boxplot_PC1_scores.pdf as Figure 2e; 
boxplot_PC2_scores.pdf as Figure 2f 

* results on simulated double region perturbed data (df_doubleRegion)

