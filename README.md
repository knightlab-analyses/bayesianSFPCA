# bayesianSFPCA
#### R code for the manuscript of Bayesian SFPCA 

### required R packages: 
parallel(v.3.4.3), rstan(v.2.18.2), loo(v.2.0.0.9000), Matrix(v.1.2.14), ggplot2(v.3.1.0), splines(v.3.4.3), bayesplot(v.1.6.0)

***
#### Reproduce results for simulation I 
simulate single and double region pertrubed data (datsim_splinectomeR.R)

results on simulated single region perturbed data (dfm_singleRegion)
* observed data (sfpca_dfm.R)

  - Observed_spaghetti_group.pdf as Figure 1a; 

  - Observed_patient3.pdf as Figure 2a; 

  - Observed_patient5.pdf as Figure 2b
* results from using splinectomeR (figure_spline.R)

  - splinectome_dfm_curves.pdf as Figure 1b
* results from applying Bayesian SFPCA (sfpca_dfm.R)

  - FPCs_mean_PC1.pdf as Figure 2c;
 
  - FPCs_mean_PC2.pdf as Figure 2d; 

  - boxplot_PC1_scores.pdf as Figure 2e;
 
  - boxplot_PC2_scores.pdf as Figure 2f 

results on simulated double region perturbed data (df_doubleRegion)
* observed data (sfpca_df.R)

  - Observed_spaghetti_group.pdf as Figure 1c; 

  - Observed_patient4.pdf as Figure 3a;

  - Observed_patient8.pdf as Figure 3b
* results from using splinectomeR (figure_spline.R) 

splinectome_df_curves.pdf as Figure 1d)
* results from applying Bayesian SFPCA (sfpca_df.R)

  - FPCs_mean_PC1.pdf as Figure 3c; 

  - FPCs_mean_PC2.pdf as Figure 3d; 

  - boxplot_PC1_scores.pdf as Figure 3e;

  - boxplot_PC2_scores.pdf as Figure 3f

***
#### Reproduce results for simulation II
* simulated data with 100 total samples and 0 % missing values (N100_M0)

  - R code for simulation: simulate_N100_M0.R;

  - R code for applying Bayesian SFPCA: sfpca_N100_M0.R;

  - Figure 4a (Mean_trueVsestimated.pdf); 

  - Figure 4b (FPCs_trueVsestimated.pdf);

  - Figure 4c (fpcScores_scatterplots.pdf)

* simulated data with 100 total samples and 80 % missing values (N100_M80)

  - R code for simulation: simulate_N100_M80.R;

  - R code for using splinectomeR: figure_splinectome_N100_M80.R; 

  - R code for applying Bayesian SFPCA: sfpca_N100_M80.R;

  - Figure 4d (Mean_trueVsestimated.pdf); 

  - Figure 4e (FPCs_trueVsestimated.pdf);

  - Figure 4f (fpcScores_scatterplots.pdf) 

* simulated data with 50 total samples and 80 % missing values (N50_M80)

  - R code for simulation: simulate_N50_M80.R;

  - R code for applying Bayesian SFPCA: sfpca_N50_M80.R;

  - Figure 5a (Mean_trueVsestimated.pdf);

  - Figure 5b (FPCs_trueVsestimated.pdf);

  - Figure 5c (fpcScores_scatterplots.pdf)

* simulated data with 25 total samples and 80 % missing values (N25_M80)

  - R code for simulation: simulate_N25_M80.R;

  - R code for using splinectomeR: figure_splinectome_N25_M80.R;

  - R code for applying Bayesian SFPCA: sfpca_N25_M80.R;

  - Figure 5d (Mean_trueVsestimated.pdf);

  - Figure 5e (FPCs_trueVsestimated.pdf);

  - Figure 5f (fpcScores_scatterplots.pdf)

* simulated data with 10 total samples and 50 % missing values (N10_M50)

  - R code for simulation: simulate_N10_M50.R;

  - R code for applying Bayesian SFPCA: sfpca_N10_M50.R;

  - Figure 6a (Mean_trueVsestimated.pdf);

  - Figure 6b (FPCs_trueVsestimated.pdf);

  - Figure 6c (fpcScores_scatterplots.pdf)

* simulated data with 10 total samples and 80 % missing values (N10_M80)

  - R code for simulation: simulate_N10_M80.R;

  - R code for using splinectomeR: figure_splinectome_N10_M80.R;

  - R code for applying Bayesian SFPCA: sfpca_N10_M80.R;

  - Figure 6d (Mean_trueVsestimated.pdf);

  - Figure 6e (FPCs_trueVsestimated.pdf);

  - Figure 6f (fpcScores_scatterplots.pdf)

***
#### Reproduce results for ECAM data application
* Observed data (Shannon_modelSelection.R)

  - Figure 7a: Observed_spaghetti.pdf

  - Figure 7b: Observed_shannon_abx.pdf

  - Figure 7c: Observed_shannon_delivery.pdf

  - Figure 7d: Observed_shannon_diet.pdf

* Results from splinectomeR (figure_splinectome.R)

  - Figure 7e: splinectome_delivery_curves.pdf

  - Figure 7f: splinectome_diet_curves.pdf

* Results from Bayesian SFPCA (model selection: Shannon_modelSelection.R + model fitting: sn_3pc_1n.R)

  - Figure 8a: mean_spaghetti_transform.pdf

  - Figure 8b: FPCs_mean_PC1.pdf

  - Figure 8c: FPCs_mean_PC2.pdf

  - Figure 8d: FPCs_mean_PC3.pdf

  - Figure 9a: boxplot_abx_PC3_scores.pdf

  - Figure 9b: boxplot_delivery_PC1_scores.pdf

  - Figure 9c: boxplot_delivery_PC2_scores.pdf

  - Figure 9d: boxplot_diet_PC1_scores.pdf

  - Figure 9e: boxplot_diet_PC2_scores.pdf

  - Figure 9f: boxplot_diet_PC3_scores.pdf



