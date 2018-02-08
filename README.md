# BAYSPLINE
Bayspline is a Bayesian calibration for the alkenone paleothermometer, as published by Tierney & Tingley (2018), an accepted manuscript at Paleoceanography & Paleoclimatology (http://doi.org/10.1002/2017PA003201). This repository currently contains Matlab code to run the prediction model, i.e., to predict sea-surface temperatures from UK'37 data. The main script is Predict_Uk.m and it requires the posterior distribution of model paramaters, bayes_posterior_v2.mat, to run. The code uses Matlab's parallel processing ability to do Metropolis-Hastings estimation and takes anywhere from 20-90 seconds to run depending on the computer and the size of the input data. A forward model will be added shortly, and a Python version of BAYSPLINE is under development. 
