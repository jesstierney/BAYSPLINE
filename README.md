# BAYSPLINE
Bayspline is a Bayesian calibration for the alkenone paleothermometer, as published by Tierney & Tingley (2018), an accepted manuscript at Paleoceanography & Paleoclimatology (http://doi.org/10.1002/2017PA003201). This repository currently contains Matlab code to run the prediction model, i.e., to predict sea-surface temperatures from UK'37 data. The main script is UK_Predict.m and it requires the posterior distribution of model paramaters, bayes_posterior_v2.mat, to run. The code uses Matlab's parallel processing ability to do Metropolis-Hastings estimation and takes anywhere from 20-90 seconds to run depending on the computer and the size of the input data. 

NOTE: If your data are from within the following regions, the SSTs predicted are seasonal rather than annual averages (see TT2018 for a full discussion and a map of the highlighted areas):

NORTH ATLANTIC (above 48N): August-September-October
NORTH PACIFIC (above 45N): June-July-August
MEDITERRANEAN & BLACK SEA: November-May

There is also a forward model function, UK_forward.m. It models UK37 from SSTs, using the default posterior file.