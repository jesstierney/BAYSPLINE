# BAYSPLINE

BAYSPLINE is a Bayesian calibration for the alkenone paleothermometer. It uses basis splines to fit the subtle, but significant, loss of sensitivity of the U<sup>K'</sup><sub>37</sub> index to temperature at high SSTs. When using this model, please cite:

Tierney, J.E. & Tingley, M.P. (2018) BAYSPLINE: A New Calibration for the Alkenone Paleothermometer. Paleoceanography and Paleoclimatology 33, 281-301, [https://doi.org/10.1002/2017PA003201]. 

## Quick guide

To predict sea-surface temperatures from U<sup>K'</sup><sub>37</sub> data, use the script `UK_Predict.m`. The function requires only two arguments, your data ("uk", a vector of U<sup>K'</sup><sub>37</sub> values) and a prior value for the standard deviation of the time series ("pstd", in degrees Celsius). The recommended value for pstd is 5. There is an optional third argument to use a different posterior parameter file, this is for use with some the data assimilation packages we have in development, so you can ignore it. This function requires only the posterior distribution of model paramaters, `bayes_posterior_v2.mat`, to run. The code uses Matlab's parallel processing ability to do Metropolis-Hastings estimation and takes anywhere from 20-90 seconds to run depending on the number of cores on your computer and the size of the input data.

NOTE: If your data are from within the following regions, the SSTs predicted are seasonal rather than annual averages (see TT2018 for a full discussion and a map of the highlighted areas):

NORTH ATLANTIC (above 48˚N): August-September-October

NORTH PACIFIC (above 45˚N): June-July-August

MEDITERRANEAN & BLACK SEA: November-May

To predict U<sup>K'</sup><sub>37</sub> values from SSTs, use `UK_forward.m`. This function assumes that you are giving it the "correct" SSTs, i.e., if you are in the North Atlantic region defined above, it is up to the user to give it August-October SSTs. You can find out whether your location has seasonal production by the using the checkSeasonality.m script, which will output the months of production for a given list of latitudes and longitudes. The forward model function limits possible values to between 0 and 1, which are the bounds of the proxy.

`DemoScript.m` demonstrates how to use BAYSPLINE step-by-step. It loads U<sup>K'</sup><sub>37</sub> data from a .csv file (`DemoUKData.csv`), runs BAYSPLINE, and then includes some tips on how to export the results as new .csv files.