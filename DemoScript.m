% DEMO script: How to use BAYSPLINE to calibrate your data
%
% Step 1: You need to load your UK'37 data into Matlab. If you want you can
% just copy and paste it in as a vector, but a more streamlined approach
% would be to simply load it up from a .txt, .csv, .xlsx file. Matlab's
% function "readtable" can read data from all of these formats. Here, I'll
% use that function to load up my UK'37 from the Gulf of Aden, which is
% stored in the file, 'DemoUKData.csv'.
myData = readtable('DemoUKData.csv');
% Now your data is loaded into the workspace as a table. 
%
% Step 2: Run BAYSPLINE to get SST. I'm going to use a prior standard
% deviation of 5˚C here, because these data have high UK'37 values.
% Practically speaking, 5˚C works pretty well for most applications. If you
% have lower UK'37 values though 10˚C is a reasonable prior as well. These
% priors are all pretty "uninformative", i.e., they don't influence the
% posterior very much. If you think about it, a prior of 5˚C gives you a
% 95% confidence interval that spans -10 to +10˚C - that's pretty wide.
output = UK_predict(myData.uk37,5);
% The code above runs BAYSPLINE, and the results will be in the 'output'
% structure. Note that you don't need the optional third argument
% ('bayes').
% When you run this code, Matlab will fire up the parallel computing
% function to take advantage of however many cores you have on your
% computer. The more data you have in your file, the longer it will take.
% When it is done, two plots will appear: prior vs. posterior, and the
% calibrated time series with error bars.
% 
% The contents of the output are described in UK_predict.m. output.SST contains
% the 2.5%, 50%, and 97.5% confidence levels which is also what is plotted
% automatically by BAYSPLINE. output.ens contains the full 1000 member
% ensemble. You might want to this to do statistical analyses.
%
% To make a simple quick plot of your data with the error bars, you can do:
figure(1); clf;
plot(myData.ageBP,output.SST);
% The middle line is the median.
% 
% Step 3: Saving your results
% You can save your results in a .mat file by doing:
save('calibratedData.mat','myData','output');
% Or maybe you prefer to export the data to a .csv file. We can append the
% SST results to the myData table and then write it to file.
myData.SST_2_5 = output.SST(:,1);
myData.SST_50 = output.SST(:,2);
myData.SST_97_5 = output.SST(:,3);

writetable(myData,'DemoCalibratedData.csv');
% Now you have a new .csv file that you can use to plot the results in your
% program of choice. You can export the full ensemble as well,
% but since it has 1000 columns, it might be easier to save it in its own
% file. To do this you first have to make the ensemble a table, then write
% it.
FullEnsemble = array2table(output.ens);
writetable(FullEnsemble,'DemoFullEnsemble.csv');
% The resulting file has 1000 columns with possible SST values.