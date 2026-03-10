% Setup Matlab
clear; close all; clc;
addpath("../GMM_estimation_v3")

% Setup and low Colombia data
pctiles_raw     = csvread("../../Data/WorldBankData/Clean data/colombia_exp_percentiles.csv",1,0);
county_names    = readtable("../../Data/WorldBankData/Clean data/colombia_exp_names.csv");

% Initialize results
r2              = zeros(10,4);
B              = zeros(10,4);
SE              = zeros(10,4);
Est_baseline    = zeros(10,2);

% Loop over top 10 countries
for i = 1:10
    if i == 1
        plot = 1;
    else
        plot = 0;
    end
    [r2i, est_basei, B_vec, SE_vec] = QuantilesColombia(pctiles_raw,i,plot);
    r2(i,:) = r2i;
    B(i,:) = B_vec;
    SE(i,:) = SE_vec;
    Est_baseline(i,:) = est_basei;
end

mean_r2 = mean(r2);

% Display results
fprintf('%10s:\t    Baseline\tLog-Normal\tLog-Normal(selection)\t    Pareto\n', "Country/R2");
for i = 1:10
    fprintf(' %10s:\t %10.3f\t%10.3f\t%10.3f\t%10.3f\n', string(county_names{i, 1}), r2(i, 1), r2(i, 2), r2(i, 4), r2(i, 3));
end
fprintf(' %10s:\t %10.3f\t%10.3f\t%10.3f\t%10.3f\n', "Mean",mean_r2(1), mean_r2(2), mean_r2(4), mean_r2(3));


% Display results
fprintf('%10s:\t    Baseline\tLog-Normal\tLog-Normal(selection)\t    Pareto\n', "Country/B(se)");
for i = 1:10
    fprintf(' %10s:\t %10.3f\t%10.3f\t%10.3f\t%10.3f\n', string(county_names{i, 1}), B(i, 1), B(i, 2), B(i, 4), B(i, 3));
    fprintf(' %10s:\t (%10.3f)\t(%10.3f)\t(%10.3f)\t(%10.3f)\n', string(county_names{i, 1}), SE(i, 1), SE(i, 2), SE(i, 4), SE(i, 3));
end
% fprintf(' %10s:\t %10.3f\t%10.3f\t%10.3f\t%10.3f\n', "Mean",mean_r2(1), mean_r2(2), mean_r2(4), mean_r2(3));





% Write results to CSV
output_filename = '../../Output/QQ_Colombia_R2_Summary.csv';
fileID = fopen(output_filename, 'w');
fprintf(fileID, '"Country/R2","Baseline","Log-Normal","Log-Normal (Selection)","Pareto"\n');
for i = 1:size(r2, 1)
    fprintf(fileID, '"%s",%f,%f,%f,%f\n', string(county_names{i, 1}), r2(i, 1), r2(i, 2), r2(i, 4), r2(i, 3));
end
fprintf(fileID, '"Mean",%f,%f,%f,%f\n', mean_r2(1), mean_r2(2), mean_r2(4), mean_r2(3));
fclose(fileID);
