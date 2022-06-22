% Configuration for noise study

% Set which parts of the script to run ------------------------------------
system_definition = false;
DeePC = false;
plotting = true;

% Values for parameter search ---------------------------------------------
seeds = [0 1; 2 3; 4 5; 6 7];
noise_levels = [0.00075 0.0015 0.003 0.006 0.012 0.024];
lambda_gs = [0 3 5 7 10 15 20 30 40 50 60 100 150 200 250 300 400 500 600 700 800];
% seeds = [0 1; 0 1];
% noise_levels = [0.00075 0.00075];
% lambda_gs = [10 10];

% file names
file_save = 'errors3.mat';
file_load = 'errors3.mat';