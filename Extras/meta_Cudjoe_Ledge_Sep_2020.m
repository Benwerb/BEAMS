% meta_bbay_sed_BEAMS.m

% for BEAMS deployment over sediments at Hog Reef

% path where aquadopp file is saved, including file name; no extension
aqdop_path = [strippath(cd,2),'Sensor data\aquadopp\MARBMS04'];
% path for SeapHOx Data
SP_path = [strippath(cd,2),'Sensor data\SeapHOx\U_902BMS_trimmed.txt'];
% path to save finished .mat file
mat_path = cd;
% .mat file name
mat_name = 'Cudjoe_Sep_2020_BEAMS.mat';
plot_folder = 'Plots';
sitename = 'Cudjoe_Sep_2020';

cor_cutoff = 70; % correlation cutoff for aquadopp data
min2avg = 15; % minutes to average for aquadopp
R2cutoff = 0.7; % cutoff for R2 to be used for ustar calculations.
zmin = 0.2;
zmax = 1.1; % maximum height above benthos to use for fit
Sal_est = 36; % estimate for salinity; used in denisty and viscosity
use_def_Sal = false; % if Mcat wasn't working, it will replace all salinity with Sal_est above
daterange = [datenum('06-30-2020 14:59:58') datenum('10 -26-2020 08:29:31')];
% pumpz = [0.75 0.50];
pumpz = [0.79, 0.33];
U0range = [0.002 0.02]; % velocity range to calculate NEP/NEC
ipumps = [1 2]; %pumps to use to calculate NEP & NEC for Mcgillis 2011
PAR_bin = [0:25:100 200:100:600]; % for PI curve

Eoint_25 = -0.420214442432741; % based on discrete sample; calib_SP.m
Eoext_25 = -1.437804503682005; % 
TA0 = 2310; 
Q = 0.8:0.1:1.2; % Photosynthetic quotient for TA calculations=============