%% DEPENDENCIES
% Clone the following repositories and have in MATLAB path: 
% https://github.com/CarbonLab/import_functions/tree/main/lab_instruments/pHETSOS
% https://github.com/yui-takeshita/seawater_equilib

% Download the following functions and have in MATLAB path
% https://www.marine.csiro.au/datacentre/projects/seawater/seawater_ver1_2d.zip
% https://www.mathworks.com/matlabcentral/fileexchange/23441-polar-cartesian-into-degree-north-reference

%%
clear all;
close all;

% load metadata for specific deployment

% -------------------------------------------------------------------------
% % Hog sediment deployment Bermuda July 2015
% run([cd,'\Bermuda\2015\Andersson_July\mfiles for analysis\Calibrated and processed BEAMS data\meta_sed_hog_BEAMS.m'])
% -------------------------------------------------------------------------
% Hog Reef deployment Bermuda July 2015
% run([cd,'\Bermuda\2015\Andersson_July\mfiles for analysis\Calibrated and processed BEAMS data\meta_hog_reef_BEAMS.m'])
% -------------------------------------------------------------------------
% % Baileys Bay Sediment Deployment Bermuda July 2015
% run([cd,'\Bermuda\2015\Andersson_July\mfiles for analysis\Calibrated and processed BEAMS data\meta_bbay_sed_BEAMS.m']) 
% -------------------------------------------------------------------------
% % Baileys Bay Seagrass (alongside EC) deplyment, Bermuda July 2015
% % DON'T RUN THIS YET. NO AQUADOPP DATA; GOT USTAR FROM WADE.
% run([cd,'\Bermuda\2015\Andersson_July\mfiles for analysis\Calibrated and processed BEAMS data\meta_bbay_sgrss_BEAMS.m'])
% -------------------------------------------------------------------------
% Mote Coral Nursey Deployment, Aug 2018 First Deploymeht
% run([cd,'\Coral_nursery_Platz\2018_Aug\mfiles for analysis\Calibrated_BEAMS_data_first_dep\meta_bbay_Mote_Aug_2018_firstdep.m']);
% -------------------------------------------------------------------------
% % Mote Coral Nursey Deployment, Aug 2018 Second Deployment
% run([cd,'\Coral_nursery_Platz\2018_Aug\mfiles for analysis\Calibrated_BEAMS_data\meta_bbay_Mote_Aug_2018.m']);
% -------------------------------------------------------------------------
% Mote Coral Nursey Deployment, March 2019 Second Deployment
% run([cd,'\Coral_nursery_Platz\2019_Mar\mfiles for analysis\Calibrated_BEAMS_data\meta_bbay_Mote_Mar_2019.m']);
% -------------------------------------------------------------------------
% Mote Coral Nursey Deployment, March 2019 Second Deployment
% run([cd,'\Coral_nursery_Platz\2020_Jul\Cudjoe Ledge_\mfiles for analysis\Calibrated_BEAMS_data\meta_Cudjoe_Ledge_Jul_2020.m']);
% -------------------------------------------------------------------------
% Mote Coral Nursey Deployment, Jul 2020, Marker 32
% run([cd,'\Coral_nursery_Platz\2020_Jul\Marker 32\mfiles for analysis\Calibrated_BEAMS_data\meta_Marker32_Jul_2020.m']);
% -------------------------------------------------------------------------
% Mote Coral Nursey Deployment, Sep 2020, Marker 32
% run([cd,'\Coral_nursery_Platz\2020_Sep\Marker 32\mfiles for analysis\Calibrated_BEAMS_data\meta_Marker32_Sep_2020.m']);
% -------------------------------------------------------------------------
% Mote Coral Nursey Deployment, Sep 2020, Marker32
% run([cd,'\Coral_nursery_Platz\2020_Sep\Marker 32\mfiles for analysis\Calibrated_BEAMS_data\meta_Marker32_Sep_2020.m']);
% -------------------------------------------------------------------------
% Mote Coral Nursey Deployment, Sep 2020, Cudjoe
% run([cd,'\Coral_nursery_Platz\2020_Sep\Cudjoe Ledge\mfiles for analysis\Calibrated_BEAMS_data\meta_Cudjoe_Ledge_Sep_2020.m']);
% -------------------------------------------------------------------------
% B1_Cheeca_26June2024
% run([cd, 'C:\Users\bwerb\Documents\GitHub\BEAMS\testdata'])

%% Meta information
% Filepaths and save locations
% Base directory  
base_dir = 'G:\.shortcut-targets-by-id\18QB8oWsbcr1m5sL5U67c2irRD1encSYU\Deployments_2024\B1_Cheeca_26June2024';

% Filepaths and save locations  
SPfilepath = fullfile(base_dir, 'BEAMS', '20240626.txt');  
aqdop_path = fullfile(base_dir, 'ADP', '626A1305.csv');
% localsavepath = 'C:\Users\bwerb\Documents\GitHub\BEAMS\testdata';
mat_name = fullfile(base_dir, 'B1_Cheeca_26June2024B1_BEAMS.mat');  


% SP Extract info
use_def_Sal = false; % defaulting to false
startdate = datetime(2024,6,26,14,30,0,TimeZone="UTC");
enddate = datetime(2024,7,2,14,0,0,TimeZone="UTC");
daterange = [datenum(startdate) datenum(enddate)];
% daterange = SPraw.SDN(1):SPraw.SDN(end); % Default for now, change later!

% Calc TA Gradient info
TA0 = 2250; % Assumed background total alkalinity for calc_TA_gradient
Q = 0.8:0.1:1.2; % Photosynthetic quotient for calc_TA_gradient

% ustar_from_aquadopp inital estimates
    % Just an estimate to give the model a starting place.  
    % This could be from bottle samples or meta data. Use if Mcat fails.
Sal_est = 36; 

% Aquadopp info
hasAquadopp = true; % Aquadopp data included? True or False
delimit = ';';
min2avg = 15; % minutes to average for aquapdopp

% BEAMS struct info
pumpz = [0.79, 0.33]; % Where can this data be found in meta?
ipumps = [1 2]; % pumps to use to calculate NEP & NEC for Mcgillis 2011

% % U0range = [0.002 0.02]; % velocity range to calculate NEP/NEC
%% Extract and parse SeapHOx data (ARM board)
disp('Start Extracting SeapHOx data...');
GetMfetData(SPfilepath);
% Change variable names
SPraw = data; clearvars data; % rename data to SPraw

DeploymentInfo = meta_parser(meta); % Parse k0, k2, and pump info from meta

% Match names with old parser for continuity
oldNames = {'TC_cont', 'pSal', 'TC_opt', 'TC_MCat', 'O2satper'};
newNames = {'TC', 'PSAL', 'OPT_TC', 'MCAT_TC', 'O2SATPER'};
SPraw = renamevars(SPraw, oldNames, newNames); 

% Add SDN column
SPraw.SDN = datenum(SPraw.("MM/DD/YYYY HH:MM:SS")); % BEAMS should always output column in this format using GetMFETData

% Remove '#' and convert to double
SPraw.SampNum = str2double(erase(SPraw.SampNum, '#'));

% Split Pump into multiple columns and extract values
tokens = regexp(SPraw.Pump, 'P(\d+) (\d+)/(\d+)', 'tokens');

% Convert cell arrays into numeric arrays
Pump = cellfun(@(x) str2double(x{1}{1}), tokens);        % Extract pump number
CycleNumber = cellfun(@(x) str2double(x{1}{2}), tokens); % Extract cycle number
Number_of_Samples_per_depth = cellfun(@(x) str2double(x{1}{3}), tokens);      % Extract total cycles

% Add parsed data as new columns in the table
SPraw.pumpnum = Pump;
SPraw.CycleNumber = CycleNumber;
SPraw.Number_of_Samples_per_depth = Number_of_Samples_per_depth;

% Find number of pumps from parsed data
DeploymentInfo.npumps = nanmax(SPraw.pumpnum); 

% replace psal with default salinity if chosen 
if use_def_Sal == true
    SPraw.PSAL(:) = Sal_est;
end

% calculate oxygen in umol/kg
SPraw.DOXY = calcO2sat(SPraw.OPT_TC, SPraw.PSAL).*SPraw.O2SATPER./100;

% Calculate pHext total from vrse, temperature, salinity, k0, and k2
SPraw.pHext = calc_pHext(SPraw.Vrse, SPraw.TC, SPraw.PSAL, DeploymentInfo.k0ext, DeploymentInfo.k2ext);

disp('Parsing data to different pump heights, calcualting TA...');
% parse into 2 pump heights
[SP, SPraw] = parse_to_pumpheights_ARM_2pump_V2(SPraw, DeploymentInfo, daterange, 4);
% pHext is the vrse we get from the DSD
SP.pH = SP.pHext;

% Calculate the TA gradient
SP = calc_TA_gradient_V2(SP, TA0, Q, DeploymentInfo.npumps);

% Save files
if exist(mat_name, 'file')
    save(mat_name, 'SPraw', 'SP', '-append');  % Append if file exists
else
    save(mat_name, 'SPraw', 'SP');  % Create new file if it doesn't exist
end

%% Averaging Aquadopp and ustar calculations start here
% ******************* Extract data ******************************

% extract aquadopp from raw correlation and uvw data
disp('Opening Raw Aquadopp data...');
AD = aquadopp2mat(aqdop_path, delimit, daterange);
disp('Raw Aquadopp Data extracted');
disp('% of data removed; 1st row = % cutoff, 2,3,4th row is u, v, w');

% [AD.cor_cutoffs; AD.per_reject]
% ***************************************************************
% avrerage aquadopp data 

% Average the aquadop data
ADavg = average_aquadopp(AD, min2avg);
disp('Aquadopp data averaged');

% ***************************************************************
% Is there a reason to predifine these values or can I just use:
zmin = min(ADavg.bin_depth); % BW
zmax = max(ADavg.bin_depth); % BW
if(hasAquadopp)
    ADavg = ustar_from_aquadopp(ADavg, [zmin zmax], Sal_est);
    % % save ADavg.
else
    % need Drag Coefficient, and point velocity measurement.
    
end

% save AD and ADavg to matfile
if exist(mat_name, 'file')
    save(mat_name, 'AD', 'ADavg', '-append');  % Append if file exists
else
    save(mat_name, 'AD', 'ADavg');  % Create new file if it doesn't exist
end

disp('Friction velocity calcualted; Aqaudopp data saved');

%% Create BMS and BMSbin structures, and calculate fluxes (NEC&NEP)

% ==================Create BMS============================
disp('Calculating NEP and NEC...');
BMS = SP;
% X variable to bin Aquadopp data to.
X = BMS.SDN;
% Bin aquadopp vector variables (Pressure, U0, direction, ustar; U0 and
% direction are measured at 1m, extracted in ustar_from_aquadopp.m)
BMS.PRES_AD  = bin_data_to_X(ADavg.SDN, ADavg.Pres, X);
BMS.U0       = bin_data_to_X(ADavg.SDN, ADavg.U0, X); % BW commented
% BMS.DIR      = bin_data_to_X(ADavg.SDN, ADavg.direction1m, X); BW commented
BMS.ustar    = bin_data_to_X(ADavg.SDN, ADavg.ustar, X);
% BMS.ustar    = bin_data_to_X(ADavg.SDN, ADavg.uv(3,:).*0.12, X)'; %for Hog Reef BEAMS
BMS.ustar_rm = bin_data_to_X(ADavg.SDN, ADavg.ustar_runmean, X);
% calculate density
BMS.DENS = sw_dens(SP.PSAL, SP.TC, BMS.PRES_AD);
BMS.pumpz = pumpz;
% BMS.ustar_rm = BMS.U0.*0.12;


% calcualte NEP & NEC
BMS = calc_NEP_NEC_BEAMS(BMS, ipumps, 'ustar');

% some manual QC; delete data in date ranges listed in meta file
% still under construction.
% if(~isempty(bad_data_SDN))
%     % still need to code it
%     
% end

disp('Binning BEAMS data...');
% ***********Create BMSbin**************************************
% X variable to bin data to
X = floor(nanmin(BMS.SDN)):1/24:ceil(nanmax(BMS.SDN));

BMSbin.SDN = X;
% bin data hourly. Vector variables first
BMSbin.PRES_AD  = bin_data_to_X_GF(BMS.SDN, BMS.PRES_AD, X);
BMSbin.U0       = bin_data_to_X_GF(BMS.SDN, BMS.U0, X);
% BMSbin.DIR      = bin_data_to_X_GF(BMS.SDN, BMS.DIR, X); % BW commented
BMSbin.ustar    = bin_data_to_X_GF(BMS.SDN, BMS.ustar, X);
BMSbin.ustar_rm = bin_data_to_X_GF(BMS.SDN, BMS.ustar_rm, X);
BMSbin.NEP      = bin_data_to_X_GF(BMS.SDN, BMS.NEP, X);
BMSbin.NEP_WM   = bin_data_to_X_GF(BMS.SDN, BMS.NEP_WM, X);

% variables from 2 differnet heights. 3 is commented out
BMSbin.TC(1,:)       = bin_data_to_X_GF(BMS.SDN, BMS.TC(1,:), X);
BMSbin.TC(2,:)       = bin_data_to_X_GF(BMS.SDN, BMS.TC(2,:), X);
% BMSbin.TC(3,:)       = bin_data_to_X_GF(BMS.SDN, BMS.TC(3,:), X);
BMSbin.DOXY(1,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DOXY(1,:), X);
BMSbin.DOXY(2,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DOXY(2,:), X);
% BMSbin.DOXY(3,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DOXY(3,:), X);
BMSbin.pH(1,:)       = bin_data_to_X_GF(BMS.SDN, BMS.pH(1,:), X);
BMSbin.pH(2,:)       = bin_data_to_X_GF(BMS.SDN, BMS.pH(2,:), X);
% BMSbin.pH(3,:)       = bin_data_to_X_GF(BMS.SDN, BMS.pH(3,:), X);
BMSbin.PSAL(1,:)     = bin_data_to_X_GF(BMS.SDN, BMS.PSAL(1,:), X);
BMSbin.PSAL(2,:)     = bin_data_to_X_GF(BMS.SDN, BMS.PSAL(2,:), X);
% BMSbin.PSAL(3,:)     = bin_data_to_X_GF(BMS.SDN, BMS.PSAL(3,:), X);
BMSbin.O2SATPER(1,:) = bin_data_to_X_GF(BMS.SDN, BMS.O2SATPER(1,:), X);
BMSbin.O2SATPER(2,:) = bin_data_to_X_GF(BMS.SDN, BMS.O2SATPER(2,:), X);
% BMSbin.O2SATPER(3,:) = bin_data_to_X_GF(BMS.SDN, BMS.O2SATPER(3,:), X);
BMSbin.Pres(1,:)     = bin_data_to_X_GF(BMS.SDN, BMS.Pres(1,:), X);
BMSbin.Pres(2,:)     = bin_data_to_X_GF(BMS.SDN, BMS.Pres(2,:), X);
% BMSbin.Pres(3,:)     = bin_data_to_X_GF(BMS.SDN, BMS.Pres(3,:), X);  
BMSbin.DENS(1,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DENS(1,:), X);
BMSbin.DENS(2,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DENS(2,:), X);
% BMSbin.DENS(3,:)     = bin_data_to_X_GF(BMS.SDN, BMS.DENS(3,:), X);  
% No PAR data for AOML
% BMSbin.PAR(1,:)      = bin_data_to_X_GF(BMS.SDN, BMS.PAR(1,:), X); % BW no PAR
% BMSbin.PAR(2,:)      = bin_data_to_X_GF(BMS.SDN, BMS.PAR(2,:), X); % BW no PAR
% BMSbin.PAR(3,:)      = bin_data_to_X_GF(BMS.SDN, BMS.PAR(3,:), X);  

% Now TA and NEC
for ii = 1:size(BMS.TA,1)
   BMSbin.TA(ii,1,:)    = bin_data_to_X_GF(BMS.SDN, squeeze(BMS.TA(ii,1,:)), X); 
   BMSbin.TA(ii,2,:)    = bin_data_to_X_GF(BMS.SDN, squeeze(BMS.TA(ii,2,:)), X); 
%    BMSbin.TA(ii,3,:)    = bin_data_to_X_GF(BMS.SDN, squeeze(BMS.TA(ii,3,:)), X);      
   BMSbin.NEC(ii,:)     = bin_data_to_X_GF(BMS.SDN, BMS.NEC(ii,:), X);
   BMSbin.NEC_WM(ii,:)  = bin_data_to_X_GF(BMS.SDN, BMS.NEC_WM(ii,:), X);
end
disp('Finished Metabolism Calculations');

if exist(mat_name, 'file')
    save(mat_name, 'BMS', 'BMSbin', '-append');  % Append if file exists
else
    save(mat_name, 'BMS', 'BMSbin');  % Create new file if it doesn't exist
end

% save(mat_name, 'BMS', 'BMSbin', '-append');