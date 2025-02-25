% rawAD2BEAMS.m

% template file to go from raw Aquadopp data to mat file used for BEAMS
% calculations.

% Extracts velocity field data from Aquadopp from raw aquadopp data. Uses a
% correlation cutoff. 
% 
% 



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
run([cd,'\Coral_nursery_Platz\2020_Sep\Cudjoe Ledge\mfiles for analysis\Calibrated_BEAMS_data\meta_Cudjoe_Ledge_Sep_2020.m']);




hasAquadopp = true;

%% Averaging Aquadopp and ustar calculations start here

% ******************* Extract data ******************************

% extract aquadopp from raw correlation and uvw data
disp('Opening Raw Aquadopp data...');
AD = aquadoppraw2mat(aqdop_path, cor_cutoff, daterange);
disp('Raw Aquadopp Data extracted');
disp('% of data removed; 1st row = % cutoff, 2,3,4th row is u, v, w');
[AD.cor_cutoffs; AD.per_reject]
% ***************************************************************
% avrerage aquadopp data 
ADavg = average_aquadopp(AD, min2avg);
disp('Aquadopp data averaged');
% save AD and ADavg to matfile
cd(mat_path);
save(mat_name, 'AD', 'ADavg');
% clear AD to free up some memory
clear AD;
% ***************************************************************
if(hasAquadopp)
    ADavg = ustar_from_aquadopp(ADavg, [zmin zmax], Sal_est);
    % % save ADavg.
else
    % need Drag Coefficient, and point velocity measurement.
    
end

    
cd(mat_path);


save(mat_name, 'ADavg', '-append');

disp('Friction velocity calcualted; Aqaudopp data saved');

%% Extract and parse SeapHOx data (ARM board)

disp('Start Extracting SeapHOx data...');
% for Kaneohe Bay deployment, 3 pumps (LL) or 2 pumps (V1)
% SPraw = parse_pHoxGFdata_V1_LL(SP_path);
% % for 2 pumps, 
% SPraw = parse_pHoxGFdata_V1(SP_path);
% ARM_V1 includes PAR as well. Use for deployments after July 2015
SPraw = parse_pHOxGFdata_ARM_V5(SP_path);

% replace psal with default salinity if chosen 
if(use_def_Sal)
    SPraw.PSAL(:) = Sal_est;
end
% calculate oxygen in umol/kg
SPraw.DOXY = calcO2sat(SPraw.OPT_TC, SPraw.PSAL).*SPraw.O2SATPER./100;
% Calculate pH from internal reference
SPraw.pHint = calc_dfet_pHint(SPraw.Vint, SPraw.DFET_TC, Eoint_25);
SPraw.pHext = calc_dfet_pHext(SPraw.Vext, SPraw.DFET_TC, SPraw.PSAL, Eoext_25);

disp('Parsing data to different pump heights, calcualting TA...');
% parse into 3 pump heights
[SP, SPraw] = parse_to_pumpheights_ARM_2pump_V2(SPraw, daterange,4);
% select int or ext to use as pH for TA and NCC calculations. 
SP.pH = SP.pHint;
SP = calc_TA_gradient_V2(SP, TA0, Q, 2);
cd(mat_path);
save(mat_name, 'SPraw', 'SP', '-append');


%% Create BMS and BMSbin structures, and calculate fluxes (NEC&NEP)

% ==================Create BMS============================
disp('Calculating NEP and NEC...');
BMS = SP;
% X variable to bin Aquadopp data to.
X = BMS.SDN;
% Bin aquadopp vector variables (Pressure, U0, direction, ustar; U0 and
% direction are measured at 1m, extracted in ustar_from_aquadopp.m)
BMS.PRES_AD  = bin_data_to_X(ADavg.SDN, ADavg.Pres, X);
BMS.U0       = bin_data_to_X(ADavg.SDN, ADavg.U0, X);
BMS.DIR      = bin_data_to_X(ADavg.SDN, ADavg.direction1m, X);
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
BMSbin.DIR      = bin_data_to_X_GF(BMS.SDN, BMS.DIR, X);
BMSbin.ustar    = bin_data_to_X_GF(BMS.SDN, BMS.ustar, X);
BMSbin.ustar_rm = bin_data_to_X_GF(BMS.SDN, BMS.ustar_rm, X);
BMSbin.NEP      = bin_data_to_X_GF(BMS.SDN, BMS.NEP, X);
BMSbin.NEP_WM   = bin_data_to_X_GF(BMS.SDN, BMS.NEP_WM, X);

% variables from 3 differnet heights
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
BMSbin.PAR(1,:)      = bin_data_to_X_GF(BMS.SDN, BMS.PAR(1,:), X);
BMSbin.PAR(2,:)      = bin_data_to_X_GF(BMS.SDN, BMS.PAR(2,:), X);
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

% save('Hog_reef_BEAMS_ustar2U0_012.mat', 'BMS', 'BMSbin');

save(mat_name, 'BMS', 'BMSbin', '-append');











