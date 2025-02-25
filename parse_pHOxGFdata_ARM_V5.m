function [SP, varnames] = parse_pHOxGFdata_ARM_V5(filename)

% [SP_data, varnames] = parse_SPdata_SF_V1(filename)
%
% Parses SeapHOx data for program version SeaFET_Controllver_V2. The text file
% has pH estimates on the fly, however, does not have TC calculated from 
% Vtherm yet. This assumes oxygen, and Microcat data=a is there as well. 
%
% Yui Takeshita
% Scripps Institution of Oceanography
% Created 10/23/2013




varnames = {'SDN'
    'VBATT_ON'
    'VBATT_OFF'
    'Vtherm'
    'Vint'
    'Vext'
    'Vk'
    'Visobatt'
    'Cont_TC'
    'DFET_TC'
    'Vpar'
    'pHint_prelim'
    'pHext_prelim'
    'Ik'
    'Ib'
    'pumpnum'
    'sampnum'
    'OPT_MN'
    'OPT_SN'
    'M_DOXY'
    'O2SATPER'
    'OPT_TC'
    'Dphase'
    'Bphase'    
    'MCAT_TC'
    'PSAL'
    'Pres'
    'MCAT_SDN'};

% get file ID
fid = fopen(filename);


% extract meta data
fstr = fgetl(fid);

% look for empty line; that's where data starts
while(~isempty(fstr))
    
    % bandaid solution to ignore last header which doesn't have '\t'
    if(~strncmpi(fstr,'"Deploy mode',12))
        % parse it into two strings, using tab as delimiter
        trex = textscan(fstr,'%s%s', 'delimiter', '\t');
        % extract meta data name and value
        metaname = trex{1}{1};
        metaval  = trex{2}{1};
        clear trex
        % remove " " for metaname
        metaname(strfind(metaname,'"')) = [];
        metaname(strfind(metaname,'(')) = [];
        metaname(strfind(metaname,')')) = [];
        metaname(strfind(metaname,' ')) = '_';
        metaname(strfind(metaname,'-')) = '_';
        SP.(metaname) = metaval;
    end
    % read next line
    fstr = fgetl(fid);
    
end

[trex] = textscan(fid, ['%s%s%f%f%f%f%f%f%f%f',...
                      '%f%f%f%f%f%f%f%f%f%f',...
                      '%f%f%f%f%f%f%f%f%f%f', ...
                      '%f%f%f%f%f%s%s'], ...
                      'delimiter', '\t', ...
                      'collectoutput', 1, ...
                      'headerlines', 1);           
% Close file ID
fclose(fid);   



% Declare variables
for var = 1:length(varnames)
    SP.(varnames{var}) = [];
end


% Divide output of textscan
% columns in data are floats
SDN = trex{1}; SDN = SDN(:,2);
data = trex{2};
MCATSDN = trex{3};
% ================ build data structure=====================
% ADC data
SP.SDN             = datenum(SDN);
SP.VBATT_ON        = data(:,1);
SP.VBATT_OFF       = data(:,2);
SP.Vtherm          = data(:,3);
SP.Vint            = data(:,4);
SP.Vext            = data(:,5);
SP.Vk              = data(:,6);
SP.Visobatt        = data(:,7);
SP.Cont_TC         = data(:,8);
SP.DFET_TC         = data(:,9);
SP.Vpar            = data(:,11);
SP.pHint_prelim    = data(:,12);
SP.pHext_prelim    = data(:,13);
SP.pumpnum         = data(:,14);
SP.sampnum         = data(:,15);
SP.Ik              = data(:,16)./0.00015;
SP.Ib              = data(:,17)./0.001;
% optode data
SP.OPT_MN          = data(:,18);
SP.OPT_SN          = data(:,19);
SP.M_DOXY          = data(:,20);
SP.O2SATPER        = data(:,21);
SP.OPT_TC          = data(:,22);
SP.Dphase          = data(:,23);
SP.Bphase          = data(:,24);
SP.MCAT_TC         = data(:,30);
SP.Pres            = data(:,32); % pressure sensor on SBE37. 
SP.PSAL            = data(:,33);
% see if MCAT SDN has both date and time
% make array of NaNs first. in case there is no MCAT data due to not being
% plugged in, or malfunction.

% presize MCAT_SDN vector
SP.MCAT_SDN = NaN(size(SP.MCAT_TC));
% find the NANs.
inan = strcmp(MCATSDN(:,1),'NAN');
% if array is not all NANs (the case for a functioning microcat), then use
% datenum
if(~isempty(find(~inan,1)))
    SP.MCAT_SDN(~inan) = datenum(MCATSDN(~inan,1)) + mod(datenum(MCATSDN(~inan,2)),1);
end
% Calculate PAR based on factory calibration; 1mV = 5 umol m-2  s-1
SP.PAR = SP.Vpar.*1000.*5; 
% write variable names
SP.vars = varnames;

return
