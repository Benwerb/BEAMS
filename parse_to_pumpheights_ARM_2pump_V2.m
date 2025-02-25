function [SP, SPraw] = parse_to_pumpheights_ARM_2pump_V2(SPraw, daterange, rmvfirst2npoints)

% original function by Yui Takeshita at SIO.
%
% Modified: 3/29/2020 to read in meta data that is included in SPraw.
% Specifically, number of samples per depth, and sampling period so it
% knows how many samples to extract, and how long to average
%
% included rmvfist2npoints input, which is the number of samples to not
% include in the averages, standard deviations, outliers, etc for each
% measurement height. 
%
%
% Yui Takeshita
% MBARI



% added 3/29/2020
nreps = str2double(SPraw.Number_of_Samples_per_depth);
npumps = nanmax(SPraw.pumpnum);
sampinterval = str2double(SPraw.Sampling_period_s);
min2avg = sampinterval * nreps / 60; 


% **** DELETE DATA OUT OF DATE RANGE AND START PUMP 1 AND END PUMP 3 ****
% cutoff data outside of datarange
iSDN = inrange(SPraw.SDN, daterange);
% vars = fieldnames(SPraw);
vars = SPraw.vars;

vars = vertcat(vars, {'DOXY'; 'PAR'; 'pHint'; 'pHext'});

for v = 1:length(vars)  
    SPraw.(vars{v})(~iSDN) = [];
end

% now find first case of pump 1, and last case of pump n, and delete data 
istart = find(SPraw.pumpnum == 1 & SPraw.sampnum == 1, 1, 'first');
iend = find(SPraw.pumpnum == npumps & SPraw.sampnum == nreps, 1, 'last');
for v = 1:length(vars)
    if(iend~=length(SPraw.SDN))
        SPraw.(vars{v})(iend+1:end) = [];
    end
    SPraw.(vars{v})(1:istart-1) = [];
end

% ************************************************************************

% Variables to parse into differnet pump heights

vars = {'SDN', 'pHint', 'pHext', 'TC', 'DOXY', 'PAR', 'PSAL', 'Pres', 'O2SATPER', 'DFET_TC', 'OPT_TC', 'MCAT_TC'};
SPraw.TC = SPraw.OPT_TC;
% SPraw.pH = SPraw.pHint_prelim;

% find indicies for pump start and pump end at each height
ip1start = []; ip1end = [];
ip2start = []; ip2end = [];
ip3start = []; ip3end = [];

for i = 1:length(SPraw.pumpnum) 
    if(SPraw.pumpnum(i) == 1 && SPraw.sampnum(i) == 1)
        ip1start = [ip1start; i];
    elseif(SPraw.pumpnum(i) == 1 && SPraw.sampnum(i) == nreps)
        ip1end = [ip1end; i];
    elseif(SPraw.pumpnum(i) == 2 && SPraw.sampnum(i) == 1)
        ip2start = [ip2start; i];
    elseif(SPraw.pumpnum(i) == 2 && SPraw.sampnum(i) == nreps)
        ip2end = [ip2end; i];
%     elseif(SPraw.pumpnum(i) == 3 && SPraw.sampnum(i) == 1)
%         ip3start = [ip3start; i];
%     elseif(SPraw.pumpnum(i) == 3 && SPraw.sampnum(i) == 20)
%         ip3end = [ip3end; i];
    end
end



% loop through each variable
for v = 1:length(vars)
%     % presize matrix
%     SPraw.([vars{v},'1raw']) = [];
%     SPraw.([vars{v},'2raw']) = [];
% %     SPraw.([vars{v},'3raw']) = [];

    % presize vector variables
    SPraw.([vars{v},'1avg']) = NaN(length(ip1start),1);
    SPraw.([vars{v},'2avg']) = NaN(length(ip1start),1);
    SPraw.([vars{v},'1std']) = NaN(length(ip1start),1);
    SPraw.([vars{v},'2std']) = NaN(length(ip1start),1);
    SPraw.([vars{v},'1nsamp']) = NaN(length(ip1start),1);
    SPraw.([vars{v},'2nsamp']) = NaN(length(ip1start),1);
    SPraw.([vars{v},'1outlier']) = zeros(length(ip1start),1);
    SPraw.([vars{v},'2outlier']) = zeros(length(ip1start),1);
    % matrix for raw data
    SPraw.([vars{v},'1raw']) = NaN(length(ip1start),nreps);
    SPraw.([vars{v},'2raw']) = NaN(length(ip1start),nreps);

    for j = 1:length(ip1start)
        % filter out outliers using Thompson Tau method
        % Omit first 2 points after pump switches
        vec1in = SPraw.(vars{v})(ip1start(j)+rmvfirst2npoints:ip1end(j)); 
        vec2in = SPraw.(vars{v})(ip2start(j)+rmvfirst2npoints:ip2end(j));
%         vec3in = SPraw.(vars{v})(ip3start(j)+2:ip3end(j)); 
        vec1out = removeoutliers(vec1in); 
        vec2out = removeoutliers(vec2in);
%         vec3out = removeoutliers(vec3in);
        % check to see if an outlier was removed
        if(length(vec1out)~= length(vec1in))
            SPraw.([vars{v},'1outlier'])(j) = 1;
        end
        if(length(vec2out)~= length(vec2in))
            SPraw.([vars{v},'2outlier'])(j) = 1;
        end
%         if(length(vec3out)~= length(vec3in))
%             SPraw.([vars{v},'3outlier'])(j) = 1;
%         end        

        % save all raw data into m x n matrix, where m = length of each new
        % meausrement set, n = reps at each height 3/29/2020
        SPraw.([vars{v},'1raw'])(j,:) = SPraw.(vars{v})(ip1start(j):ip1end(j));
        SPraw.([vars{v},'2raw'])(j,:) = SPraw.(vars{v})(ip2start(j):ip2end(j));

        % save raw data with first 2 points removed
%         SPraw.([vars{v},'1raw']) = [SPraw.([vars{v},'1raw']); vec1in];
%         SPraw.([vars{v},'2raw']) = [SPraw.([vars{v},'2raw']); vec2in];   
%         SPraw.([vars{v},'3raw']) = [SPraw.([vars{v},'3raw']); vec3in];
        % save average of variable, and its std. dev.
        SPraw.([vars{v},'1avg'])(j) = nanmean(vec1out);
        SPraw.([vars{v},'2avg'])(j) = nanmean(vec2out);
%         SPraw.([vars{v},'3avg'])(j) = nanmean(vec3out);
        SPraw.([vars{v},'1std'])(j) = nanstd(vec1out);
        SPraw.([vars{v},'2std'])(j) = nanstd(vec2out);
%         SPraw.([vars{v},'3std'])(j) = nanstd(vec3out);
        SPraw.([vars{v},'1nsamp'])(j) = length(vec1out);
        SPraw.([vars{v},'2nsamp'])(j) = length(vec2out);
%         SPraw.([vars{v},'3nsamp'])(j) = length(vec3out);
    end
end

% now interpolate through
% X variable to interpolate onto
% SPraw.SDNint = SPraw.SDN(1):datenum(0,0,0,0,min2avg,0):SPraw.SDN(end);
% instead of clean 0:15:30:45 minutes, use the actual measurement time for
% better accuracy. 3/29/2020
SPraw.SDNint = sort(vertcat(SPraw.SDN1avg, SPraw.SDN2avg))';
interp_method = 'linear';

for v = 1:length(vars)
    SPraw.([vars{v},'1int']) = interp1(SPraw.SDN1avg, SPraw.([vars{v},'1avg']), SPraw.SDNint, interp_method);
    SPraw.([vars{v},'2int']) = interp1(SPraw.SDN2avg, SPraw.([vars{v},'2avg']), SPraw.SDNint, interp_method);
%     SPraw.([vars{v},'3int']) = interp1(SPraw.SDN3avg, SPraw.([vars{v},'3avg']), SPraw.SDNint, interp_method);         
end

% now make SP structure, where variable are 3xn matricies

% extract desired variables into matrix

for v = 2:length(vars)
    topvar = [vars{v},'1int'];
%     midvar = [vars{v},'2int'];
    botvar = [vars{v},'2int'];

%     SP.(vars{v}) = [SPraw.(topvar)(:)'; vec2row(NaN(size(SPraw.(topvar)))); SPraw.(botvar)(:)'];    
    SP.(vars{v}) = [SPraw.(topvar)(:)'; SPraw.(botvar)(:)'];    
    
    % now put together std data
    SP.([vars{v},'std']) = NaN(size(SPraw.SDNint));
    
    SP.([vars{v},'std'])(1:2:end) = SPraw.([vars{v},'1std']);
    SP.([vars{v},'std'])(2:2:end) = SPraw.([vars{v},'2std']);
end

SP.SDN = SPraw.SDNint;


return        