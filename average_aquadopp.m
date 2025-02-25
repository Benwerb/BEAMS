function ADavg = average_aquadopp(AD, mins_to_avg)

% ADavg = average_aquadopp(AD, mins_to_avg)
%
% Input is structure AD (created from aquadopp2mat), and mins_to_avg, which
% is now many minutes to average. date_range is a vector of size 2, where
% [start_date end_date], in SDN format.
%
% Yui Takeshita
% Scripps Insitution of Oceanography
% Created 3/12/2014
%

% % declare end time
% istart = find(AD.SDN > date_range(1), 1, 'first');
% iend = find(AD.SDN < date_range(2), 1, 'last');

% 
avg_length = 1/24/60*mins_to_avg; % mins to average in SDN.
samp_freq = AD.SDN(3) - AD.SDN(2);
num_avg = floor(avg_length/samp_freq); % number of measurements to average
ncol_new = floor(length(AD.SDN)/num_avg);


% Start building AD structure (averaged)
for vvar = 1:length(AD.vector_vars)
    ADavg.(AD.vector_vars{vvar}) = NaN(1,ncol_new);
end
for mvar = 1:length(AD.matrix_vars)
    ADavg.(AD.matrix_vars{mvar}) = NaN(size(AD.u,1),ncol_new);
end
ADavg.bin_depth = AD.bin_depth;

for i = 1:ncol_new
    ii = (i-1)*num_avg+1; % indext to start averaging
    % Average vector variables
    for vvar = 1:length(AD.vector_vars)
        ADavg.(AD.vector_vars{vvar})(i) = nanmean(AD.(AD.vector_vars{vvar})(ii:ii+num_avg-1));
    end    
    % Average matrix variables
    for mvar = 1:length(AD.matrix_vars)
        ADavg.(AD.matrix_vars{mvar})(:,i) = nanmean(AD.(AD.matrix_vars{mvar})(:,ii:ii+num_avg-1),2);
    end
end

return















