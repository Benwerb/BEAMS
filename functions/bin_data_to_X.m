function [Ybin, Ystd, Ybinmedian, nperbin] = bin_data_to_X(X,Y,Xbin)


% bin_data_to_X(X,Y,Xbin)

% X and Xbin must be sorted. 

% Xbin are the bins. must be sorted. (low to high)
% X is the independent variable for data to be averaged. MUST BE SORTED.
% Y is the dependent variable to be averaged.


% Author: Yui Takeshita
% Scripps Institution of Oceanography
% Created: 12/29/2014
%
% 6/2/2022
% added extra output for n counts per bin
% added output for median bins

% start binning data to averaged aquadopp date number
dXbin = diff(Xbin);

Ybin = NaN(size(Xbin));
Ystd = Ybin;
Ybinmedian = Ybin; % added 6/2/2022 YT
nperbin = Ybin; % added 6/2/2022 YT
for i = 1:length(Xbin)
    % get indicies to average
    if(i == 1)
        % first point
        iavg = inrange(X, [Xbin(i)-(dXbin(i)/2), Xbin(i)+(dXbin(i)/2)]);
    elseif(i == length(Xbin))
        % last point
        iavg = inrange(X, [Xbin(i)-(dXbin(i-1)/2), Xbin(i)+(dXbin(i-1)/2)]);
    else
        % everything else
        iavg = inrange(X, [Xbin(i)-(dXbin(i-1)/2), Xbin(i)+(dXbin(i)/2)]);
    end
    Ybin(i) = nanmean(Y(iavg));
    Ystd(i) = nanstd(Y(iavg));
    Ybinmedian(i) = nanmedian(Y(iavg)); % added 6/2/2022 YT
    nperbin(i) = length(find(iavg)); % added 6/2/2022 YT
end


