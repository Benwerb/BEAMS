function Ybin= bin_data_to_X_GF(X,Y,Xbin)


% bin_data_to_X_GF(X,Y,Xbin)

% X and Xbin must be sorted. 

% Xbin are the bins. must be sorted. (low to high)
% X is the independent variable for data to be averaged. MUST BE SORTED.
% Y is the dependent variable to be averaged.


% Author: Yui Takeshita
% Scripps Institution of Oceanography
% Created: 12/29/2014

% start binning data to averaged aquadopp date number
dXbin = diff(Xbin);

Ybin = NaN(size(Xbin));
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
        dX = min([dXbin(i-1), dXbin(i)]);
        iavg = inrange(X, [Xbin(i)-(dX/2), Xbin(i)+(dX/2)]);
    end
    Ybin(i) = nanmean(Y(iavg));
end


