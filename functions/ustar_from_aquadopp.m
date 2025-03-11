function [ADavg] = ustar_from_aquadopp(ADavg, zrange, Sal_est)

% [ADavg] = ustar_from_aquadopp(ADavg, zrange, Sal_est);
%
% Calculates friction velocity from current profile measured by HR
% Aquadopp. ADavg is the structure containing Aquadopp data, averaged over
% X minutes (created by average_aquadopp.m), zrange is the depth range to
% use for the ustar calculation [zmin zmax], and Sal_est is the estimate
% for salinity (used for viscosity calculation). 
%
% Calculates drag coefficient based on U0 at 1 m, running average of ustar,
% and at the end viscosity and Reynolds number. Calculated values are saved
% into the ADavg structure.
%
% Created by: Yui Takeshita
% MBARI
% Version Dec-6-2018


% calculate ustar by fitting ustar, d, and z0 for 

% inital guesses for fit. ustar is estimated from an drag coefficient from
% Takeshita et al. 2016. It may need to be adjusted for less turbulent
% systems. 
dguess = 0.1;
z0guess = 0.001;
% a = [ustarguess, dguess, z0guess];
lb = [0, 0, 0.0000001]; % lower bound for ustar, d, and z0. 
ub = [1, 0.8, 0.03];    % upper bound for ustar, d, and z0.
iz = inrange(ADavg.bin_depth, zrange); % index for depth range to use for fit
% index for 1 m above benthos. This will be used for U0. 
i1m = find(ADavg.bin_depth <= 1, 1, 'last');    

% presize matrix
ADavg.ustar = NaN(size(ADavg.SDN));
ADavg.d = ADavg.ustar;
ADavg.z0 = ADavg.ustar;
ADavg.vel_fit = NaN(size(ADavg.uv)); % velocity profile for fit

% loop through each profile
for i = 1:length(ADavg.SDN)
    
    % use data points that are not NaNs.
    inonan = ~isnan(ADavg.uv(:,i));
    % if current profile has more than 3 points in the desired range, do fit:
    if(length(find(inonan&iz))>5)

        % ommit lowest point; large error
        inonan(1) = false;
        vel = ADavg.uv(inonan&iz,i);
        z = ADavg.bin_depth(inonan&iz);
        
        ustarguess = nanmean(ADavg.uv(i1m-5:i1m+5,i))*0.06; % sqrt(0.019), which is drag coefficient from Takeshita et al. 2016
        if(isnan(ustarguess))
            ustarguess = 0.01;
        end
        
        % *****************************************************************
        % USE LSQCURVEFIT IF YOU HAVE OPTIMIZATION TOOLBOX (MBARI DOES NOT)
        % *****************************************************************
        
        a = [ustarguess, dguess, z0guess];
        % use lsqcurvefit to fit current profile data to law of the wall.
        options = optimset('Display','none', 'diagnostics', 'off');
        % can't use lsqcurvefit because MBARI doesn't have optimization
        % toolbox
                [c, resnorm, resid, ~]  = lsqcurvefit(@lawofwall, a, z, vel, lb, ub, options);

        % *****************************************************************
        % USE LMFSOLVE INSTEAD IF YOU DON'T HAVE OPTIMIZATION TOOLBOX
        % *****************************************************************
        %
        % LMFsolve can be downloaded from:
        % https://www.mathworks.com/matlabcentral/fileexchange/16063-lmfsolve-m-levenberg-marquardt-fletcher-algorithm-for-nonlinear-least-squares-problems
        % !!!!!!!!!! only 14 ratings; not thoroughly tested !!!!!!!!!!!!
        % 
        % My attempt at creating a cost fxn when d < 0. Can't set lower
        % bound for LMFsolve, so this was an issue. nt an issue for
        % lsqcurvefit. 
%         p = @(x) z-x(2);
% %         lawall = @(x) vel - x(1)./0.41.*log((z-x(2))./x(3)) + abs((p(x)<0)*10000) + (x(2) < 0).*1000000;
% 
%         % include cost functions for when displacement height > minimum height.
%         l = @(x) vel - x(1)./0.41.*log((z-x(2))./x(3));% + abs((p(x)<0)*10000) + (x(2) < 0).*1000000;
%         
%         % fxn to use if you want to adjust d
%         lawall = @(x) l(x);
%         % fxn to use if you want to set d. In this case, z should be
%         % adjusted so that z = z' - d, where z' is the original z.
%         lawall2 = @(x) vel - x(1)./0.41.*log((z)./x(2));
%     
%         % the guesses for the respective fxns above
%         x0 = [ustarguess, dguess, z0guess];      
%         x0_ = [ustarguess, z0guess];
%         % LMF solve set up for the 2 versions of the fit
%         %         [c, resnorm, niter] = LMFsolve(lawall, x0, 'display', 0);
%         [c2, resnorm, niter2] = LMFsolve(lawall2, x0_, 'display', 0);
%         % save coefficients 
%         % c2 = [c2(1); c2(2); c2(3)]; % use this for lawall
%         c2 = [c2(1); 0; c2(2)]; % use this for lawall2
%         c = c2;        
        % **********************LMFSOLVE ENDS HERE*************************
        
        % save values
        ADavg.ustar(i) = c(1);
        ADavg.d(i) = c(2);
        ADavg.z0(i) = c(3);
        ADavg.vel_fit(inonan&iz,i) = lawofwall(c,z);
    else
        % if there wasn't enough data in the profile, make it all NaN.
        ADavg.ustar(i) = NaN;
        ADavg.d(i) = NaN;
        ADavg.z0(i) = NaN;
        ADavg.vel_fit(:,i) = NaN(size(ADavg.vel_fit(:,i)));
    end
end

% running average for ustar; window is 2m+1 (m = 2nd input), so ~1 hour.
% This can be adjusted based on your sampling interval.
ADavg.ustar_runmean = runmean(ADavg.ustar(:),3,1,'edge')';

% Some bonus calculations. Not sure how accurate or useful nu and Re are.
% Calcualte drag coefficient (at 1m), viscocity, and reynolds number

% BW commented out becasue missing ADavg.direction
ADavg.U0 = ADavg.uv(i1m,:);
% % ADavg.direction1m = ADavg.direction(i1m,:);
% % 
% % ADavg.Cd = NaN(size(ADavg.ustar));
% % ADavg.nu = ADavg.Cd;
% % ADavg.Re = NaN(length(ADavg.bin_depth), length(ADavg.Cd)); 
% % ADavg.Cd = (ADavg.ustar.^2)./ADavg.uv(i1m,:).^2; %use U0 from 1 m
% % ADavg.nu = SW_Viscosity(ADavg.TC, 'C', ones(size(ADavg.TC)).*Sal_est, 'ppt')./(gsw_rho(ones(size(ADavg.TC)).*Sal_est, ADavg.TC, ADavg.Pres));
% % ADavg.Re = (ADavg.bin_depth*ADavg.ustar(:)')./(ones(length(ADavg.bin_depth),1)*ADavg.nu);


return


