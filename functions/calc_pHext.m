% *************************************************************************
% calc_pH_ext
%
% Function to calcualte pH_ext. For details, see Martz et al. 2010. Code
% based on Bresnahan et al. 2014.
% *************************************************************************

function [pHext_tot, pHext_free] = calc_pHext(Vext, tempC, salt, k0, k2)

% [pHext_tot, pHext_free] = calc_dfet_pHext(Eext, tempC, salt, Eoext25)
%
% Calculates pH using Honeywell Durafet and external reference electrode
% (Chloride Ion Selective Electrode) in seawater
%
% Code from Bresnahan et al. 2014, Methods in Oceanography
%

% Univ gas constant, Faraday constant,
R = 8.3145; F = 96487;
% See Martz et al. 2010 for greater detail
tempK = tempC+273.15; % Convert temp from C to K
S_T = (R*tempK)/F*log(10); % Nernst temp dependence
Z = 19.924.*salt./(1000-1.005.*salt); % Ionic strength, Dickson et al. 2007
SO4_tot = (0.14/96.062).*(salt./1.80655); % Total conservative sulfate
cCl = 0.99889./35.453.*salt./1.80655; % Conservative chloride
mCl = cCl.*1000./(1000-salt.*35.165/35); % mol/kg-H2O
K_HSO4 = exp(-4276.1./tempK+141.328-23.093.*log(tempK)...
    +(-13856./tempK+324.57-47.986.*log(tempK)).*Z.^0.5...
    +(35474./tempK-771.54+114.723.*log(tempK)).*Z-2698./tempK.*Z.^1.5...
    +1776./tempK.*Z.^2+log(1-0.001005.*salt)); % Bisulfate equilibrium const., Dickson et al. 2007

DHconst = 0.00000343.*tempC.^2+0.00067524.*tempC+0.49172143; % Debye-Huckel, Khoo et al. 1977
log10gamma_HCl = 2*(-DHconst.*sqrt(Z)./(1+1.394*sqrt(Z))+(0.08885-0.000111*tempC).*Z);
pHext_free = -(((k0+k2*(tempC))-Vext)-S_T.*(log10(mCl)+log10gamma_HCl))./S_T; % mol/kg-H2O
pHext_free = pHext_free-log10((1000-salt.*35.165/35)/1000); % mol/kg-sw
pHext_tot = pHext_free-log10(1+SO4_tot./K_HSO4); % mol/kg-sw

end