function [TA_model, delTA, DIC_model, C] = barnes_pHOxy_model(pH, DOXY, TC, Sal, pre_TA, Q, delOXY_fraction)

% [TA_model, delTA, C] = barnes_pHOxy_model(pH, DOXY, TC, Sal, pre_TA, Q, delOXY_fraction)
% 
% Models TA based on pH, O2, Temperature, and Salinity data. Model based on
% Barnes 1983. This function expects pH, DOXY, TC, and Sal to be an array of the same size, and
% pre_TA, Q, and delOXY_fraction to be a scalar. pH is on the total scale at in situ conditions,
% DOXY is O2 [umol/kg], TC is temperature in Celsius, and Sal is practical salinity. Q is the
% assumed metabolic quotient. 
% 
% It assumes that the starting TA is pre_TA, and then calculates changes in TA using pH, DOXY, TC, 
% and Sal. So for example, if you want to calculate the deltaTA for the gradient flux, the inputs
% for pH, DOXY, TC, and Sal should be a 2x1 (or 1x2) array, with measurements from z=1 and z=2. If
% it was an incubation, then pH, DOXY, TC, and Sal should be a 1xn array, where n is the number of
% measurements made during the incubation. 
%
%
% delOXY_fraction is the fraction (0-1) of the deltaO2 that is considered
% to be metabolic. so if delOXY_fraction is 0.8, that means 80% of the
% change in oxygen is assumed to be due to metabolism, and 20% due to gas
% exchange/advection. 
% 
% 
% Yui Takeshita
% MBARI 
% 

% *************CALCULATE EQUILIBRIUM CONSTANTS****************************
TK = TC+273.15; %Convert temperatureTC 
[KC1, KC2] = calcKC_sw(TK, Sal);
%Get KB (Dickson 1990)
KB = calcKB_sw(TK, Sal);
KF = calcKF_sw(TK, Sal);
KW = calcKW_sw(TK, Sal);
KS = calcKS_sw(TK, Sal);
% Total ion concentrations
TB = get_Tion(Sal,'B'); % Total Boron concentration, umol/kg-sol
TF = get_Tion(Sal,'F'); % Total Fluoride concentration, umol/kg-sol
TS = get_Tion(Sal,'S'); % Total sulfate concentration, umol/kg-sol
% ************************************************************************

% =========== Calculate Model Parameters ================
H_total = 10.^-pH; %mol/kg-sol.
H_free  = H_total./(1+TS./KS);
BA = TB.*KB./(KB + H_total).*1e6; % Borate Alkalinity. [B(OH)4-]. 
OH = KW./H_total.*1E6; % Hydroxide alkalinity. also [OH-]. [umol/kg-sol]
FA = (H_free.*TF)./(KF + H_free)*1e6; %HF alkalinity.
K = (H_total.*KC1 + KC1.*KC2 + H_total.^2)./(H_total.*KC1 + 2.*KC1.*KC2);
delO2 = zeros(size(DOXY)); 
delO2(2:end) = diff(DOXY).*delOXY_fraction;
% ========================================================

% Brakceted part in equation 2. It is DIC/CA, where CA = carbonate
% alkalinity.

% Initialize model
% Estimate initial Alk. 
delTA = NaN(size(TC));
TA_model = delTA;
DIC_model = delTA;
% delTA_ = delTA;
a = delTA;
b = delTA;
c = delTA;
d = delTA;
TA_model(1)  = pre_TA; %initial guess at TA. based on pH/pCO2. subject to change.
DIC_model(1) = (pre_TA - BA(1) - OH(1) + FA(1) + H_total(1)*1e6).*K(1);


%ignore first point. 
for it=2:length(TC)
    
    %Reinitialize model (TA) every day, at midnight.
%     if(rem(SDN(it),1) == 0)
%         TA_model(it-1) = TA(it-1);
%     end
%      it_ = it + i_start - 1;
    
    %Modelled change in TA. Based on equation 9, Barnes 1983.
    a(it) = (delO2(it).*Q)./ (K(it)-0.5);
%     a(it) = -(delDIC(it_)./Q_prime) ./ (K(it_)-0.5);
%     b(it) = (K(it)-K(it-1)).*pre_TA./ (K(it)-0.5);
    b(it) = (K(it)-K(it-1)).*TA_model(it-1)./(K(it)-0.5);
    c(it) = - K(it).*(BA(it) + OH(it) - FA(it) - H_total(it)*1e6)./ (K(it)-0.5);
    d(it) = K(it-1).*(BA(it-1)+OH(it-1) - FA(it-1) - H_total(it-1)*1e6)./ (K(it)-0.5);
%     c(it) = - K(it).*(BA(it) + OH(it))./ (K(it)-0.5);
%     d(it) = K(it-1).*(BA(it-1)+OH(it-1))./ (K(it)-0.5);
    
    %Sum the components, and calculate the change in alkalinity.
    delTA(it) = -(a(it) + b(it) + c(it) + d(it));
%     delTA(it) = TA_model(it) - TA_model(1);

    %Add alkalinity to model
    TA_model(it) = TA_model(it-1) + delTA(it);
    % Add DIC from model
    DIC_model(it) = (TA_model(it) - BA(it) - OH(it) + FA(it) + H_total(it)).*K(it);
       
end

C(:,1) = a;
C(:,2) = b;
C(:,3) = c;
C(:,4) = d;

