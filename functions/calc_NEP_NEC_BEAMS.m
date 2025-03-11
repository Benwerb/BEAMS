function BMS = calc_NEP_NEC_BEAMS(BMS, ipumps, ustarvar)

% BMS = calc_NEP_NEC_BEAMS(BMS, ipumps)
%
% BMS is the BEAMS structure, ans ipumps is which pumps to use for
% calculating NEP/NEC using equations from McGillis 2011. Top pump should
% be ipumps(1), and bototm pump shoudl be ipumps(2). Adds to the BMS
% strcuture NEP_WM, NEC_WM, NEP, and NEC categories. NEP and NEC are
% benthic metabolism rates calculated using a fit.
%
% Author: Yui Takeshita
% Created 11/16/2015
% Last modified: 11/16/2015


% ***********Calculate metabolism based on McGillis 2011***************
% NEP calculations
karman = 0.41;
% dt = diff(BMS.SDN); dt = [NaN; dt]'; % time steps in days
deltaO2 = BMS.DOXY(ipumps(1),:) - BMS.DOXY(ipumps(2),:); %gradient from top and bottom pump
BMS.NEP_WM = -BMS.DENS(1,:).*BMS.(ustarvar).*karman.*deltaO2./(log(BMS.pumpz(ipumps(1))/BMS.pumpz(ipumps(2)))).*3600./1000;
% NEC calculations. loop through the 5 different Q scenarios. 
BMS.NEC_WM = NaN(size(BMS.TA,1), length(BMS.SDN));
for ii = 1:size(BMS.TA,1)
    deltaTA = vec2row(squeeze(BMS.TA(ii,ipumps(1),:) - BMS.TA(ii,ipumps(2),:))); % gradient from top and bottom pump
    BMS.NEC_WM(ii,:) = BMS.DENS(1,:).*BMS.(ustarvar).*karman.*deltaTA./(log(BMS.pumpz(ipumps(1))/BMS.pumpz(ipumps(2)))).*3600./1000./2;
end
% *********************************************************************

% Presize NEP and NEC
BMS.NEP = NaN(size(BMS.SDN));
BMS.NEC = NaN(5,length(BMS.SDN));

% ========== Calculate metabolism based on fitting method =============
for ii = 1:length(BMS.SDN)
    % equation is the form: c1*exp(c2(z))+c3
    
    % ================= fitting parameters ====================
    c1_0 = 0;  % flux
    % fit 
    c2_0_DOXY = BMS.DOXY(1,ii);       % initial guess DO
    if(isnan(c2_0_DOXY)); c2_0_DOXY = nanmean(BMS.DOXY(:,ii)); end;
    c_DOXY = [c1_0 c2_0_DOXY];
    lb_DOXY = [c1_0-0.05 c2_0_DOXY-10];
    ub_DOXY = [c1_0+0.05 c2_0_DOXY+10];
    % ==========================================================

    iuseO2 = ~isnan(BMS.DOXY(:,ii));
    DOXYz = BMS.DOXY(iuseO2,ii);

    if(sum(~isnan(BMS.DOXY(:,ii))) <= 1 || isnan(BMS.U0(ii)))
        BMS.NEP(ii) = NaN;
        BMS.NEP_R2(ii) = NaN;
        BMS.NEC(:,ii) = NaN;
        BMS.NEC_R2(:,ii) = NaN;
    else
        options = optimset('Display','none', 'diagnostics', 'off');
        [c_fit_DOXY, resnorm, ~, flags(1,ii)]  = lsqcurvefit(@conc_prof_BL, c_DOXY, ...
            BMS.pumpz(iuseO2), vec2row(DOXYz), lb_DOXY, ub_DOXY, options);
        BMS.NEP(ii) = c_fit_DOXY(1)*3600;
        c2_fit(1,ii) = c_fit_DOXY(2);
        BMS.NEP_R2(ii) = 1 - (resnorm./(sum((DOXYz-nanmean(DOXYz)).^2)));

        % if statment to see if there are TA values
        if(sum(~isnan(BMS.TA(1,:,ii))) > 1)
            for j = 1:size(BMS.TA,1)

                % =========== Fitting parameter for TA =============
                c2_0_TA = BMS.TA(j,1,ii);
                if(isnan(c2_0_TA)); c2_0_TA = nanmean(BMS.TA(j,:,ii)); end;
                c_TA = [c1_0 c2_0_TA];
                lb_TA = [c1_0-0.05 c2_0_TA-10];
                ub_TA = [c1_0+0.05 c2_0_TA+10];    
                % ==================================================

                iuseTA = ~isnan(BMS.TA(j,:,ii));
%                 iuseTA = [1 3]; 
%                 TAz = [2205+trex.dTA(j,ii); 2205];
                TAz = BMS.TA(j,iuseTA,ii)';
                [c_fit_TA, resnorm, ~, flags(1,ii)]  = lsqcurvefit(@conc_prof_BL, c_TA, ...
                    BMS.pumpz(iuseTA), vec2row(TAz), lb_TA, ub_TA, options);
                BMS.NEC(j,ii) = -c_fit_TA(1)*3600./2;
                c2_fit(2,ii) = c_fit_TA(2);
                BMS.NEC_R2(j,ii) = 1 - (resnorm./(sum((TAz-nanmean(TAz)).^2)));
            end
        end
%         swelldates = [datenum('9/12/2014') datenum('9/13/2014')];
%         NCP(:,inrange(SDN, swelldates)) = NaN;
%         c2_fit(inrange(SDN, swelldates)) = NaN;
%         NCC(:,inrange(SDN, swelldates)) = NaN;
    end
end

    % nested function used to fit chemical gradient
    function C = conc_prof_BL(c,z)

    % Nested function for law of the wall.
    % a(1) = ustar (friction velocity)
    % a(2) = d (displacement height)
    % a(3) = 
    
    % z is the height of the bins.
    d = 0.1; % set displacement height to 0.1 m. 
    C = -(c(1)/(BMS.(ustarvar)(ii)*0.41))*log(z/d)+c(2);
    end

end



