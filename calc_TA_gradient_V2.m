function [SP] = calc_TA_gradient_V2(SP, TA0, Q, npumps)

% [SP] = calc_TA_gradient(SP, TA0, Q)

% Calculates TA graident for pHOxyGF system. 

dTA = NaN(length(Q),length(SP.SDN));
TA = NaN(length(Q), npumps, length(SP.SDN));




for j = 1:length(Q)
    for ii = 1:length(SP.SDN)

        [~,dTA_] = barnes_pHOxy_model(SP.pH(end:-1:1,ii), ...
        SP.DOXY(end:-1:1,ii), SP.TC(end:-1:1,ii), nanmean(SP.PSAL(end:-1:1,ii)), TA0, Q(j), 1);
        % uncomment if you want to use constant salinity of 34.
%         [~,dTA_] = barnes_pHOxy_model(SP.pH(end:-1:1,ii), ...
%         SP.DOXY(end:-1:1,ii), SP.TC(end:-1:1,ii), [34;34;34], TA0, Q(j), 1);
        
        if(npumps == 3)
            TA(j,3,ii) = TA0;
            TA(j,2,ii) = TA0 + dTA_(2);
            TA(j,1,ii) = TA(j,2,ii) + dTA_(3);
        elseif(npumps==2)
            TA(j,2,ii) = TA0;
            TA(j,1,ii) = TA0 + dTA_(2);
        end

    end            
end

SP.TA = TA;



