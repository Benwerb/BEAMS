function [I] = get_I_SW(Sal)

% Returns the ionic strength of SW as a function fo salinity.
% From Dickson 2007, Best practices.
% Units are in mol/kg-H2O

I = 19.924.*Sal./(1000-1.005*Sal);   % mol/kg-H2O.


return
