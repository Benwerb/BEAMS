function F = lawofwall(a,z)

% Function for law of the wall.
% a(1) = ustar (friction velocity)
% a(2) = d (displacement height)
% a(3) = z0 (roughness scale)

% z is the height of the bins.


F = a(1)./0.41.*log((z-a(2))./a(3));

