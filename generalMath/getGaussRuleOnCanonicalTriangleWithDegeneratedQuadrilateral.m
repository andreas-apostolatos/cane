function [GP, GW] = getGaussRuleOnCanonicalTriangleWithDegeneratedQuadrilateral ...
    (numGP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
%% Function documentation
%
% Returns the Gauss point coordinates and weights for the integration of a
% function over the canonical triangle using the concept of the degenerated
% quadrilateral. The rule is constructed symmetrically, that is, the
% tensor product rule is chosen homogeneously for both parametric 
% directions of the integration space.
%
%       Input :
%       numGP : Number of quadrature points
%
%      Output :
%          GP : array containing the Gauss point coordinates in a list
%               GP = [zeta1_1 zeta2_1
%                     ...     ...
%                     zeta1_n zeta2_n]
%          GW : array containing the Gauss point weights in a vector
%
% Function layout :
%
% 0. Read input
%
% 1. Initialize output arrays
%
% 2. Issue the Gauss point coordinates and weights over the canonical 1d space
%
% 3. Compute the Gauss point coordinates and weights over the canonical quadrilateral
%
% 4. Transform the Gauss points into the canonical triangle by degenerating the quadrilateral
%
%% Function main body

%% 0. Read input

% Number of Gauss points in one dimension
noGP1D = sqrt(numGP);

% Check input
if isinteger(noGP1D)
    error('Inappropriate number of Gauss points requested');
end

% Initialize Gauss point counter
counterGP = 1;

%% 1. Initialize output arrays
GP = zeros(numGP,2);
GW = zeros(numGP,1);

%% 2. Issue the Gauss point coordinates and weights over the canonical 1d space
[GP1D,GW1D] = getGaussPointsAndWeightsOverUnitDomain(noGP1D);

%% 3. Compute the Gauss point coordinates and weights over the canonical quadrilateral
for iEtaGP = 1:noGP1D
    for iXiGP = 1:noGP1D
        GP(counterGP,1) = GP1D(1,iXiGP);
        GP(counterGP,2) = GP1D(1,iEtaGP);
        GW(counterGP,1) = GW1D(1,iXiGP)*GW1D(1,iEtaGP);
        counterGP = counterGP + 1;
    end
end

%% 4. Transform the Gauss points into the canonical triangle by degenerating the quadrilateral
for iGP = 1:numGP
    GP(iGP,1) = (1 + GP(iGP,1))*(1 - GP(iGP,2))/4;
    GP(iGP,2) = (1 + GP(iGP,2))/2;
    detJ = abs((1 - GP(iGP,2))/4);
    GW(iGP,1) = GW(iGP,1)*detJ;
end

end
