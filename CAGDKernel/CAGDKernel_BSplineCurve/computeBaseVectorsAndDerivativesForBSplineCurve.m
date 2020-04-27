function dG = computeBaseVectorsAndDerivativesForBSplineCurve ...
    (xiSpan, p, CP, mixedDerivOrder, dR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the base vectors and their derivatives up to mixed derivative
% order mixedDerivOrder for a B-Spline curve.
%
%          Input :
%         xiSpan : The knot span indices of the surface parameter where the
%                  base vectors are computed
%              p : The polynomial orders in xi-direction for the B-Spline 
%                  basis
%             CP : The set of Control Point coordinates and weights
%             dR : The set of basis functions and their derivatives
%
%         Output :
%             dG : The array containing the base vector G1 and its
%                  derivatives sorted as:
%                  % dG^0/dxi^0 d^2G/dxi^2 ... d^mixedDerivOrderG/dxi^mixedDerivOrder
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the knot span contributions in xi-direction
% ->
%    1i. Loop over all the derivatives
%    ->
%        1i.1. Get the Control Point indices
%
%        1i.2. Get the Control Point indices
%    <-
%
%   1ii.  Update counter for the basis
% <-
%
%% Function main body

%% 0. Read input

% Check input
if length(dR(1,:)) < mixedDerivOrder
    error('The number of requested derivatives for the base vector is higher than the number of derivatives provided by the basis functions dR');
end

% Initialize counters
counterBasis = 1;

% Initialize output array
noCoordinates = 3;
dG = zeros(noCoordinates,mixedDerivOrder);

%% 1. Loop over all the knot span contributions
for iBasis = 0:p
    %% 1i. Loop over all the derivatives
    for iDeriv = 1:mixedDerivOrder
        %% 1i.1. Get the Control Point indices
        xiIndex = xiSpan - p + iBasis;
        
        %% 1i.2. Get the Control Point indices
        for iCoord = 1:noCoordinates
            dG(iCoord,iDeriv) = dG(iCoord,iDeriv) + dR(counterBasis,iDeriv)*CP(xiIndex,iCoord);
        end
    end

    %% 1ii.  Update counter for the basis
    counterBasis = counterBasis + 1;
end

end
