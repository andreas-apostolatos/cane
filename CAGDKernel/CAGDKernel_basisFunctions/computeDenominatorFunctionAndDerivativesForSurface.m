function dF = computeDenominatorFunctionAndDerivativesForSurface(dN,xiSpan,p,etaSpan,q,CP,nDeriv)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the denominator function dF = Sum_i N_i(u,v) * w_i at the
% parametric location where the BSpline basis functions and their
% derivatives have been evaluated. Source reference:
%
% Les Piegl and Wayne Tiller, The NURBS Book. Springer-Verlag, Berlin 1995
% p. 72.
%
%          Input :
%             dN : The B-Spline basis functions and their derivatives
% xiSpan,etaSpan : The knot span indices in the u-,v- directions
%            p,q : The polynomial degees of the basis in u-,v- directions
%             CP : The set of the control point coordinates and weights
%         nDeriv : The number of derivatives to be computed
%
%         Output :
%             dF : The denominator function dF = Sum_i N_i(u,v) * w_i
%
% Function main body
%
% 0. Read input
%
% 1. Compute the denominator function iteratively
%
%% Function main body

%% 0. Read input

% Initialize output
dF = zeros((nDeriv+1)*(nDeriv+2)/2,1);

% Counter for the derivatives
counterDrvs = 1;

% Counter for the basis functions
counterBasis = 1;

%% 1. Compute the denominator function iteratively

% Loop over all partial derivatives with respect to v-direction
for j=0:nDeriv
    % Loop over all the partial derivatives with respect to u-direction
    for i=0:nDeriv-j
        % Loop over all the contributions from each basis function in v-direction
        for l=0:q
            % Loop over all the contributions from each basis function in u-direction
            for k=0:p
                % Get the Control Point indices
                xiIndex = xiSpan - p + k;
                etaIndex = etaSpan - q + l;
                
                % Get the derivative index
                indexDrv = computeIndexForBSplineBasisFunctionsAndDerivatives(nDeriv,i,j);
                
                % Sum up the contribution from each basis function
                dF(counterDrvs,1) = dF(counterDrvs,1) + dN(counterBasis,indexDrv)*CP(xiIndex,etaIndex,4);
                
                % Update the counter for the basis functions
                counterBasis = counterBasis + 1;
            end
        end
        % Reset the counter for the basis
        counterBasis = 1;
        
        % Update the counter for the derivatives
        counterDrvs = counterDrvs + 1;
    end
end

end
