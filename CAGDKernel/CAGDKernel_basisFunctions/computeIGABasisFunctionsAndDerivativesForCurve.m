function dR = computeIGABasisFunctionsAndDerivativesForCurve...
    (knotSpanIndex, p, xi, Xi, CP, isNURBS, numDeriv)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns an array containing the B-Spline/NURBS basis functions and their
% derivatives up to order numDeriv at the chosen parametric location.
%
%         Input :
% knotSpanIndex : Knot span index
%             p : The polynomial degree of the B-Spline curve
%            xi : The curve parameter where the basis functions are to be
%                 evaluated
%            Xi : The knot vector of the B-Spline curve
%            CP : The set of Control Point coordinates and weights of the 
%                 B-Spline curve
%       isNURBS : Flag on the whether the basis is a B-Spline or a NURBS
%      numDeriv : The number of derivatives to be computed
%
%        Output :
%            dR : Array containing the B-Spline/NURBS basis functions and 
%                 their derivatives up to order numDeriv. Element dR(i,j), 
%                 i=1,...,p+1 and j=1,...,numDeriv+1 returns the (j+1)-th 
%                 derivative of the i-th non identically zero basis 
%                 function at the knot span with index knotSpanIndex
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the BSpline basis functions and their derivatives at xi
%
% 2. Compute the denominator function and its derivatives at xi
%
% 3. Loop over all the basis functions
%
%    3i. Loop over all the derivatives
%
%        3i.2. Compute the Control point index
%
%        3i.3. Compute the product of the derivatives of the basis functions with the Control Point weights
%
%        3i.4. Loop over all the involved derivatives
%
%        3i.5. Divide by the denominator function
%
%% Function main body

%% 0. Read input

% Initialize the output array
dR = zeros(p+1,numDeriv+1);

%% 1. Compute the BSpline basis functions and their derivatives at xi
dN = computeBSplineBasisFunctionsAndDerivativesForCurve(knotSpanIndex,p,xi,Xi,numDeriv);
if ~isNURBS
    dR = dN;
    return;
end

%% 2. Compute the denominator function and its derivatives at xi
dF = computeDenominatorFunctionAndDerivativesForCurve(dN,knotSpanIndex,p,CP,numDeriv);

%% 3. Loop over all the basis functions
for i = 0:p
    %% 3i. Loop over all the derivatives
    for j = 0:numDeriv
        %% 3i.2. Compute the Control point index
        index = knotSpanIndex-p+i;
        
        %% 3i.3. Compute the product of the derivatives of the basis functions with the Control Point weights
        v = dN(i+1,j+1)*CP(index,4);
        
        %% 3i.4. Loop over all the involved derivatives
        for k=1:j
            v = v - nchoosek(j,k)*dF(k+1,1)*dR(i+1,j-k+1);
        end
        
        %% 3i.5. Divide by the denominator function
        dR(i+1,j+1) = v/dF(1,1);
    end
end

end

