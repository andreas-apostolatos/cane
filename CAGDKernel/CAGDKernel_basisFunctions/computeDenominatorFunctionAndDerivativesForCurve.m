function dF = computeDenominatorFunctionAndDerivativesForCurve(dN,knotSpanIndex,p,CP,nDeriv)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the denominator functon for the NURBS basis functions namely
% Sum_i (N_i * CP_i) where N_i is the set of the B-Spline basis functions
% as well as its first n derivatives needed for the computation of the 
% NURBS basis functions parametrizing a curve in 2D/3D. Source reference:
%
% Les Piegl and Wayne Tiller, The NURBS Book. Springer-Verlag: Berlin 1995
% p. 72.
%
%         Input :
%            dN : The set of the non-zero BSpline basis functions and their
%                 (n+1)-first derivatives
% knotSpanIndex : Knot span index
%             p : The polynomial degree of the curve
%            CP : The set of Control Point coordinates and weights
%        nDeriv : The number of derivatives to be computed for the
%                 denominator function
%   
%        Output : 
%            dF : The denominator function dF = Sum_i (N_i * CP_i) sorted 
%                 in an array dF(j) where j=1,...,n+1 is the (j-1)-th 
%                 derivative of the denominator function
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the derivatives
%
%    1i. Loop over all the basis functions
%
%        1i.1. Get the control point index
%
%        1i.2. Compute the denominator function iteratively
%
%% Function main body

%% 0. Read input

% Initialize output array
dF = zeros(nDeriv+1,1);

%% 1. Loop over all the derivatives
for i=1:nDeriv+1
    %% 1i. Loop over all the basis functions
    for j=0:p
        %% 1i.1. Get the Control Point index
        index = knotSpanIndex-p+j;
        
        %% 1i.2. Compute the denominator function iteratively
        dF(i,1) = dF(i,1) + dN(j+1,i)*CP(index,4);
    end
end

end

