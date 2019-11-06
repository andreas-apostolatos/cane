function dN = computeBSplineBasisFunctionsAndDerivativesForSurface...
    (xiKnotSpan,p,xi,Xi,etaKnotSpan,q,eta,Eta,mixedDerivOrder)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-Spline basis functions and their derivatives for the 
% parametrization of a surface.
%
% The basis functions array dN(Nxi,Neta,derivs) is sorted as:
%
% dN(:,:) = [dN(1,1) dN(2,1) ... dN(xiNoBasisFct,1) dN(1,2) dN(2,2) ...
%            dN(xiNoBasisFct,2) ... dN(1,etaNoBasisFct) dN(2,etaNoBasisFct) 
%            ... dN(xiNoBasisFct,etaNoBasisFct)]'
%
% where dN denotes the derivatives. The derivatives themselves are sorted
% into the second dimension of the array which given the maximum order of 
% the derivatives n = mixedDerivOrder is given by:
%
% dN(derivs) = [N dN/dxi d^2N/dxi^2 ... d^nN/dxi^n 
%               dN/deta d^2N/deta*dxi ... d^(n-1)N/deta*dxi^(n-1) 
%               dN^2/deta^2 d^3N/deta^2*dxi ... d^(n-2)N/deta^2*dxi^(n-2)
%               ...                         ... ...
%               dN^(n-1)/deta^(n-1) d^(n-1)N/deta^(n-1)*dxi 
%               dN^n/deta^n]
%
% To get the i-th derivative with respect to xi-direction and j-th
% derivative with respect to eta-direction use the index 
%
% derivIndex = (mixedDerivOrder - j)*(mixedDerivOrder - j + 1) + i
%
% and to get the k-th basis function with its derivatives themselves use
% the index
%
% derivIndex = 
%            computeIndexForBSplineBasisFunctionsAndDerivatives(nDeriv,i,j)
%
% namely in total dN(basisFncIndex,derivIndex)
%
%           Input :
%      xiKnotSpan : The knot span index in xi-direction
%     etaKnotSpan : The knot span index in eta-direction
%               p : The polynomial degree in xi-direction
%               q : The polynomial degree in eta-direction
%          Xi,Eta : The knot vectors in xi-,eta- direction
%              xi : The surface parameter in xi-direction
%             eta : The surface parameter in eta-direction
% mixedDerivOrder : The maximum order for the mixed derivative
%
%          Output :
%              dN : The array containing the B-Spline basis functions and
%                   their derivatives sorted as explained above
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the B-Spline basis functions and their derivatives at each parametric direction xi, eta
%
% 2. Compute the B-Spline basis functions and their derivatives for a surface up to maximum mixed derivative order mixedDerivOrder
%
%% Function main body

%% 0. Read input

% Compute the number of the basis functions affecting the current knot span
noBasisFcts = (p+1)*(q+1);

% Compute the number of the derivatives to be computed
noDrvs = (mixedDerivOrder+1)*(mixedDerivOrder+2)/2;

% Initialize output array
dN = zeros(noBasisFcts,noDrvs);

% Initialize counter for the basis functions
counterBasis = 1;

% Initialize counter for the derivatives of the basis functions
counterDerivs = 1;

% Find the minimum between the derivative order and the polynomial degree :

% in xi-direction
dxi = min(mixedDerivOrder,p);

% in eta-direction
deta = min(mixedDerivOrder,q);

%% 1. Compute the B-Spline basis functions and their derivatives at each parametric direction xi, eta

% Compute the B-Spline basis functions along xi
dNxi = computeBSplineBasisFunctionsAndDerivativesForCurve(xiKnotSpan,p,xi,Xi,mixedDerivOrder);

% Compute the B-Spline basis functions along eta
dNeta = computeBSplineBasisFunctionsAndDerivativesForCurve(etaKnotSpan,q,eta,Eta,mixedDerivOrder);

%% 2. Compute the B-Spline basis functions and their derivatives for a surface up to maximum mixed derivative order mixedDerivOrder
for etaBasis = 1:q+1
    for xiBasis = 1:p+1
        for l = 1:deta+1
%             dd = min(mixedDerivOrder-l+1,dxi);
            dd = mixedDerivOrder-l+1;
            for k=1:dd+1
                % Compute the derivatives of the B-Spline basis functions
                dN(counterBasis,counterDerivs) = dNxi(xiBasis,k)*dNeta(etaBasis,l);
                
                % Update the derivative index
                counterDerivs = counterDerivs + 1;
            end
        end
        
        % Reset the derivative index to one
        counterDerivs = 1;
        
        % Update the basis functions index
        counterBasis = counterBasis + 1;
    end
end

end

