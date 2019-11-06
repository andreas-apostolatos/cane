function dR = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,mixedDerivOrder)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the IGA basis functions (B-Spline or NURBS) and their derivatives 
% for the parametrization of a surface.
%
% The basis functions array dR(Rxi,Reta,derivs) is sorted as:
%
% dR(:,:) = [dR(1,1) dR(2,1) ... dR(xiNoBasisFct,1) dR(1,2) dR(2,2) ...
%            dR(xiNoBasisFct,2) ... dR(1,etaNoBasisFct) dR(2,etaNoBasisFct) 
%            ... dR(xiNoBasisFct,etaNoBasisFct)]'
%
% where dR denotes the derivatives. The derivatives themselves are sorted
% into the second dimension of the array which given the maximum order of 
% the derivatives n = mixedDerivOrder is given by:
%
% dR(derivs) = [R dR/dxi d^2R/dxi^2 ... d^nR/dxi^n 
%               dR/deta d^2R/deta*dxi ... d^(n-1)R/deta*dxi^(n-1) 
%               dR^2/deta^2 d^3R/deta^2*dxi ... d^(n-2)R/deta^2*dxi^(n-2)
%               ...                         ... ...
%               dR^(n-1)/deta^(n-1) d^(n-1)R/deta^(n-1)*dxi 
%               dR^n/deta^n]
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
% namely in total dR(basisFncIndex,derivIndex)
%
%           Input :
%          xiSpan : The knot span index in xi-direction
%         etaSpan : The knot span index in eta-direction
%               p : The polynomial degree in xi-direction
%               q : The polynomial degree in eta-direction
%          Xi,Eta : The knot vectors in xi-,eta- direction
%              xi : The surface parameter in xi-direction
%             eta : The surface parameter in eta-direction
%              CP : The set of the control point coordinates and weights
%         isNURBS : Flag on wheter the basis is a NURBS or a B-Spline
% mixedDerivOrder : The maximum order for the mixed derivative
%
%          Output :
%              dR : The array containing the B-Spline/NURBS basis functions 
%                   and their derivatives sorted as explained above
%
% % Function layout :
%
% 0. Read input
%
% 1. Compute the B-Spline basis functions and their derivatives
%
% 2. Compute the denominator function and its derivatives
%
% 3. Compute the NURBS basis functions and their derivarives iteratively
%
%% Function main body

%% 0. Read input

% Number of basis functions
noBasisFnc = (p+1)*(q+1);

% Number of derivatives
noDerivs = (mixedDerivOrder+1)*(mixedDerivOrder+2)/2;

% Initialize output array
dR = zeros(noBasisFnc,noDerivs);

%% 1. Compute the B-Spline basis functions and their derivatives
dN = computeBSplineBasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,mixedDerivOrder);
if ~isNURBS
    dR = dN;
    return;
end

%% 2. Compute the denominator function and its derivatives
dF = computeDenominatorFunctionAndDerivativesForSurface(dN,xiSpan,p,etaSpan,q,CP,mixedDerivOrder);

%% 3. Compute the NURBS basis functions and their derivarives iteratively

% Initialize the counter of the basis functions
counterBasis = 1;

% Loop over all the basis functions in eta-direction
for etaBasis=0:q
    % Loop over all the basis functions in xi-direction
    for xiBasis=0:p
        % Loop over all the derivatives in eta-direction
        for l=0:mixedDerivOrder
            % Loop over all the derivatives in xi-direction
            for k=0:mixedDerivOrder-l
                % Get the Control Point indices
                xiIndex = xiSpan - p + xiBasis;
                etaIndex = etaSpan - q + etaBasis;
                
                % Get the index for the B-Spline basis functions
                % derivatives
                indexBSplineDrvs = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder,k,l);
                
                % Store the temporary value
                eta = dN(counterBasis,indexBSplineDrvs)*CP(xiIndex,etaIndex,4);
                
                for j=1:k
                    % Compute the index of the derivatives to the B-Spline
                    % basis functions
                    indexDeriv = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder,k-j,l);
                    
                    % Update the termporary value
                    indexDrvDF = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder,j,0);
                    
                    % Update the parametric coordinate eta
                    eta = eta - getBinomialCoefficients(k,j)*dF(indexDrvDF,1)*dR(counterBasis,indexDeriv);
                end
                
                for i=1:l
                    % Compute the index of the derivatives to the B-Spline
                    % basis functions
                    indexDeriv = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder,k,l-i);
                    
                    % Update the termporary value
                    indexDrvDF = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder,0,i);
                    
                    % Update the parametric coordinate eta
                    
                    eta = eta - getBinomialCoefficients(l,i)*dF(indexDrvDF,1)*dR(counterBasis,indexDeriv);
                    
                    % Initialize second temporary value
                    v2 = 0;
                    
                    for j=1:k
                        % Compute the index of the derivatives to the B-Spline
                        % basis functions
                        indexDeriv = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder,k-j,l-i);

                        % Update the second termporary value
                        indexDrvDF = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder,j,i);
                        
                        % Update the scalar v2
                        v2 = v2 + getBinomialCoefficients(k,j)*dF(indexDrvDF,1)*dR(counterBasis,indexDeriv);
                    end
                    
                    % Update the first temporary variable
                    eta = eta - getBinomialCoefficients(l,i)*v2;
                end
                
                % Update the value of the basis function derivatives
                dR(counterBasis,indexBSplineDrvs) = eta/dF(1,1);
            end
        end
        
        % Update the counter of the basis
        counterBasis = counterBasis + 1;
    end
end

end
