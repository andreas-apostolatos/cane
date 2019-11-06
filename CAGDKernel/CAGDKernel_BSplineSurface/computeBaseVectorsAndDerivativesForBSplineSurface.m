%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dGXi,dGEta] = computeBaseVectorsAndDerivativesForBSplineSurface...
    (xiSpan,p,etaSpan,q,CP,mixedDerivOrder,dR)
%% Function documentation
%
% Returns the base vectors and their derivatives up to mixed derivative
% order mixedDerivOrder. Since GXi = dX/dxi and GEta = dX/deta, there is a
% symmetry in higher derivatives expressed as:
%
% d^(l+1 + m+1)X/dxi^(l+1)/deta^(m+1) = d^(l + m+1)GXi/dxi^l/deta^(m+1)
% d^(l+1 + m+1)X/dxi^(l+1)/deta^(m+1) = d^(l+1 + m)GEta/dxi^(l+1)/deta^m
% ==> d^(l+1 + m)GEta/dxi^(l+1)/deta^m = d^(l + m+1)GXi/dxi^l/deta^(m+1)
%
% Therefore all the higher order derivatives can be stored into the array 
% dGXi which contains all the derivatives of base vector GXi.
%
%          Input :
% xiSpan,etaSpan : The knot span indices of the surface parameter where the
%                  base vectors are computed
%            p,q : The polynomial orders in xi-,eta-direction for the
%                  B-Spline basis
%             CP : The set of Control Point coordinates and weights
%             dR : The set of basis functions and their derivatives
%
%         Output :
%           dGXi : The array containing the base vector G1 and its
%                  derivatives sorted as:
%                  dGXi = [GXi                   dGXi/dxi        d^2GXi/d^2xi      ... d^(n-1)GXi/d^(n-1)xi  d^nGXi/d^nxi
%                          dGXi/deta             d^2GXi/deta/dxi d^3Gxi/deta/dxi^2 ... d^nGXi/deta/d^(n-1)xi
%                          ...                   ...             ...               ...
%                          dGXi^(n-1)/deta^(n-1) dGXi^n/deta^(n-1)/dxi 
%                          dGXi^n/deta^n]
%          dGEta : The array containing the base vector G2 and its 
%                  derivatives sorted as:
%                  dGEta = [GEta dGEta/deta d^2GEta/d^2eta ... d^(n-1)GEta/d^(n-1)eta  d^nGEta/d^neta]
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the entries of the base vectors and their derivatives iteratively
%
%% Function main body

%% 0. Read input

% Initialize counters
counterBasis = 1;
counterBaseVct1 = 1;
counterBaseVct2 = 1;

% Initialize output array
noCoordinates = 3;
nXiBaseVctAndDrvs = (mixedDerivOrder+1)*(mixedDerivOrder+2)/2;
nEtaBaseVctAndDrvs = mixedDerivOrder+1;
dGXi = zeros(noCoordinates,nXiBaseVctAndDrvs);
dGEta = zeros(noCoordinates,nEtaBaseVctAndDrvs);

% Check input
% nDrvsBasis = length(dR(1,:));
% if nDrvsBasis < nBaseVctAndDrvs-1
%     error('Not enough basis function derivatives were provided for the computation of the derivatives for the base vectors');
% end

%% 1. Compute the entries of the base vectors and their derivatives iteratively

% Loop over all the knot span contributions in eta-direction
for etaBasis = 0:q
    % Loop over all the knot span contributions in xi-direction
    for xiBasis = 0:p
        % Loop over all the derivatives in eta-direction
        for l = 0:mixedDerivOrder
            % Loop over all the derivatives in xi-direction
            for k = 0:mixedDerivOrder-l
                % Get the Control Point indices
                xiIndex = xiSpan - p + xiBasis;
                etaIndex = etaSpan - q + etaBasis;

                % Case l == 0 :: derivative with respect to eta is zero
                % thus all resulting vectors are the k-th derivatives 
                % of base vector G1 w.r.t. xi

                % Case k == 0 :: derivative with respect to xi is zero
                % thus all resulting vectors are the l-th derivatives 
                % of base vector G2 w.r.t. eta

                % Case l > 0 || k > 0 :: This will be the l-th w.r.t xi
                % and (k-1)-th w.r.t. eta derivative of G1 base vector
                % which equals the (l-1)-th w.r.t. xi and k-th w.r.t.
                % eta derivative of G2 base vector
                
                indexBasisFctDrv = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder+1,k+1,l);
                for coord = 1:3
                    dGXi(coord,counterBaseVct1) = dGXi(coord,counterBaseVct1) + dR(counterBasis,indexBasisFctDrv)*CP(xiIndex,etaIndex,coord);
                end
                counterBaseVct1 = counterBaseVct1 + 1;
                    
                if k == 0
                    indexBasisFctDrv = computeIndexForBSplineBasisFunctionsAndDerivatives(mixedDerivOrder+1,k,l+1);
                    for coord = 1:3
                        dGEta(coord,counterBaseVct2) = dGEta(coord,counterBaseVct2) + dR(counterBasis,indexBasisFctDrv)*CP(xiIndex,etaIndex,coord);
                    end
                    counterBaseVct2 = counterBaseVct2 + 1;
                end
            end
        end
        % Reset the index of the derivatives to the basis functions as well
        % as the counters for both base vectors
        counterBaseVct1 = 1;
        counterBaseVct2 = 1;
        
        % Update counter for the basis
        counterBasis = counterBasis + 1;
    end
end

end

