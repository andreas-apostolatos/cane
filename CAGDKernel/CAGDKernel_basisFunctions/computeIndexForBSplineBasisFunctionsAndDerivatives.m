function index = computeIndexForBSplineBasisFunctionsAndDerivatives ...
    (mixedDerivOrder, xiDerivOrder, etaDerivOrder)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the index one has to use in the columns of the output to function
% computeBSplineBasisFunctionsAndDerivativesForSurface so that to get the
% desirable derivative of the B-Spline basis functions.
%
%           Input :
% mixedDerivOrder : Maximum mixed derivative order of the basis functions
%    xiDerivOrder : Requested derivative in xi-direction
%   etaDerivOrder : Requested derivative in eta-direction
%
%          Output :
%           index : The index to be plugged in into the columns of dN so
%                   that one gets the xiDerivOrder-th derivative with
%                   respect to xi-direction and etaDerivOrder-th derivative
%                   with respect to eta-direction namely dN(:,index)
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the derivative index
%
%% Function main body

%% 0. Read input
if xiDerivOrder+etaDerivOrder>mixedDerivOrder
    error('Requested mixed derivative order exceeds the maximum mixed derivative order of the basis functions');
end

%% 1. Compute the derivative index
index = etaDerivOrder*(2*(mixedDerivOrder+1) - etaDerivOrder + 1)/2 + xiDerivOrder+1;

end
