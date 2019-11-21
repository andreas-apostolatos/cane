function [dN,isInside] = computeBilinearBasisFunctionsAndFirstDerivatives(u,v)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the bilinear shape functions corresponding to a quadrilateral
% element in a counterclock-wise fashion:
%
%                      eta
%                      / \
%         N4-(-1,+1)    |    N1-(+1,+1)
%              _________|_________
%              |        |        |
%              |        |        |
%              |        |        |
%        -------------------------------> xi
%              |        |        |
%              |        |        |
%              |        |        |
%              _________|_________
%         N3-(-1,-1)    |    N2-(+1,-1)
%                       |
%
%    Input :
%      u,v : The parameters on the parametric domain
%
%   Output :
%       dN : 4x3 array for which dN(:,1) contains the basis functions
%            themselves, dN(:,2) the derivatives of the basis functions 
%            with respect to u-direction and dN(:,3) the derivatives of the 
%            basis functions with respect to v-direction.
% isInside : Flag on whether the point where the basis functions are to be
%            computed is inside or outside the quadrilateral
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the bilinear basis functions and their parametric derivatives
%
% 2. Check if the given point is inside the canonical quadrilateral
%
%% Function main body

%% 0. Read input

% Initialize the output array
dN = zeros(4,3);

% Initialize output array
isInside = true;

%% 1. Compute the bilinear basis functions and their parametric derivatives

% Assign the values of the functions at the parametric locations
dN(:,1) = [(1+u)*(1+v)/4
           (1+u)*(1-v)/4
           (1-u)*(1-v)/4;
           (1-u)*(1+v)/4];

% Derivatives dN/du
dN(:,2) = [(1+v)/4
           (1-v)/4
           -(1-v)/4
           -(1+v)/4];

% Derivatives dN/dv
dN(:,3) = [(1+u)/4
           -(1+u)/4
           -(1-u)/4
           (1-u)/4];
       
%% 2. Check if the given point is inside the canonical quadrilateral
if dN(1,1) < -1 || dN(1,1) > 1 || dN(2,1) < -1 || dN(2,1) > 1 || ...
        dN(3,1) < -1 || dN(3,1) > 1 || dN(4,1) < -1 || dN(4,1) > 1
    isInside = false;
end

end
