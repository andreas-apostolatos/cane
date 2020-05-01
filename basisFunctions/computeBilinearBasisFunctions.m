function [N, isInside] = computeBilinearBasisFunctions(u, v)
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
%        N : A 4x1 array containing the values of the basis functions at 
%            the given parametric locations
% isInside : If the given parametric coordinates are inside the canonical
%            quadrilateral
%
% Function layout :
%
% 1. Compute the bilinear basis functions
%
% 2. Check if the given point is inside the canonical quadrilateral
%
%% Function main body

% Initialize flag
isInside = true;

%% 1. Compute the bilinear basis functions
N = [(1 + u)*(1 + v)/4
	 (1 + u)*(1 - v)/4
     (1 - u)*(1 - v)/4
     (1 - u)*(1 + v)/4];

%% 2. Check if the given point is inside the canonical quadrilateral
if N(1,1) < -1 || N(1,1) > 1 || N(2,1) < -1 || N(2,1) > 1 || ...
        N(3,1) < -1 || N(3,1) > 1 || N(4,1) < -1 || N(4,1) > 1
    isInside = false;
end

end
