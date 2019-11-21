function [zeta,N,isInside] = computePointProjectionAndBasisFunctionsOnLinearTriangle(P,P1,P2,P3)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the canonical coordinates (zeta1, zeta2) in the canonical space
% of a linear triangle when projecting a point onto it as well as the basis
% functions of the linear triangle in the canonical space and a flag
% indicating whether the projection was found inside the triangle.
%
%       Input :
%           P : The point to be projected
%    P1,P2,P3 : The vertices of the triangle
%
%      Output : 
%        zeta : = [zeta1, zeta2] the parametric coordinates of the 
%               projected point in the canonical space of the triangle
%           N : Array [N1;N2;N3] containing the basis functions at the 
%               parametric location of the projected point in the canonical 
%               space of the triangle in the same order as the vertices of 
%               the triangle
%    isInside : Flag on whether the projection of the point has been found
%               inside the triangle
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the canonical coordinates
%
% 2. Clamp the canonical coordinates given a tolerance
%
% 3. Check the validity of the canonical coordinates
%
% 4. Compute the Constant Strain Triangle (CST) basis functions if the canonical coordinates are found in range
%           
%% Function main body

%% 0. Read input

% Define a tolerance
tol = 1e-3;

% Number of parametric coordinates
noCoord = 2;

% Initialize output variables
isInside = true;
N = zeros(3,1);

%% 1. Compute the canonical coordinates
zeta = [(P1 - P3)'*(P1 - P3) (P2 - P3)'*(P1 - P3)
        (P1 - P3)'*(P2 - P3) (P2 - P3)'*(P2 - P3)]\[(P - P3)'*(P1 - P3)
                                                    (P - P3)'*(P2 - P3)];

%% 2. Clamp the canonical coordinates given a tolerance
for iCoord = 1:noCoord
    if zeta(iCoord,1) < 0 && zeta(iCoord,1) >= -tol
        zeta(iCoord,1) = 0;
    end
    if zeta(iCoord,1) > 1 && zeta(iCoord,1) <= 1 + tol
        zeta(iCoord,1) = 1;
    end
end
                                                
%% 3. Check the validity of the canonical coordinates
for iCoord = 1:noCoord
    if zeta(iCoord,1) < 0 || zeta(iCoord,1) > 1
        isInside = false;
        break;
    end
end

%% 4. Compute the Constant Strain Triangle (CST) basis functions if the canonical coordinates are found in range
if isInside
    N(1:2,1) = zeta(:,1);
    N(3,1) = 1 - sum(N(1:2,1));
end

end
