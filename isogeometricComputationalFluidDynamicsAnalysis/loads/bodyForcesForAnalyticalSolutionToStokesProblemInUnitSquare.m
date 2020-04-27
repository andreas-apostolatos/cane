function [bx, by] = bodyForcesForAnalyticalSolutionToStokesProblemInUnitSquare ...
    (x, y)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
%   Input :
%     x,y : The locations on the Cartesian domain (unit square)
%
%  Output :
% [bx,by] : The vector of the body forces at the given Cartesian location
%
%% Function main body

% Compute the x-coordinate of the body forces at (x,y)
bx = (12 - 24*y)*x^4 + (-24 + 48*y)*x^3 + (-48*y + 72*y^2 - 48*y^3 + 12)*x^2 + ...
    (-2 + 24*y - 72*y^2 + 48*y^3)*x + 1 - 4*y + 12 * y^2 - 8*y^3;

% Compute the y-coordinate of the body forces at (x,y)
by = (8 - 48*y + 48*y^2)*x^3 + (-12 + 72*y - 72*y^2)*x^2 + ...
    (4 - 24*y + 48*y^2 - 48*y^3 + 24*y^4)*x - 12*y^2 + 24*y^3 - 12*y^4;

end
