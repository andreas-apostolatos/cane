function lm = fillUpLagrangeMultipliers(p, Xi, CP, isNURBS)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Takes all nessecary arguments that a membrane needs and assigns them to
% the structure
%
%   Input :
%       p : polynomial degrees of the plane surface
%      Xi : knot vectors in xi,eta-directions
%      CP : set of control points and weights
% isNURBS : Flag on whether the geometrical basis is a NURBS or a B-Spline
%
%  Output :
%      lm : structure containing all above information
%
% Function layout :
%
% 1. Assign NURBS parameters
%
% 2. Get the number of Control Points
%
%% Function main body

%% 1. Assign NURBS parameters
lm.p = p;
lm.Xi = Xi;
lm.CP = CP;
lm.isNURBS = isNURBS;

%% 2. Get the number of Control Points
lm.noCPs = length(lm.CP(:, 1));

end
