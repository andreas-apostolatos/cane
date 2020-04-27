function checkInputForBSplineCurve(p, numKnots, numCPs)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% checks compatibility of input parameters for a NURBS line
%
%    Input :
%        p : polynomial order
% numKnots : number of knots
%   numCPs : number of control points
%
% Output :
%          messages on the compatibility
%
%% Function main body
if (numCPs + p + 1 ~= numKnots)
  error('U, p and Control points dont match!')
end

end
