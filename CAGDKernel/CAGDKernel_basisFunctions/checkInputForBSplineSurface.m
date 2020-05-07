function checkInputForBSplineSurface ...
    (p, numKnots_xi, numCPs_xi, q, numKnots_eta, numCPs_eta)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% checks compatibility of input parameters for a NURBS surface
%
%                    Input :
%                      p,q : polynomial degrees
% numKnots_xi,numKnots_eta : number of knots in xi,eta-direction
%     numCPs_xi,numCPs_eta : number of control points in xi,eta-direction
%
%                   Output :
%                            message on the compatibility
%
%% Function main body
if (numCPs_xi + p + 1 ~= numKnots_xi)
  error('Xi, p and Control points dont match!')
end
if (numCPs_eta + q + 1 ~= numKnots_eta)
  error('Eta, q and Control points dont match!')
end

end
