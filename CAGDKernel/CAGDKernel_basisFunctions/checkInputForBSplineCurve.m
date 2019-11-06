function checkInputForBSplineCurve(p,m,n)
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
%  Input :
%    p,q : polynomial degrees
%  mu,mv : number of knots in u,v-direction
%  nu,nv : number of control points in u,v-direction
%
% Output :
%          messages on the compatibility
%
%% Function main body

if (n+p+1 ~= m)
  error('U, p and Control points dont match!')
end

end
