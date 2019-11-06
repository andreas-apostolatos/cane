function checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta)
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
%    Input :
%      p,q : polynomial degrees
% mxi,meta : number of knots in u,v-direction
%    nxi,neta : number of control points in u,v-direction
%
% Output :
%   messages on the compatibility
%
%% Function main body

if (nxi+p+1 ~= mxi)
  error('Xi, p and Control points dont match!')
end
if (neta+q+1 ~= meta)
  error('Eta, q and Control points dont match!')
end

end
