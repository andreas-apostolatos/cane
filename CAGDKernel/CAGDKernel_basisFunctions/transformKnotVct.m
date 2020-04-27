function Xi = transformKnotVct(Xi, targetDomain)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the transformed to the targetDomain knot vector Xi.
%
%        Input :
%           Xi : The knot vector to be transformed
% targetDomain : Vector [targetDomainStart targetDomainEnd] indicating the
%                target space to where knot vector Xi needs to be
%                transformed
%
%       Output :
%           Xi : The transformed knot vector
%
%% Function main body
XiStart = Xi(1);
XiEnd = Xi(end);
for i = 1:length(Xi)
    Xi(i) = targetDomain(1,2) + (Xi(i) - XiEnd)*(targetDomain(1,1) - targetDomain(1,2))/(XiStart - XiEnd);
end

end
