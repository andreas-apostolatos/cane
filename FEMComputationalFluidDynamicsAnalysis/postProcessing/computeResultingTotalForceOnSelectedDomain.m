function F = computeResultingTotalForceOnSelectedDomain ...
    (analysis, nodesDomain, FComplete, parameters)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% Return total resulting force acting on the selected boundary.
%
%              Input :
%           analysis : Analysis type and number of dimensions
%        nodesDomain : Global numbering of the nodes on the selected domain
%          FComplete : The complete force vector
%         parameters : Flow parameters
%   
%             Output :
%                  F : Vector of noSpatialDimensions containing the 
%                      components of the total resulting force along each 
%                      spatial dimension over the selected domain
%
%% Function main body
F = zeros(analysis.noSpatialDimensions,1);
DOFsDomain = analysis.noFields * nodesDomain - analysis.noSpatialDimensions;
for k = 1:analysis.noSpatialDimensions
    activeDOFsDomain = DOFsDomain + (k-1);
    F(k,1) = -parameters.rho * sum(FComplete(activeDOFsDomain(1:end)));
end

end