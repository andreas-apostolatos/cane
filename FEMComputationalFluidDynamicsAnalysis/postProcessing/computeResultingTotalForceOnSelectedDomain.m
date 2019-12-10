function F = computeResultingTotalForceOnSelectedDomain(analysis, nodesDomain, FComplete, parameters)
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
% Calculate force on a single body
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
    
% Reshape the vector to [Vx, Vy, p]
FComplete = -reshape(FComplete, [analysis.noFields,                     ...
             length(FComplete) / analysis.noFields])';

% Select only the forces on the domain
Force = FComplete(nodesDomain, 1:analysis.noSpatialDimensions);

F = zeros(analysis.noSpatialDimensions,1);
for k = 1:analysis.noSpatialDimensions
    F(k,1) = sum(Force(:,k)) * parameters.rho;
end
    
end