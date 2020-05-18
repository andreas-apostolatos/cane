function F = computeForceVctOnSelectedDomain ...
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
% Return the force vector on the nodes of the selected boundary.
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

DOFs = reshape(vertcat(vertcat(analysis.noFields*nodesDomain' - (analysis.noFields - 1), ...
    analysis.noFields*nodesDomain' - (analysis.noFields - 2)), ...
    analysis.noFields*nodesDomain' - (analysis.noFields - 3)), 1, [])';
F = - parameters.rho*FComplete(DOFs, 1);

end