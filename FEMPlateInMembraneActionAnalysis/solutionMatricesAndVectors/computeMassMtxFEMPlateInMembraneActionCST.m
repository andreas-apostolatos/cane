function massMtx = computeMassMtxFEMPlateInMembraneActionCST ...
    (propAnalysis, strMsh, propParameters, propGaussInt)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the mass matrix corresponding to a plate in membrane action
% problem using the page-wise construction for a sparse matrix. This
% function expects a spatial discretization using the Constant Strain
% Triangle (CST).
%
%               Input :
%        propAnalysis : Structure containing general information on the 
%                       analysis,
%                           .type : Analysis type
%              strMsh : The nodes and the elements of the underlying mesh
%      propParameters : Structure containing information on the material 
%                       parameters,
%                           .nue : Density
%                             .E : Young's modulus
%                           .nue : Poisson ratio
%                             .t : Thickness
%        propGaussInt : Structure containing information on the numerical 
%                       integration,
%                           .type : 'default', 'user'
%                     .domainNoGP : Number of Gauss Points for the domain 
%                                   integration
%                   .boundaryNoGP : Number of Gauss Points for the boundary 
%                                   integration            
%            
%              Output :
%             massMtx : The mass matrix
%
% Function layout :
%
% 0. Read input
%
% 1. Create the element freedom tables for all elements at once
%
% 2. Get the coordinates of the nodes in a matrix form
%
% 3. Numnerical quadrature
%
% 4. Loop over all the quadrature points
% ->
%    4i. Transform the Gauss Point location from the parameter to the physical space
%
%   4ii. Compute the basis functions and their derivatives at the Gauss Point
%
%  4iii. Compute the determinant of the Jacobian for the transformation from the physical to the parameter spce
%
%   4iv. Form the basis functions matrix at the Gauss Point page-wise
%
%    4v. Compute the element mass matrices at the Gauss point page-wise
% <-
%
% 5. Assemble to the global mass matrix
%
%% Functions main body

%% 0. Read input

% Total number of elements in the mesh
numElmnts = length(strMsh.elements(:, 1));

% Number of nodes in the mesh
numNodes = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
numDOFs = 2*numNodes;

% Number of nodes at the element level
numNodesEl = 3;

% Number of DOFs per node
numDOFsPerNode = 2;

% Number of DOFs at the element level
numDOFsEl = 2*numNodesEl;

% Initialize output
% massMtx = zeros(numDOFs,numDOFs);

massMtxEl = zeros(numElmnts, numDOFsEl, numDOFsEl);

%% 1. Create the element freedom tables for all elements at once
EFT = zeros(numDOFsEl, numElmnts);
for iEFT = 1:numNodesEl
    for iDOFsPerNode = 1:numDOFsPerNode - 1
        EFT(numDOFsPerNode*iEFT, :) = numDOFsPerNode*strMsh.elements(:, iEFT+1)';
        EFT(numDOFsPerNode*iEFT - (numDOFsPerNode - iDOFsPerNode), :) = ...
            EFT(numDOFsPerNode*iEFT, :) - (numDOFsPerNode - iDOFsPerNode);
    end
end

%% 2. Get the coordinates of the nodes in a matrix form
nodes1 = strMsh.nodes(strMsh.elements(:, 2), 2:end);
nodes2 = strMsh.nodes(strMsh.elements(:, 3), 2:end);
nodes3 = strMsh.nodes(strMsh.elements(:, 4), 2:end);

%% 3. Numnerical quadrature
if strcmp(propGaussInt.type, 'default')
    numGP = 2;
elseif strcmp(propGaussInt.type, 'user')
    numGP = propGaussInt.domainNoGP;
end
[GP, GW] = getGaussRuleOnCanonicalTriangle(numGP);

%% 4. Loop over all the quadrature points
for iGP = 1:numGP
    %% 4i. Transform the Gauss Point location from the parameter to the physical space
    xGP = GP(iGP, 1)*nodes1 + GP(iGP, 2)*nodes2 + (1 - GP(iGP, 1) - GP(iGP, 2))*nodes3;
    
    %% 4ii. Compute the basis functions and their derivatives at the Gauss Point
    [dN, area] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (nodes1, nodes2, nodes3, xGP(:, 1, :), xGP(:, 2, :));
    
    %% 4iii. Compute the determinant of the Jacobian for the transformation from the physical to the parameter spce
    detJxxi = 2*area;
    
	%% 4iv. Form the basis functions matrix at the Gauss Point page-wise
    N = zeros(numElmnts, 2, numDOFsEl);
    for i = 1:numNodesEl
        N(:, 1, 2*i - 1) = dN(:, i, 1);
        N(:, 2, 2*i) = dN(:, i, 1);
    end
    
    %% 4v. Compute the element mass matrices at the Gauss point page-wise
    massMtxEl = massMtxEl + ...
           pstimes(pmtimes(ptranspose(N), N)*propParameters.rho*GW(iGP), detJxxi);
end

%% 5. Assemble to the global mass matrix
[massMtx] = assembleSparseMatricies(EFT, numDOFs, numDOFsEl, massMtxEl);

end