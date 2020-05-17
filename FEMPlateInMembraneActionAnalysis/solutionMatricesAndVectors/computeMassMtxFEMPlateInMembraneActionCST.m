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
% 1. Get the index of the nodes of each element in the array of nodes
%
% 2. Create the element freedom tables for all elements at once
%
% 3. Get the coordinates of the nodes in a matrix form
%
% 4. Numnerical quadrature
%
% 5. Loop over all the quadrature points
% ->
%    5i. Transform the Gauss Point location from the parameter to the physical space
%
%   5ii. Compute the basis functions and their derivatives at the Gauss Point
%
%  5iii. Compute the determinant of the Jacobian for the transformation from the physical to the parameter spce
%
%   5iv. Form the basis functions matrix at the Gauss Point page-wise
%
%    5v. Compute the element mass matrices at the Gauss point page-wise
% <-
%
% 6. Assemble to the global mass matrix
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
massMtxEl = zeros(numElmnts, numDOFsEl, numDOFsEl);

%% 1. Get the index of the nodes of each element in the array of nodes
[~, idx] = ismember(strMsh.elements(:, 2:numNodesEl + 1), strMsh.nodes(:, 1));

%% 2. Create the element freedom tables for all elements at once
EFT = zeros(numDOFsEl, numElmnts);
for iEFT = 1:numNodesEl
    for iDOFsPerNode = 1:numDOFsPerNode - 1
        EFT(numDOFsPerNode*iEFT, :) = numDOFsPerNode*idx(:, iEFT)';
        EFT(numDOFsPerNode*iEFT - (numDOFsPerNode - iDOFsPerNode), :) = ...
            EFT(numDOFsPerNode*iEFT, :) - (numDOFsPerNode - iDOFsPerNode);
    end
end

%% 3. Get the coordinates of the nodes in a matrix form
nodesI = strMsh.nodes(idx(:, 1), 2:end);
nodesJ = strMsh.nodes(idx(:, 2), 2:end);
nodesK = strMsh.nodes(idx(:, 3), 2:end);

%% 4. Numnerical quadrature
if strcmp(propGaussInt.type, 'default')
    numGP = 2;
elseif strcmp(propGaussInt.type, 'user')
    numGP = propGaussInt.domainNoGP;
end
[GP, GW] = getGaussRuleOnCanonicalTriangle(numGP);

%% 5. Loop over all the quadrature points
for iGP = 1:numGP
    %% 5i. Transform the Gauss Point location from the parameter to the physical space
    xGP = GP(iGP, 1)*nodesI + GP(iGP, 2)*nodesJ + (1 - GP(iGP, 1) - GP(iGP, 2))*nodesK;
    
    %% 5ii. Compute the basis functions and their derivatives at the Gauss Point
    [dN, area] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (nodesI, nodesJ, nodesK, xGP(:, 1, :), xGP(:, 2, :));
    
    %% 5iii. Compute the determinant of the Jacobian for the transformation from the physical to the parameter spce
    detJxxi = 2*area;
    
	%% 5iv. Form the basis functions matrix at the Gauss Point page-wise
    N = zeros(numElmnts, 2, numDOFsEl);
    for i = 1:numNodesEl
        N(:, 1, 2*i - 1) = dN(:, i, 1);
        N(:, 2, 2*i) = dN(:, i, 1);
    end
    
    %% 5v. Compute the element mass matrices at the Gauss point page-wise
    massMtxEl = massMtxEl + ...
           pstimes(pmtimes(ptranspose(N), N)*propParameters.rho*GW(iGP), detJxxi);
end

%% 6. Assemble to the global mass matrix
[massMtx] = assembleSparseMatricies(EFT, numDOFs, numDOFsEl, massMtxEl);

end