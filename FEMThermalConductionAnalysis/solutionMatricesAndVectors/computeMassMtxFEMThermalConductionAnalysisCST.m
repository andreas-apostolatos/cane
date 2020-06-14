function massMtx = computeMassMtxFEMThermalConductionAnalysisCST ...
    (propAnalysis, thrMsh, parameters, gaussInt)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documentation
%
% Returns the mass matrix corresponding to a thermal conduction problem in
% 2D with the pagewise computation of the element stiffness matrices and
% their assembly to a sparse system matrix. This functions expects a
% discretization of the temperature field with the Constant Strain Triangle
% (CST).
%
%               Input :
%        propAnalysis : Information on the analysis type
%                           .type : Analysis type
%              thrMsh : The nodes and the elements of the underlying mesh
%          parameters : The parameters the physical field
%            gaussInt : On the spatial integration
%                           .type : 'default', 'user'
%                     .domainNoGP : Number of Gauss Points for the domain 
%                                   integration
%                   .boundaryNoGP : Number of Gauss Points for the boundary 
%                                   integration            
%            
%              Output :
%             massMtx : The mass matrix of the system
%
% Function layout :
%
% 0. Read input
%
% 1. Create the element freedom tables for all elements at once using the index of the nodes for each element in the array of nodes
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
%    4vi. Compute the stiffness matrix at the Gauss point and add the contribution
% <-
%
% 5. Assemble to the global mass matrix
%
%% Functions main body

%% 0. Read input

% Number of nodes in the mesh
numNodes = length(thrMsh.nodes(:, 1));

% Number of DOFs in the mesh
numDOFs = numNodes;

% Number of nodes at the element level
numNodesEl = 3;

% Number of DOFs at the element level
numDOFsEl = numNodesEl;

% Total number of elements in the mesh
numElmnts = length(thrMsh.elements(:, 1));

% Initialize output
massMtxEl = zeros(numElmnts, numDOFsEl, numDOFsEl);

%% 1. Create the element freedom tables for all elements at once using the index of the nodes for each element in the array of nodes
[~, EFT] = ismember(thrMsh.elements(:, 2:numNodesEl + 1), thrMsh.nodes(:, 1));

%% 2. Get the coordinates of the nodes in a matrix form
nodesI = thrMsh.nodes(EFT(:, 1), 2:end);
nodesJ = thrMsh.nodes(EFT(:, 2), 2:end);
nodesK = thrMsh.nodes(EFT(:, 3), 2:end);

%% 3. Numnerical quadrature
if strcmp(gaussInt.type, 'default')
    numGP = 1;
elseif strcmp(gaussInt.type, 'user')
    numGP = gaussInt.domainNoGP;
end
[GP, GW] = getGaussRuleOnCanonicalTriangle(numGP);

%% 4. Loop over all the quadrature points
for iGP = 1:numGP
    %% 4i. Transform the Gauss Point location from the parameter to the physical space
    xGP = GP(iGP, 1)*nodesI + GP(iGP, 2)*nodesJ + (1 - GP(iGP, 1) - GP(iGP, 2))*nodesK;
    
    %% 4ii. Compute the basis functions and their derivatives at the Gauss Point
    [dN, area] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (nodesI, nodesJ, nodesK, xGP(:, 1, :), xGP(:, 2, :));

    %% 4iii. Compute the determinant of the Jacobian for the transformation from the physical to the parameter spce
    detJxxi = 2*area;
        
	%% 4iv. Form the basis functions matrix at the Gauss Point page-wise
    N = zeros(numElmnts, 1, numDOFsEl);
    for i = 1:numNodesEl
        N(:, 1, i) = dN(:, i, 1);
    end
    
    %% 4vi. Compute the stiffness matrix at the Gauss point and add the contribution
    massMtxEl = massMtxEl + ...
        parameters.cp*parameters.rho*pstimes(pmtimes(ptranspose(N), N)*GW(iGP), detJxxi);
end

%% 5. Assemble to the global mass matrix
[massMtx] = assembleSparseMatricies(EFT', numDOFs, numDOFsEl, massMtxEl);

end