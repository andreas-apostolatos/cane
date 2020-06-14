function [K, F, minElEdgeSize] = ...
    computeStiffMtxAndLoadVctFEMThermalConductionAnalysisCST ...
    (propAnalysis, u, uSaved, uDot, uDotSaved, uMeshALE, ...
    precompStiffMtx, precomResVct, DOFNumbering, thrMsh, F, ...
    loadFactor, computeBodyForces, propStrDynamics, t, ...
    propParameters, propInt)
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
% Returns the stiffness matrix and the flux (load) vector corresponding to 
% the 2D thermal conduction analysis using the Constant Strain Triangle
% (CST) for the temperature field discretization.
%
%             Input :
%      propAnalysis : Structure containing general information about the 
%                     analysis,
%                           .type : The analysis type
%                 u : The discrete solution field of the previous nonlinear
%                     iteration
%            uSaved : The discrete solution field of the previous time step
%              uDot : The time derivative of the discrete solution field of 
%                     the previous nonlinear iteration
%         uDotSaved : The time derivative of the discrete solution field of 
%                     the previous time step
%          uMeshALE : Dummy variable for this function
%   precompStiffMtx : constant part of the stiffness matrix which can be
%                     precomputed
%      precomResVct : Constant part of the residual vector which can be
%                     precomputed
%      DOFNumbering : The global numbering of the DOFs
%            thrMsh : Nodes and elements in the mesh
%                 F : The global load vector corresponding to surface
%                     tractions
%        loadFactor : Load factor, nessecary for nonlinear computations 
%                     (dummy variable for this function)
% computeBodyForces : Function handle to body force vector computation
%   propStrDynamics : Structure containing information on the time
%                     integration regarding the structural dynamics
%    propParameters : Problem specific technical parameters
%           propInt : Structure containing information on the quadrature
%
%            Output :
%                 K : The master stiffness matrix of the system
%                 F : The updated load vector accounting also for the body
%                     forces
%     minElEdgeSize : The minimum element edge size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Create the element freedom tables for all elements at once using the index of the nodes for each element in the array of nodes
%
% 2. Get the coordinates of the nodes in a matrix form
%
% 3. Get the minimum element edge size
%
% 4. Numerical quadrature
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
%    5v. Form the B-Operator matrix for the plate in membrane action problem page-wise
%
%   5vi. Compute the stiffness matrix at the Gauss point and add the contribution
% <-
%
% 6. Assemble to the global system matrix
%
%% Functions main body

%% 0. Read input

% Number of nodes in the mesh
numNodes = length(thrMsh.nodes(:, 1));

% Number of DOFs in the mesh
numDOFs = numNodes;

% Number of nodes at the element level
numNodesEl = 3;

% Number of degrees of freedom per element
numDOFsEl = numNodesEl;

% Total number of elements in the mesh
numElmnts = length(thrMsh.elements(:, 1));

% Initialize output arrays
stiffMtxEl = zeros(numElmnts, numDOFsEl, numDOFsEl);

%% 1. Create the element freedom tables for all elements at once using the index of the nodes for each element in the array of nodes
[~, EFT] = ismember(thrMsh.elements(:, 2:numNodesEl + 1), thrMsh.nodes(:, 1));

%% 2. Get the coordinates of the nodes in a matrix form

% define function to calculate euclidean norm
euclideanNorm = @(nodes) sqrt(nodes(:, 1, 1).^2 + nodes(:, 2, 1).^2 + ...
    nodes(:, 3, 1).^2);

% Get the nodes of the mesh
nodesI = thrMsh.nodes(EFT(:, 1), 2:end);
nodesJ = thrMsh.nodes(EFT(:, 2), 2:end);
nodesK = thrMsh.nodes(EFT(:, 3), 2:end);

% get element sizes
h = min( [euclideanNorm(nodesI - nodesJ) euclideanNorm(nodesI - nodesK) ...
          euclideanNorm(nodesJ - nodesK)], [], 2);

%% 3. Get the minimum element edge size
minElEdgeSize = min(h);

%% 4. Numerical quadrature
if strcmp(propInt.type, 'default')
    numGP = 1;
elseif strcmp(propInt.type, 'user')
    numGP = propInt.domainNoGP;
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
    N = zeros(numElmnts, 1, numDOFsEl);
    for i = 1:numNodesEl
        N(:, 1, i) = dN(:, i, 1);
    end
    
    %% 5v. Form the B-Operator matrix for the plate in membrane action problem page-wise
    B = zeros(numElmnts, 2, numDOFsEl);
    for i = 1:numNodesEl
        B(:, 1, i) = dN(:, i, 2);
        B(:, 2, i) = dN(:, i, 3);
    end
    
    %% 5vi. Compute the stiffness matrix at the Gauss point and add the contribution
    stiffMtxEl = stiffMtxEl + propParameters.k*pstimes(pmtimes(ptranspose(B), B)*GW(iGP), detJxxi);
end

%% 6. Assemble to the global system matrix
K = assembleSparseMatricies(EFT', numDOFs, numDOFsEl, stiffMtxEl);

end