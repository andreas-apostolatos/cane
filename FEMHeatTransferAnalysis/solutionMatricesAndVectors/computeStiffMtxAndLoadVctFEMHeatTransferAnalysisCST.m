function [K, F, minElEdgeSize] = ...
    computeStiffMtxAndLoadVctFEMHeatTransferAnalysisCST ...
    (propAnalysis, u, uSaved, uDot, uDotSaved, uMeshALE, ...
    precompStiffMtx, precomResVct, DOFNumbering, strMsh, F, ...
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
% the 2D heat transfer analysis using the Constant Strain Triangle
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
%            strMsh : Nodes and elements in the mesh
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
% 1. Create the element freedom tables for all elements at once
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
numNodes = length(strMsh.nodes(:, 1));

% Number of DOFs in the mesh
numDOFs = numNodes;

% Number of nodes at the element level
numNodesEl = 3;

% Number of degrees of freedom per element
numDOFsEl = numNodesEl;

% Total number of elements in the mesh
numElmnts = length(strMsh.elements(:, 1));

% Initialize output arrays
stiffMtxEl = zeros(numElmnts, numDOFsEl, numDOFsEl);

%% 1. Create the element freedom tables for all elements at once
EFT = zeros(numDOFsEl, numElmnts);
for iEFT = 1:numNodesEl
        EFT(iEFT, :) = strMsh.elements(:, iEFT)';
end

%% 2. Get the coordinates of the nodes in a matrix form

% define function to calculate euclidean norm
euclideanNorm = @(nodes) sqrt(nodes(:, 1, 1).^2 + nodes(:, 2, 1).^2 + ...
    nodes(:, 3, 1).^2);

% Get the nodes of the mesh
nodes1 = strMsh.nodes(strMsh.elements(:, 1), :);
nodes2 = strMsh.nodes(strMsh.elements(:, 2), :);
nodes3 = strMsh.nodes(strMsh.elements(:, 3), :);

% get element sizes
h = min( [euclideanNorm(nodes1 - nodes2) euclideanNorm(nodes1 - nodes3) ...
          euclideanNorm(nodes2 - nodes3)], [], 2);

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
    xGP = GP(iGP, 1)*nodes1 + GP(iGP, 2)*nodes2 + (1 - GP(iGP, 1) - GP(iGP, 2))*nodes3;
    
    %% 5ii. Compute the basis functions and their derivatives at the Gauss Point
    [dN, area] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (nodes1, nodes2, nodes3, xGP(:, 1, :), xGP(:, 2, :));
        
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
    stiffMtxEl = stiffMtxEl + propParameters.k * pstimes(pmtimes(ptranspose(B), B)*GW(iGP), detJxxi);
end

%% 6. Assemble to the global system matrix
K = assembleSparseMatricies(EFT, numDOFs, numDOFsEl, stiffMtxEl);

end