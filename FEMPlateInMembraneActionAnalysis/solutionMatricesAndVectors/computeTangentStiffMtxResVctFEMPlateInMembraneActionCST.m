function [tanMtx, resVct, minElEdgeSize] = ...
    computeTangentStiffMtxResVctFEMPlateInMembraneActionCST ...
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
% Returns the tangent stiffness matrix and residue vector corresponding 
% to the geometrically nonlinear plate in membrane action analysis using 
% the Constant Strain Triangle (CST) for the displacement field 
% discretization.
%
%               Input :
%        propAnalysis : Structure containing general information on the
%                       analysis,
%                           .type : Analysis type
%                   u : The discrete solution field of the current time  
%                       step
%              uSaved : The discrete solution field of the previous time 
%                       step
%                uDot : The time derivative of the discrete solution field 
%                       of the current time step
%           uDotSaved : The time derivative of the discrete solution field 
%                       of the previous time step
%            uMeshALE : Dummy array for this function
%     precompStiffMtx : constant part of the stiffness matrix which can be
%                       precomputed
%        precomResVct : Constant part of the residual vector which can be
%                       precomputed
%        DOFNumbering : The global numbering of the DOFs
%              strMsh : The nodes and the elements of the underlying mesh
%                   F : Global load vector corresponding to surface tractions
%          loadFactor : The load factor for the nonlinear steps
% computeBodyForceVct : Function handle to body force vector computation
%     propStrDynamics : Transient analysis parameters:
%                         .method : Time integration method
%                             .T0 : Start time of the simulation
%                           .TEnd : End time of the simulation
%                             .nT : Number of time steps
%                             .dt : Time step (numeric or adaptive)
%                   t : The current time of the transient simulation
%      propParameters : The parameters the physical field
%            gaussInt : Structure containing information on the numerical 
%                       integration
%                           .type : 'default', 'user'
%                     .domainNoGP : Number of Gauss Points for the domain 
%                                   integration
%                   .boundaryNoGP : Number of Gauss Points for the boundary 
%                                   integration            
%
%          Output :
%          tanMtx : The tangential stiffness matrix of the system
%          resVct : The residual vector
%   minElEdgeSize : The minimum element size in the mesh
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
% 5. Compute the material matrices for each element
%
% 6. Loop over all the quadrature points
% ->
%    6i. Transform the Gauss Point location from the parameter to the physical space
%
%   6ii. Compute the basis functions and their derivatives at the Gauss Point
%
%  6iii. Compute the determinant of the Jacobian for the transformation from the physical to the parameter spce
%
%   6iv. Compute the deformation gradient tensor F at the Gauss point
%  
%    6v. Form the B-Operator matrix corresponding to the current configuration at the Gauss point page-wise
%
%   6vi. Compute the material stifness matrix at the Gauss point page-wise
%
%  6vii. Compute the Green-Lagrange strain tensor (page-wise) and re-arrange it to the Voigt notation
%
% 6viii. Compute the Cauchy stress tensor (page-wise) out of the Voigt Green-Lagrange strain vector
%
%   6ix. Compute the element geometric stiffness matrix at the Gauss point page-wise
%
%    6x. Compute the internal residual K(u)*u at the Gauss point page-wise
%
%   6xi. Form the basis functions matrix at the Gauss Point page-wise
%
%  6xii. Compute the body force vector at the Gauss point page-wise
% <-
%
% 7. Assemble the global system matrices and vectors
%
% 8. Update the force vector with the body force contributions
%
%% Functions main body

%% 0. Read input

% Number of nodes in the mesh
numNodes = length(strMsh.nodes(:, 1));

% Number of DOFs in the mesh
numDOFs = 2*numNodes;

% Number of nodes at the element level
numNodesEl = 3;

% Number of DOFs per node
numDOFsPerNode = 2;

% Number of degrees of freedom per element
numDOFsEl = numDOFsPerNode*numNodesEl;

% Total number of elements in the mesh
numElmnts = length(strMsh.elements(:, 1));

% Initialize output arrays
tanMtxMatEl = zeros(numElmnts, numDOFsEl, numDOFsEl);
tanMtxGeoEl = zeros(numElmnts, numDOFsEl, numDOFsEl);
FBodyEl = zeros(numElmnts, numDOFsEl, 1);
resIntEl = zeros(numElmnts, numDOFsEl, 1);

%% 1. Create the element freedom tables for all elements at once
EFT = zeros(numDOFsEl, numElmnts);
for iEFT = 1:numNodesEl
    for counterDOFsPerNode = 1:numDOFsPerNode - 1
        EFT(numDOFsPerNode*iEFT, :) = numDOFsPerNode*strMsh.elements(:, iEFT + 1)';
        EFT(numDOFsPerNode*iEFT - (numDOFsPerNode - counterDOFsPerNode), :) = ...
            EFT(numDOFsPerNode*iEFT, :) - (numDOFsPerNode - counterDOFsPerNode);
    end
end

%% 2. Get the coordinates of the nodes in a matrix form

% define function to calculate euclidean norm
euclideanNorm = @(nodes) sqrt(nodes(:, 1, 1).^2 + nodes(:, 2, 1).^2 + ...
    nodes(:, 3, 1).^2);

% Get the nodes of the mesh
nodes1 = strMsh.nodes(strMsh.elements(:, 2), 2:end);
nodes2 = strMsh.nodes(strMsh.elements(:, 3), 2:end);
nodes3 = strMsh.nodes(strMsh.elements(:, 4), 2:end);

% get element sizes
h = min([euclideanNorm(nodes1 - nodes2) euclideanNorm(nodes1 - nodes3) ...
    euclideanNorm(nodes2 - nodes3)], [], 2);

%% 3. Get the minimum element edge size
minElEdgeSize = min(h);

%% 4. Numerical quadrature
if strcmp(propInt.type, 'default')
    numGP = 2;
elseif strcmp(propInt.type, 'user')
    numGP = propInt.domainNoGP;
end
[GP, GW] = getGaussRuleOnCanonicalTriangle(numGP);

%% 5. Compute the material matrices for each element
C = zeros(numElmnts, 3, 3);
if strcmp(propAnalysis.type, 'planeStress')
    preFactor = propParameters.E/(1-propParameters.nue^2);
    CEl = preFactor*[1              propParameters.nue 0
                     propParameters.nue 1              0
                     0              0             (1 - propParameters.nue)/2];
elseif strcmp(propAnalysis.type, 'planeStrain')
    preFactor = propParameters.E*(1 - propParameters.nue)/(1+propParameters.nue)/(1 - 2*propParameters.nue);
    CEl = preFactor*[1                                   propParameters.nue/(1 - propParameters.nue) 0
                     propParameters.nue/(1 - propParameters.nue) 1                                   0
                     0                                   0                                   (1 - 2*propParameters.nue)/2/(1 - propParameters.nue)];
else
    error('Select a valid analysis type in analysis.type');
end

% Assemble the material matrices page-wise
for iElmnts = 1:numElmnts
    C(iElmnts, :, :) = CEl;
end

%% 6. Loop over all the quadrature points
for iGP = 1:numGP
    %% 6i. Transform the Gauss Point location from the parameter to the physical space
    xGP = GP(iGP, 1)*nodes1 + GP(iGP, 2)*nodes2 + (1 - GP(iGP, 1) - GP(iGP, 2))*nodes3;
    
    %% 6ii. Compute the basis functions and their derivatives at the Gauss Point
    [dN, area] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (nodes1, nodes2, nodes3, xGP(:, 1, :), xGP(:, 2, :));
    
    %% 6iii. Compute the determinant of the Jacobian for the transformation from the physical to the parameter spce
    detJxxi = 2*area;
    
    %% 6iv. Compute the deformation gradient tensor F at the Gauss point
    
    % Get the derivates of the basis functions
    dNdX = dN(:, :, 2:3);
    
    % Get the displacement vector
    uEl = u(EFT');
    uEl = reshape(uEl,numElmnts, 2, 3);
    
    % Build an identity (unit) matrix
    I_Mtx = zeros(numElmnts, 2, 2);
    I_Mtx(:, 1, 1) = ones(numElmnts, 1);
    I_Mtx(:, 2, 2) = ones(numElmnts, 1);   
    
    % Compute the deformation gradient tensor
    FDefGrad = I_Mtx + pmtimes(uEl, dNdX);
    
    %% 6v. Form the B-Operator matrix corresponding to the current configuration at the Gauss point page-wise
    B = zeros(numElmnts, 3, numDOFsEl);
    for i = 1:numNodesEl
        B(:, 1, 2*i - 1) = pmtimes( dNdX(:, i, 1), FDefGrad(:,1,1) );
        B(:, 2, 2*i - 1) = pmtimes( dNdX(:, i, 2), FDefGrad(:,1,2) );
        B(:, 3, 2*i - 1) = pmtimes( dNdX(:, i, 1), FDefGrad(:,1,2) ) + ...
                           pmtimes( dNdX(:, i, 2), FDefGrad(:,1,1) );
                       
        B(:, 1, 2*i) = pmtimes( dNdX(:, i, 1), FDefGrad(:,2,1) );
        B(:, 2, 2*i) = pmtimes( dNdX(:, i, 2), FDefGrad(:,2,2) );
        B(:, 3, 2*i) = pmtimes( dNdX(:, i, 1), FDefGrad(:,2,2) ) + ...
                       pmtimes( dNdX(:, i, 2), FDefGrad(:,2,1) );
    end
    
    %% 6vi. Compute the material stifness matrix at the Gauss point page-wise
    tanMtxMatEl = tanMtxMatEl + ...
           pstimes(pmtimes(pmtimes(ptranspose(B), C), B)*GW(iGP), detJxxi);
    
    %% 6vii. Compute the Green-Lagrange strain tensor (page-wise) and re-arrange it to the Voigt notation
    epsilonGLTensor = 0.5*(pmtimes(ptranspose(FDefGrad), FDefGrad) - I_Mtx);
    epsilonGLVoigt = [epsilonGLTensor(:, 1, 1), epsilonGLTensor(:, 2, 2), 2*epsilonGLTensor(:, 1, 2)];
    
    %% 6viii. Compute the Cauchy stress tensor (page-wise) out of the Voigt Green-Lagrange strain vector
    stressCauchyVoigt = pmtimes(C, epsilonGLVoigt);
    
    % Assemble the Cauchy stress tensor
    stressCauchyTensor(:, 2, 2) = stressCauchyVoigt(:, 2);
    stressCauchyTensor(:, 2, 1) = stressCauchyVoigt(:, 3);
    stressCauchyTensor(:, 1, 2) = stressCauchyVoigt(:, 3);
    stressCauchyTensor(:, 1, 1) = stressCauchyVoigt(:, 1);  
   
    %% 6ix. Compute the element geometric stiffness matrix at the Gauss point page-wise
    
    % Compute the H matrix
    H = pmtimes(pmtimes(dNdX, stressCauchyTensor), ptranspose(dNdX));
    
    % Assemble the temporary geometric stifness matrix 
    tanMtxGeoElGP = zeros(numElmnts, numDOFsEl, numDOFsEl);
    for i = 1:numNodesEl
        tanMtxGeoElGP(:, 1, 2*i - 1) = H(:, 1, i);
        tanMtxGeoElGP(:, 3, 2*i - 1) = H(:, 2, i);
        tanMtxGeoElGP(:, 5, 2*i - 1) = H(:, 3, i);
        
        tanMtxGeoElGP(:, 2, 2*i) = H(:, 1, i);
        tanMtxGeoElGP(:, 4, 2*i) = H(:, 2, i);
        tanMtxGeoElGP(:, 6, 2*i) = H(:, 3, i);
    end
    
    % Compute the element geometric stiffness matrix
    tanMtxGeoEl = tanMtxGeoEl + pstimes((tanMtxGeoElGP*GW(iGP)),detJxxi);
    
    %% 6x. Compute the internal residual K(u)*u at the Gauss point page-wise
    resIntEl = resIntEl + ...
        pstimes(pmtimes(ptranspose(B), stressCauchyVoigt)*GW(iGP), detJxxi);

	%% 6xi. Form the basis functions matrix at the Gauss Point page-wise
    N = zeros(numElmnts, 2, numDOFsEl);
    for i = 1:numNodesEl
        N(:, 1, 2*i - 1) = dN(:, i, 1);
        N(:, 2, 2*i) = dN(:, i, 1);
    end
    
    %% 6xii. Compute the body force vector at the Gauss point page-wise
    bF = computeBodyForces(xGP(:, 1), xGP(:, 2), xGP(:, 3));
    FBodyEl = FBodyEl + ...
        pstimes(pmtimes(ptranspose(N), ptranspose(bF(:, :, 1:2)))*GW(iGP), detJxxi);
end

%% 7. Assemble the global system matrices and vectors
[tanMtx] = assembleSparseMatricies ...
    (EFT, numDOFs, numDOFsEl, tanMtxMatEl, tanMtxGeoEl);
[resVctInt, FBody] = assembleSparseVectors(EFT, numDOFs, numDOFsEl, resIntEl, FBodyEl);

%% 8. Update the force vector with the body force contributions
resVct = resVctInt - loadFactor*(FBody + F);

end