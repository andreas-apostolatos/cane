function [K, F, minElEdgeSize] = ...
    computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST ...
    (propAnalysis, u, uSaved, uDot, uDotSaved, precompStiffMtx, ...
    precomResVct, DOFNumbering, strMsh, F, loadFactor, computeBodyForces, ...
    propStrDynamics, parameters, propInt)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the stiffness matrix and the load vector corresponding to the
% plate in membrane action analysis using the Constant Strain Triangle
% (CST) for the displacement field discretization.
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
%   precompStiffMtx : constant part of the stiffness matrix which can be
%                     precomputed
%      precomResVct : Constant part of the residual vector which can be
%                     precomputed
%            strMsh : Nodes and elements in the mesh
%                 F : The global load vector corresponding to surface
%                     tractions
%        loadFactor : Load factor, nessecary for nonlinear computations 
%                     (dummy variable for this function)
% computeBodyForces : Function handle to body force vector computation
%   propStrDynamics : Structure containing information on the time
%                     integration regarding the structural dynamics
%        parameters : Problem specific technical parameters
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
% 4. Numnerical quadrature
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
%   6iv. Form the basis functions matrix at the Gauss Point page-wise
%
%    6v. Form the B-Operator matrix for the plate in membrane action problem page-wise
%
%   6vi. Compute the element load vector due to body forces and add the contribution
%
%  6vii. Compute the stiffness matrix at the Gauss point and add the contribution
% <-
%
% 7. Add the contribution from the Gauss Point and assemble to the global system
%
% 8. Update the force vector with the body force contributions
%
%% Functions main body

%% 0. Read input

% Number of nodes in the mesh
noNodes = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
noDOFs = 2*noNodes;

% Number of nodes at the element level
noNodesEl = 3;

% Number of DOFs per node
noDOFsPerNode = 2;

% Number of degrees of freedom per element
noDOFsEl = noDOFsPerNode*noNodesEl;

% Total number of elements in the mesh
noElmnts = length(strMsh.elements(:,1));

% Initialize output arrays
KEl = zeros(noElmnts,noDOFsEl,noDOFsEl);
FBodyEl = zeros(noElmnts,noDOFsEl,1);

%% 1. Create the element freedom tables for all elements at once
EFT = zeros(noDOFsEl,noElmnts);
for counterEFT = 1:noNodesEl
    for counterDOFsPerNode = 1:noDOFsPerNode-1
        EFT(noDOFsPerNode*counterEFT, :) = noDOFsPerNode*strMsh.elements(:,counterEFT)';
        EFT(noDOFsPerNode*counterEFT-(noDOFsPerNode-counterDOFsPerNode), :) = ...
            EFT(noDOFsPerNode*counterEFT, :) - (noDOFsPerNode-counterDOFsPerNode);
    end
end

%% 2. Get the coordinates of the nodes in a matrix form

% define function to calculate euclidean norm
euclideanNorm = @(nodes) sqrt( nodes(:,1,1).^2 + nodes(:,2,1).^2 + nodes(:,3,1).^2 );

% Get the nodes of the mesh
nodes1 = strMsh.nodes(strMsh.elements(:,1),:);
nodes2 = strMsh.nodes(strMsh.elements(:,2),:);
nodes3 = strMsh.nodes(strMsh.elements(:,3),:);

% get element sizes
h = min( [ euclideanNorm(nodes1-nodes2) euclideanNorm(nodes1-nodes3) ...
           euclideanNorm(nodes2-nodes3)],[],2);

%% 3. Get the minimum element edge size
minElEdgeSize = min(h);

%% 4. Numnerical quadrature
if strcmp(propInt.type,'default')
    noGP = 1;
elseif strcmp(propInt.type,'user')
    noGP = propInt.nGP;
end
[GP,GW] = getGaussRuleOnCanonicalTriangle(noGP);

%% 5. Compute the material matrices for each element
C = zeros(noElmnts,3,3);
if strcmp(propAnalysis.type,'planeStress')
    preFactor = parameters.E/(1-parameters.nue^2);
    CEl = preFactor*[1              parameters.nue 0
                     parameters.nue 1              0
                     0              0             (1-parameters.nue)/2];
elseif strcmp(propAnalysis.type,'planeStrain')
    preFactor = parameters.E*(1-parameters.nue)/(1+parameters.nue)/(1-2*parameters.nue);
    CEl = preFactor*[1                                 parameters.nue/(1-parameters.nue) 0
                     parameters.nue/(1-parameters.nue) 1                                 0
                     0                                 0                                 (1-2*parameters.nue)/2/(1-parameters.nue)];
end
for loop = 1:noElmnts
    C(loop,:,:) = CEl;
end

%% 6. Loop over all the quadrature points
for counterGP = 1:noGP
    %% 6i. Transform the Gauss Point location from the parameter to the physical space
    xGP = GP(counterGP,1)*nodes1 + GP(counterGP,2)*nodes2 + (1-GP(counterGP,1)-GP(counterGP,2))*nodes3;
    
    %% 6ii. Compute the basis functions and their derivatives at the Gauss Point
    [dN,area] = computeCST2DBasisFunctionsAndFirstDerivatives...
            (nodes1,nodes2,nodes3,xGP(:,1,:),xGP(:,2,:));
        
    %% 6iii. Compute the determinant of the Jacobian for the transformation from the physical to the parameter spce
    detJxxi = 2*area;
        
	%% 6iv. Form the basis functions matrix at the Gauss Point page-wise
    N = zeros(noElmnts,2,noDOFsEl);
    for i = 1:noNodesEl
        N(:,1,noDOFsPerNode*i-noDOFsPerNode+1) = dN(:,i,1);
        N(:,2,noDOFsPerNode*i-noDOFsPerNode+2) = dN(:,i,1);
    end
    
    %% 6v. Form the B-Operator matrix for the plate in membrane action problem page-wise
    B = zeros(noElmnts,3,noDOFsEl);
    for i = 1:noNodesEl
        B(:,1,2*i-1) = dN(:,i,2);
        B(:,2,2*i) = dN(:,i,3);
        B(:,3,2*i-1) = dN(:,i,3);
        B(:,3,2*i) = dN(:,i,2);
    end
    
    %% 6vi. Compute the element load vector due to body forces and add the contribution
    bF = computeBodyForces(xGP(:,1),xGP(:,2),xGP(:,3));
    FBodyEl = FBodyEl + pstimes(pmtimes(ptranspose(N),ptranspose(bF(:,:,1:2)))*GW(counterGP),detJxxi);
    
    %% 6vii. Compute the stiffness matrix at the Gauss point and add the contribution
    KEl = KEl + pstimes(pmtimes(pmtimes(ptranspose(B),C),B)*GW(counterGP),detJxxi);
end

%% 7. Assemble to the global system matrices
% [K,FBody] = assembleSparseMatricies(EFT,noDOFs,noDOFsEl,KEl,FBodyEl);
[K] = assembleSparseMatricies(EFT,noDOFs,noDOFsEl,KEl);

%% 8. Update the force vector with the body force contributions
% F = F + FBody;

end
