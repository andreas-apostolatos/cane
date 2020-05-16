function [tanMtx, resVct, minElAreaSize] = ...
    computeTangentStiffMtxResVctFEMPlateInMembraneAction...
    (propAnalysis, u, uSaved, uDot, uDotSaved, uMeshALE, ...
    precompStiffMtx, precomResVct, DOFNumbering, strMsh, F, ...
    loadFactor, computeBodyForceVct, propStrDynamics, t, ...
    propParameters, propGaussInt)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangential stiffness matrix and residue vector corresponding 
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
%                    .noTimeSteps : Number of time steps
%                             .dt : Time step (numeric or adaptive)
%                   t : The current time of the transient simulation
%      propParameters : The parameters the physical field
%        propGaussInt : Structure containing information on the numerical 
%                       integration,
%                           .type : 'default', 'user'
%                     .domainNoGP : Number of Gauss Points for the domain 
%                                   integration
%                   .boundaryNoGP : Number of Gauss Points for the boundary 
%                                   integration            
%
%          Output :
%               K : The tangential stiffness matrix of the system
%          resVct : The residual vector
%      minElASize : The minimum element area size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Numnerical quadrature
%
% 2. Loop over all the elements in the mesh
% ->
%    2i. Get the element in the mesh
%
%   2ii. Get the nodes in the element
%
%  2iii. Create an Element Freedom Table (EFT)
%
%   2iv. Re-arrange the element displacement vector as [u1x u1y; u2x u2y; u3x u3y]
%
%    2v. Loop over the quadrature points
%    ->
%        2v.1. Transform the Gauss Point location from the parameter to the physical space
%
%        2v.2. Compute the basis functions and their derivatives at the Gauss Point
%
%        2v.3. Compute the determinant of the Jacobian transformation from the physical to the parent space
%
%        2v.4. Compute the external load vector due to body forces
%
%        2v.5. Compute the material, the geometric, the mass matrix, the internal residual and external body force vector at the Gauss point
%
%        2v.6. Assemble the local matrices and vectors to the global ones via the EFT
%    <-
% <-
% 
% 3. Compute the stiffness matrix and the residual vector
%
%% Functions main body

%% 0. Read input

% Number of nodes in the mesh
numNodes = length(strMsh.nodes(:, 1));

% Number of DOFs in the mesh
numDOFs = 2*numNodes;

% Number of nodes at the element level
numNodesEl = 3;

% Number of DOFs at the element level
numDOFsEl = 2*numNodesEl;

% Initialize the minimum element size
element1 = strMsh.elements(1, 2:end);
node1 = strMsh.nodes(element1(1, 1), 2:end);
node2 = strMsh.nodes(element1(1, 2), 2:end);
node3 = strMsh.nodes(element1(1, 3), 2:end);
minElAreaSize = 0.5*abs((node3(1, 1) - node2(1, 1))*(node1(1, 2) - node2(1, 2)) - ...
    (node1(1, 1) - node2(1, 1))*(node3(1, 2) - node2(1, 2)));

% Compute the material matrix for the given problem
if strcmp(propAnalysis.type, 'planeStress')
    preFactor = propParameters.E/(1 - propParameters.nue^2);
    C = preFactor*[1                  propParameters.nue 0
                   propParameters.nue 1              0
                   0                  0             (1 - propParameters.nue)/2];
elseif strcmp(propAnalysis.type, 'planeStrain')
    preFactor = propParameters.E*(1 - propParameters.nue)/(1 + propParameters.nue)/(1 - 2*propParameters.nue);
    C = preFactor*[1                                           propParameters.nue/(1 - propParameters.nue) 0
                   propParameters.nue/(1 - propParameters.nue) 1                                 0
                   0                                           0                                 (1 - 2*propParameters.nue)/2/(1 - propParameters.nue)];
else
    error('Select a valid analysis type in analysis.type');
end

% Initialize global master stiffness matrix
tanMtxMat = zeros(numDOFs, numDOFs);

% Initialize global geometric stiffness matrix
tanMtxGeo = zeros(numDOFs, numDOFs);

% Initialize the residual vector
resVctInt = zeros(numDOFs, 1);

% Initialize the body force vector
Fbody = zeros(numDOFs, 1);

%% 1. Numnerical quadrature
if strcmp(propGaussInt.type, 'default')
    numGP = 2;
elseif strcmp(propGaussInt.type, 'user')
    numGP = propGaussInt.domainNoGP;
end
[GP, GW] = getGaussRuleOnCanonicalTriangle(numGP);

%% 2. Loop over all the elements in the mesh
for iEl = 1:length(strMsh.elements(:, 1))
    %% 2i. Get the element in the mesh
    element = strMsh.elements(iEl, 2:end);
    
    %% 2ii. Get the nodes in the element
    node1 = strMsh.nodes(element(1, 1), 2:end);
    node2 = strMsh.nodes(element(1, 2), 2:end);
    node3 = strMsh.nodes(element(1, 3), 2:end);
    
    %% 2iii. Create an Element Freedom Table (EFT)
    EFT = zeros(numDOFsEl, 1);
    for iEFT = 1:numNodesEl
        EFT(2*iEFT - 1) = 2*element(1, iEFT) - 1;
        EFT(2*iEFT) = 2*element(1, iEFT);
    end
    
    %% 2iv. Re-arrange the element displacement vector as [u1x u1y; u2x u2y; u3x u3y]
    uElTemp = u(EFT);
    uEl = zeros(numNodesEl, 2);
    for i = 1:numNodesEl
        uEl(i, :) = (uElTemp(2*i - 1:2*i))';
    end
    
    %% 2v. Loop over the quadrature points
    for iGP = 1:numGP
        %% 2v.1. Transform the Gauss Point location from the parameter to the physical space
        XGP = GP(iGP, 1)*node1(1, :) + GP(iGP, 2)*node2(1, :) + ...
            (1 - GP(iGP, 1) - GP(iGP, 2))*node3(1, :);
        
        %% 2v.2. Compute the basis functions and their derivatives at the Gauss Point
        [dN, areaEl, isInside] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (node1, node2, node3, XGP(1, 1), XGP(1, 2));
        if ~isInside
            error('Gauss point coordinates found outside the CST triangle');
        end
        if areaEl < minElAreaSize
            minElAreaSize = areaEl;
        end
        
        %% 2v.3. Compute the determinant of the Jacobian transformation from the physical to the parent space
        detJxxi = 2*areaEl;
        
        %% 2v.4. Compute the external load vector due to body forces
        bF = computeBodyForceVct(node1(1, 1), node1(1, 2), node1(1, 3));
        
        %% 2v.5. Compute the material, the geometric, the mass matrix, the internal residual and external body force vector at the Gauss point
        [KMaterialEl, KGeometricEl, resIntEl, FBodyEl] = ...
            computeElTangentStiffMtxResVctFEMPlateInMembraneAction ...
            (uEl, dN, bF, propParameters, C, detJxxi, GW(iGP, 1));

        %% 2v.6. Assemble the local matrices and vectors to the global ones via the EFT
        
        % Assemble to the global material stiffness matrix
        tanMtxMat(EFT, EFT) = tanMtxMat(EFT, EFT) + KMaterialEl;
        
        % Assemble to the global geometric stiffness matrix
        tanMtxGeo(EFT, EFT) = tanMtxGeo(EFT, EFT) + KGeometricEl;
        
        % Assemble to the global internal residual vector
        resVctInt(EFT) = resVctInt(EFT) + resIntEl;

        % Asseble to the global external body force vector
        Fbody(EFT) = Fbody(EFT) + FBodyEl;
    end
end

%% 3. Compute the stiffness matrix and the residual vector

% Complete stiffness matrix
tanMtx = (tanMtxMat + tanMtxGeo);

% Complete residual vector
resVct = resVctInt - loadFactor*(Fbody + F);  

end