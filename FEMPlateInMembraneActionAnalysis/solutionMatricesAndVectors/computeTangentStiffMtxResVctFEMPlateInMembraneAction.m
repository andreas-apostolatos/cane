function [tanMtx,resVct,minElASize] = computeTangentStiffMtxResVctFEMPlateInMembraneAction...
    (analysis,u,uSaved,uDot,uDotSaved,uMeshALE,DOFNumbering,mesh,F,loadFactor,propStrDynamics,...
    t,parameters,computeBodyForceVct,gaussInt)
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
%            analysis : Information on the analysis type
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
%        DOFNumbering : The global numbering of the DOFs
%                mesh : The nodes and the elements of the underlying mesh
%                   F : Global load vector corresponding to surface tractions
%          loadFactor : The load factor for the nonlinear steps
%     propStrDynamics : Transient analysis parameters:
%                         .method : Time integration method
%                             .T0 : Start time of the simulation
%                           .TEnd : End time of the simulation
%                             .nT : Number of time steps
%                             .dt : Time step (numeric or adaptive)
%                   t : The current time of the transient simulation
%          parameters : The parameters the physical field
% computeBodyForceVct : Function handle to body force vector computation
%            gaussInt : On the spatial integration
%                           .type : 'default', 'user'
%                     .domainNoGP : Number of Gauss Points for the domain 
%                                   integration
%                   .boundaryNoGP : Number of Gauss Points for the boundary 
%                                   integration            
%
%          Output :
%               K : The tangential stiffness matrix of the system
%          resVct : The residual vector
%      minElASize : The minimum element size in the mesh
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
nNodes = length(mesh.nodes(:,1));

% Number of DOFs in the mesh
noDOFs = 2*nNodes;

% Number of nodes at the element level
noNodesEl = 3;

% Number of DOFs at the element level
noDOFsEl = 2*noNodesEl;


% Initialize the minimum element size
firstElementInMesh = mesh.elements(1,:);
Node1 = mesh.nodes(firstElementInMesh(1,1),:);
Node2 = mesh.nodes(firstElementInMesh(1,2),:);
Node3 = mesh.nodes(firstElementInMesh(1,3),:);
minElASize = 0.5*abs((Node3(1,1)-Node2(1,1))*(Node1(1,2)-Node2(1,2))-...
    (Node1(1,1)-Node2(1,1))*(Node3(1,2)-Node2(1,2)));

% Compute the material matrix for the given problem

% Compute the material matrix for the given problem
if strcmp(analysis.type,'planeStress')
    preFactor = parameters.E/(1-parameters.nue^2);
    C = preFactor*[1             parameters.nue 0
                   parameters.nue 1              0
                   0              0             (1-parameters.nue)/2];
elseif strcmp(analysis.type,'planeStrain')
    preFactor = parameters.E*(1-parameters.nue)/(1+parameters.nue)/(1-2*parameters.nue);
    C = preFactor*[1                                 parameters.nue/(1-parameters.nue) 0
                   parameters.nue/(1-parameters.nue) 1                                 0
                   0                                 0                                 (1-2*parameters.nue)/2/(1-parameters.nue)];
else
    error('Select a valid analysis type in analysis.type');
end

% Initialize global master stiffness matrix
tanMtxMat = zeros(noDOFs,noDOFs);

% Initialize global geometric stiffness matrix
tanMtxGeo = zeros(noDOFs,noDOFs);

% Initialize the residual vector
resVctInt = zeros(noDOFs,1);

% Initialize the body force vector
Fbody = zeros(noDOFs,1);

%% 1. Numnerical quadrature
if strcmp(gaussInt.type,'default')
    noGP = 2;
elseif strcmp(gaussInt.type,'user')
    noGP = gaussInt.domainNoGP;
end
[GP,GW] = getGaussRuleOnCanonicalTriangle(noGP);

%% 2. Loop over all the elements in the mesh
for iEl = 1:length(mesh.elements(:,1))
    %% 2i. Get the element in the mesh
    element = mesh.elements(iEl,:);
    
    %% 2ii. Get the nodes in the element
    Node1 = mesh.nodes(element(1,1),:);
    Node2 = mesh.nodes(element(1,2),:);
    Node3 = mesh.nodes(element(1,3),:);
    
    %% 2iii. Create an Element Freedom Table (EFT)
    EFT = zeros(noDOFsEl,1);
    for counterEFT = 1:noNodesEl
        EFT(2*counterEFT-1) = 2*element(1,counterEFT)-1;
        EFT(2*counterEFT) = 2*element(1,counterEFT);
    end
    
    %% 2iv. Re-arrange the element displacement vector as [u1x u1y; u2x u2y; u3x u3y]
    uElTemp = u(EFT);
    uEl = zeros(noNodesEl,2);
    for i = 1:noNodesEl
        uEl(i,:) = (uElTemp(2*i-1:2*i))';
    end
    
    %% 2v. Loop over the quadrature points
    for iGP = 1:noGP
        %% 2v.1. Transform the Gauss Point location from the parameter to the physical space
        XGP = GP(iGP,1)*Node1(1,:) + GP(iGP,2)*Node2(1,:) + ...
            (1-GP(iGP,1)-GP(iGP,2))*Node3(1,:);
        
        %% 2v.2. Compute the basis functions and their derivatives at the Gauss Point
        [dN,Area,isInside] = computeCST2DBasisFunctionsAndFirstDerivatives...
            (Node1,Node2,Node3,XGP(1,1),XGP(1,2));
        if ~isInside
            error('Gauss point coordinates found outside the CST triangle');
        end
        
        %% 2v.3. Compute the determinant of the Jacobian transformation from the physical to the parent space
        DetJxxi = 2*Area;
        
        %% 2v.4. Compute the external load vector due to body forces
        bF = computeBodyForceVct(Node1(1,1),Node1(1,2),Node1(1,3));
        
        %% 2v.5. Compute the material, the geometric, the mass matrix, the internal residual and external body force vector at the Gauss point
        [KMaterialEl,KGeometricEl,resIntEl,FBodyEl] = ...
            computeElTangentStiffMtxResVctFEMPlateInMembraneAction...
            (uEl,dN,bF,parameters,C,DetJxxi,GW);

        %% 2v.6. Assemble the local matrices and vectors to the global ones via the EFT
        
        % Assemble to the global material stiffness matrix
        tanMtxMat(EFT,EFT) = tanMtxMat(EFT,EFT) + KMaterialEl;
        
        % Assemble to the global geometric stiffness matrix
        tanMtxGeo(EFT,EFT) = tanMtxGeo(EFT,EFT) + KGeometricEl;
        
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
