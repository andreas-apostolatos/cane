function [K, F, minElEdgeSize] = ...
    computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed...
    (propAnalysis, u, uSaved, uDot, uDotSaved, uMeshALE, precompStiffMtx, ...
    precomResVct, DOFNumbering, strMsh, F, loadFactor, computeBodyForces, ...
    propStrDynamics, t, propParameters, propInt)
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
% plate in membrane action analysis using both linear triangles and
% quadrilaterals for the discretization of the meduim.
%
%             Input :
%      propAnalysis : Structure containing general information on the 
%                     analysis,
%                       .type : Analysis type
%                 u : The solution of the previous iteration
%            uSaved : The discrete solution field of the previous 
%                     time step
%              uDot : The time derivative of the solution of the previous 
%                     iteration
%         uDotSaved : The time derivative of the discrete solution 
%                     field of the previous time step
%          uMeshALE : Dummy variable for this function
%   precompStiffMtx : Constant part of the stiffness matrix which can be
%                     precomputed
%      precomResVct : Constant part of the residual vector which can be
%                     precomputed
%      DOFNumbering : Dummy variable for this function
%            strMsh : Nodes and elements in the mesh
%                 F : The global load vector corresponding to surface
%                     tractions
%        loadFactor : Load factor, nessecary for nonlinear computations 
%                     (dummy variable for this function)
% computeBodyForces : Function handle to body force vector computation
%   propStrDynamics : Dummy variable for this function
%                 t : Time instance
%    propParameters : Problem specific technical parameters
%           propInt : On the numerical quadrature
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
% 1. Choose the numnerical quadrature rule
%
% 2. Loop over all the elements in the mesh
% ->
%    2i. Get the element in the mesh and find the number of its nodes
%
%   2ii. Get the nodes in the element
%
%  2iii. Create an Element Freedome Table (EFT)
%
%   2iv. Loop over all the quadrature points
%   ->
%        2iv.1. Transform the Gauss Point location from the parameter to the physical space if the element is a constant triangle
%
%        2iv.2. Compute the basis functions and their derivatives at the Gauss Point
%
%        2iv.3. Compute the Jacobian of the transformation from the parameter to the physical space and the physical coordinates of the Gauss Point if the element is a bilinear quadrilateral
%
%        2iv.4. Form the basis functions matrix at the Gauss Point
%
%        2iv.5. Transform the derivatives from the physical to the parameter if the element is a bilinear quadrilateral
%
%        2iv.6. Form the B-Operator matrix for the plate in membrane action problem
%
%        2iv.7. Compute the element stiffness matrix at the Gauss Point and assemble to master stiffness matrix via the EFT
%
%        2iv.8. Compute the element load vector due to body forces and assemble to master load vector via the EFT
%   <-
%
%    2v. Check for minimum element edge size
% <-
%
%% Functions main body

%% 0. Read input

% Number of nodes in the mesh
noNodes = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
noDOFs = 2*noNodes;

% Initialize the minimum element edge size
firstElementInMesh = strMsh.elements(1,:);
firstElementInMesh(isnan(firstElementInMesh)) = [];
noElNodes = length(firstElementInMesh);
if noElNodes == 3
    elementType = 'linearTriangle';
elseif noElNodes == 4
    elementType = 'bilinearQuadrilateral';
else
    error('The computation of the element stiffness matrix for a %d-noded element has not been implemented',noElNodes);
end
nodes = strMsh.nodes(firstElementInMesh,:);
if strcmp(elementType,'linearTriangle')
    minElEdgeSize = min([norm(nodes(1,:)-nodes(2,:)) norm(nodes(1,:)-nodes(3,:)) norm(nodes(2,:)-nodes(3,:))]);
elseif strcmp(elementType,'bilinearQuadrilateral')
    minElEdgeSize = min([norm(nodes(1,:)-nodes(2,:)) norm(nodes(2,:)-nodes(3,:)) norm(nodes(3,:)-nodes(4,:)) norm(nodes(4,:)-nodes(1,:))]);
end

% Compute the material matrix for the given problem
if strcmp(propAnalysis.type,'planeStress')
    preFactor = propParameters.E/(1-propParameters.nue^2);
    C = preFactor*[1             propParameters.nue 0
                   propParameters.nue 1              0
                   0              0             (1-propParameters.nue)/2];
elseif strcmp(propAnalysis.type,'planeStrain')
    preFactor = propParameters.E*(1-propParameters.nue)/(1+propParameters.nue)/(1-2*propParameters.nue);
    C = preFactor*[1                                 propParameters.nue/(1-propParameters.nue) 0
                   propParameters.nue/(1-propParameters.nue) 1                                 0
                   0                                 0                                 (1-2*propParameters.nue)/2/(1-propParameters.nue)];
end

% Initialize output array
K = zeros(noDOFs,noDOFs);

%% 1. Choose the numnerical quadrature rule

% Choose the integration rule for a triangle element
if strcmp(propInt.type,'default')
    noGPTriangle = 1;
elseif strcmp(propInt.type,'user')
    noGPTriangle = propInt.noGP;
end
[GPTriangle,GWTriangle] = getGaussRuleOnCanonicalTriangle(noGPTriangle);

% Choose the integration rule for a quadrilateral element
if strcmp(propInt.type,'default')
    noGPQuadrilateral = 2;
elseif strcmp(propInt.type,'user')
    noGPQuadrilateral = propInt.nGP;
end
[GPQuad,GWQuad] = getGaussPointsAndWeightsOverUnitDomain...
    (noGPQuadrilateral);
GPQuadrilateral = zeros(length(GWQuad)^2,2);
GWQuadrilateral = zeros(1,length(GWQuad)^2);
counterGPuv = 1;
for counterGPu = 1:length(GWQuad)
    for counterGPv = 1:length(GWQuad)
        GPQuadrilateral(counterGPuv,:) = [GPQuad(1,counterGPu) GPQuad(1,counterGPv)];
        GWQuadrilateral(1,counterGPuv) = GWQuad(1,counterGPu)*GWQuad(1,counterGPv);
        counterGPuv = counterGPuv + 1;
    end
end

%% 2. Loop over all the elements in the mesh
for counterEl = 1:length(strMsh.elements(:,1))
    %% 2i. Get the element in the mesh and find the number of its nodes
    element = strMsh.elements(counterEl,:);
    element(isnan(element)) = [];
    noElNodes = length(element);
    if noElNodes == 3
        elementType = 'linearTriangle';
    elseif noElNodes == 4
        elementType = 'bilinearQuadrilateral';
    else
        error('The computation of the element stiffness matrix for a %d-noded element has not been implemented',noElNodes);
    end
    noElDOFs = 2*noElNodes;
    if strcmp(elementType,'linearTriangle')
        noGP = noGPTriangle;
        GW = GWTriangle;
    elseif strcmp(elementType,'bilinearQuadrilateral')
        noGP = noGPQuadrilateral^2;
        GW = GWQuadrilateral;
    end
    
    %% 2ii. Get the nodes in the element
    nodes = zeros(noElNodes,3);
    for counterNodesEl = 1:noElNodes
        nodes(counterNodesEl,:) = strMsh.nodes(element(1,counterNodesEl),:);
    end
    
    %% 2iii. Create an Element Freedome Table (EFT)
    EFT = zeros(noElDOFs,1);
    for counterEFT = 1:noElNodes
        EFT(2*counterEFT-1) = 2*element(1,counterEFT)-1;
        EFT(2*counterEFT) = 2*element(1,counterEFT);
    end
    
    %% 2iv. Loop over all the quadrature points
    for counterGP = 1:noGP
        %% 2iv.1. Transform the Gauss Point location from the parameter to the physical space if the element is a constant triangle
        if strcmp(elementType,'linearTriangle')
            xGP = GPTriangle(counterGP,1)*nodes(1,:) + GPTriangle(counterGP,2)*nodes(2,:) + ...
                (1-GPTriangle(counterGP,1)-GPTriangle(counterGP,2))*nodes(3,:);
        end
        
        %% 2iv.2. Compute the basis functions and their derivatives at the Gauss Point
        if strcmp(elementType,'linearTriangle')
            [dN,detJxix] = computeCST2DBasisFunctionsAndFirstDerivatives...
                (nodes(1,:),nodes(2,:),nodes(3,:),xGP(1,1),xGP(1,2));
        elseif strcmp(elementType,'bilinearQuadrilateral')
            dN = computeBilinearBasisFunctionsAndFirstDerivatives...
                (GPQuadrilateral(counterGP,1),GPQuadrilateral(counterGP,2));
        end
        
        %% 2iv.3. Compute the Jacobian of the transformation from the parameter to the physical space and the physical coordinates of the Gauss Point if the element is a bilinear quadrilateral
        if strcmp(elementType,'bilinearQuadrilateral')
            Jxix = zeros(2,2);
            xGP = zeros(1,3);
            for counterJacobian = 1:noElNodes
                Jxix(1,1) = Jxix(1,1) + dN(counterJacobian,2)*nodes(counterJacobian,1);
                Jxix(1,2) = Jxix(1,2) + dN(counterJacobian,2)*nodes(counterJacobian,2);
                Jxix(2,1) = Jxix(2,1) + dN(counterJacobian,3)*nodes(counterJacobian,1);
                Jxix(2,2) = Jxix(2,2) + dN(counterJacobian,3)*nodes(counterJacobian,2);
                xGP = xGP + dN(counterJacobian,1)*nodes(counterJacobian,:);
            end
            detJxix = det(Jxix);
        end
        
        %% 2iv.4. Form the basis functions matrix at the Gauss Point
        NMatrix = zeros(2,noElDOFs);
        for counterElNodes = 1:noElNodes
            NMatrix(1,2*counterElNodes-1) = dN(counterElNodes,1);
            NMatrix(2,2*counterElNodes) = dN(counterElNodes,1);
        end

        %% 2iv.5. Transform the derivatives from the physical to the parameter if the element is a bilinear quadrilateral
        if strcmp(elementType,'linearTriangle')
            dNdx = dN(:,2:3);
        elseif strcmp(elementType,'bilinearQuadrilateral')
            dNdx = zeros(noElNodes,2);
            for counterElNodes = 1:noElNodes
                dNdx(counterElNodes,:) = Jxix\dN(counterElNodes,2:3)';
            end
        end

        %% 2iv.6. Form the B-Operator matrix for the plate in membrane action problem
        B = zeros(3,noElDOFs);
        for counterElNodes = 1:noElNodes
           B(1,2*counterElNodes-1) = dNdx(counterElNodes,1);
           B(2,2*counterElNodes) = dNdx(counterElNodes,2);
           B(3,2*counterElNodes-1) = dNdx(counterElNodes,2);
           B(3,2*counterElNodes) = dNdx(counterElNodes,1);
        end
         
        %% 2iv.7. Compute the element stiffness matrix at the Gauss Point and assemble to master stiffness matrix via the EFT
        K(EFT,EFT) = K(EFT,EFT) + (B'*C*B)*abs(detJxix)*GW(counterGP);
        
        %% 2iv.8. Compute the element load vector due to body forces and assemble to master load vector via the EFT
%         bF = bodyForces(xGP(1,1),xGP(1,2),xGP(1,3));
%         F(EFT) = F(EFT) + NMatrix(:,1)'*bF(1:2,1)*abs(detJxix)*GW(counterGP);
    end
    
    %% 2v. Check for minimum element edge size
    if strcmp(elementType,'linearTriangle')
        h = min([norm(nodes(1,:)-nodes(2,:)) norm(nodes(1,:)-nodes(3,:)) norm(nodes(2,:)-nodes(3,:))]);
    elseif strcmp(elementType,'bilinearQuadrilateral')
        h = min([norm(nodes(1,:)-nodes(2,:)) norm(nodes(2,:)-nodes(3,:)) norm(nodes(3,:)-nodes(4,:)) norm(nodes(4,:)-nodes(1,:))]);
    end
    if h < minElEdgeSize
        minElEdgeSize = h;
    end
end

end
