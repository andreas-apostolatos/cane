function massMtx = computeMassMtxFEMThermalConductionAnalysis ...
    (propAnalysis, mesh, parameters, gaussInt)
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
% 2D.
%
%               Input :
%        propAnalysis : Information on the analysis type
%                           .type : Analysis type
%                mesh : The nodes and the elements of the underlying mesh
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
%   2iv. Loop over the quadrature points
%   ->
%        2iv.1. Transform the Gauss Point location from the parameter to the physical space
%
%        2iv.2. Compute the basis functions and their derivatives at the Gauss Point
%
%        2iv.3. Compute the determinant of the Jacobian from the physical to the parent space
%
%        2iv.4 Compute the element mass matrix at the Gauss Point
%
%        2iv.5. Assemble the local mass matrix global mass matrix via the EFT
%   <-
% <-
%
%% Functions main body

%% 0. Read input

% Number of nodes in the mesh
numNodes = length(mesh.nodes(:,1));

% Number of DOFs in the mesh
numDOFs = numNodes;

% Number of nodes at the element level
numNodesEl = 3;

% Number of DOFs at the element level
numDOFsEl = numNodesEl;

% Initialize output
massMtx = zeros(numDOFs,numDOFs);

%% 1. Numnerical quadrature
if strcmp(gaussInt.type, 'default')
    numGP = 1;
elseif strcmp(gaussInt.type, 'user')
    numGP = gaussInt.domainNoGP;
end
[GP, GW] = getGaussRuleOnCanonicalTriangle(numGP);

%% 2. Loop over all the elements in the mesh
for iElmnts = 1:length(mesh.elements(:, 1))
    %% 2i. Get the element in the mesh
    element = mesh.elements(iElmnts, 2:end);
    
    %% 2ii. Get the nodes in the element
    Node1 = mesh.nodes(element(1, 1), 2:end);
    Node2 = mesh.nodes(element(1, 2), 2:end);
    Node3 = mesh.nodes(element(1, 3), 2:end);
    
    %% 2iii. Create an Element Freedom Table (EFT)
    EFT = zeros(numDOFsEl, 1);
    for iEFT = 1:numNodesEl
        EFT(iEFT) = element(1, iEFT);
    end
    
    %% 2iv. Loop over the quadrature points
    for iGP = 1:numGP
        %% 2iv.1. Transform the Gauss Point location from the parameter to the physical space
        XGP = GP(iGP, 1)*Node1(1, :) + GP(iGP, 2)*Node2(1, :) + ...
            (1 - GP(iGP, 1) - GP(iGP, 2))*Node3(1, :);
        
        %% 2iv.2. Compute the basis functions and their derivatives at the Gauss Point
        [dN, Area, isInside] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (Node1, Node2, Node3, XGP(1, 1), XGP(1, 2));
        if ~isInside
            error('Gauss point coordinates found outside the CST triangle');
        end
        
        %% 2iv.3. Compute the determinant of the Jacobian from the physical to the parent space
        DetJxxi = 2*Area;
        
        %% 2iv.4. Compute the element mass matrix at the Gauss Point
        NMtx = [dN(1, 1), dN(2, 1), dN(3, 1)];
        
        %% 2iv.5. Assemble the local mass matrix global mass matrix via the EFT
        massMtx(EFT, EFT) = massMtx(EFT, EFT) + ...
            parameters.cp*parameters.rho*(NMtx'*NMtx)*DetJxxi*GW;
    end
end

end