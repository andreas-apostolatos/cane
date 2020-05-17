function massMtx = computeMassMtxFEMPlateInMembraneAction ...
    (analysis, strMsh, propParameters, propGaussInt)
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
% problem.
%
%               Input :
%            analysis : Information on the analysis type
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
numNodes = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
numDOFs = 2*numNodes;

% Number of nodes at the element level
numNodesEl = 3;

% Number of DOFs at the element level
numDOFsEl = 2*numNodesEl;

% Initialize output
massMtx = zeros(numDOFs,numDOFs);

%% 1. Numnerical quadrature
if strcmp(propGaussInt.type, 'default')
    numGP = 2;
elseif strcmp(propGaussInt.type, 'user')
    numGP = propGaussInt.domainNoGP;
end
[GP, GW] = getGaussRuleOnCanonicalTriangle(numGP);

%% 2. Loop over all the elements in the mesh
for iElmnts = 1:length(strMsh.elements(:, 1))
    %% 2i. Get the element in the mesh
    element = strMsh.elements(iElmnts, 2:end);
    
    %% 2ii. Get the nodes in the element
    nodeI = strMsh.nodes(element(1, 1), 2:end);
    nodeJ = strMsh.nodes(element(1, 2), 2:end);
    nodeK = strMsh.nodes(element(1, 3), 2:end);
    
    %% 2iii. Create an Element Freedom Table (EFT)
    EFT = zeros(numDOFsEl, 1);
    for iEFT = 1:numNodesEl
        EFT(2*iEFT - 1) = 2*element(1, iEFT) - 1;
        EFT(2*iEFT) = 2*element(1, iEFT);
    end
    
    %% 2iv. Loop over the quadrature points
    for iGP = 1:numGP
        %% 2iv.1. Transform the Gauss Point location from the parameter to the physical space
        XGP = GP(iGP, 1)*nodeI(1, :) + GP(iGP, 2)*nodeJ(1, :) + ...
            (1 - GP(iGP, 1) - GP(iGP, 2))*nodeK(1, :);
        
        %% 2iv.2. Compute the basis functions and their derivatives at the Gauss Point
        [dN, Area, isInside] = computeCST2DBasisFunctionsAndFirstDerivatives ...
            (nodeI, nodeJ, nodeK, XGP(1, 1), XGP(1, 2));
        if ~isInside
            error('Gauss point coordinates found outside the CST triangle');
        end
        
        %% 2iv.3. Compute the determinant of the Jacobian from the physical to the parent space
        DetJxxi = 2*Area;
        
        %% 2iv.4. Compute the element mass matrix at the Gauss Point
        NMtx = [dN(1, 1) 0        dN(2, 1) 0        dN(3, 1) 0
                0        dN(1, 1) 0        dN(2, 1) 0        dN(3, 1)];
        
        %% 2iv.5. Assemble the local mass matrix global mass matrix via the EFT
        massMtx(EFT, EFT) = massMtx(EFT, EFT) + ...
            propParameters.rho*(NMtx'*NMtx)*DetJxxi*GW(iGP, 1);
    end
end

end