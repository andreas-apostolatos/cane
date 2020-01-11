function K_element = computeElementStiffnessMatrixPlateInMembraneActionLinearCST...
    (nodes,materialProperties,analysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the element stiffness matrix corresponding to the CST (Constant-
% Strain Triangle) element. Due to the linearity of the basis functions no
% numerical integration is needed. 
%
% Closed form expression for the computation of the element stiffness 
% matrix exists:
% G.I. Tsamasfyros and E.E. Theotokoglou, The finite element method I and II, Symmetria, Athens (2000) 
%
%              Input :
%              nodes : The nodes of the CST in a counterclock-wise fashion
% materialProperties : The material properties of the structure
%           analysis : Analysis type (plane stress or plane strain)
%
%             Output :
%          K_element : The element stiffness matrix for the CST element
%
% Function layout:
%
% 0. Read input
%
% 1. Compute the element stiffness matrix
%
% 2. Check the rigid body modes of the element stiffness matrix
%
%% Function main body 

%% 0. Read input

% Get the nodes of the mesh in a counter-clockwise fashion
%
%        i ------------- k
%          \           /
%           \         /
%            \       /
%             \     /  
%              \   /
%               \ /
%                j
vertix_i = nodes(1,:);
vertix_j = nodes(2,:);
vertix_k = nodes(3,:);

% Compute the permutation symbols
% For basis function Ni:
% yjk:
yjk = vertix_j(2)-vertix_k(2);
% xkj:
xkj = vertix_k(1) - vertix_j(1);

% For basis function Nj:
% yik:
yik = -(vertix_i(2)-vertix_k(2));
% xki:
xki = -(vertix_k(1) - vertix_i(1));

% For basis function Nj:
% yij:
yij = vertix_i(2)-vertix_j(2);
% xji:
xji = vertix_j(1) - vertix_i(1);

% The area of the triangle
% Delta = det(nodes)/2;
Delta = abs((nodes(1,1)-nodes(3,1))*(nodes(2,2)-nodes(1,2))-(nodes(1,1)-nodes(2,1))*(nodes(3,2)-nodes(1,2)))/2;

%% 1. Compute the element stiffness matrix

% Compute the B-operator matrix
preFactor_B = 1/2/Delta;
Boperator = preFactor_B*[yjk 0 yik 0 yij 0
                         0 xkj 0 xki 0 xji
                         xkj yjk xki yik xji yij];
                     
% Compute the material matrix                        
if strcmp(analysis.type,'plainStress')
    preFactor_material = materialProperties.E/(1-materialProperties.nue^2);
    Dm = preFactor_material*[1 materialProperties.nue 0;
                             materialProperties.nue 1 0
                             0 0 (1-materialProperties.nue)/2];
elseif strcmp(analysis.type,'plainStrain')
    preFactor_material = materialProperties.E*(1-materialProperties.nue)/(1+materialProperties.nue)/(1-2*materialProperties.nue);
    Dm = preFactor_material*[1 materialProperties.nue/(1-materialProperties.nue) 0;
                             materialProperties.nue/(1-materialProperties.nue) 1 0
                             0 0 (1-2*materialProperties.nue)/2/(1-materialProperties.nue)];
end

% Compute the element stiffness matrix as a tensor product (B'*D*B*t*Delta see page 46 of the book)
K_element = Boperator'*Dm*Boperator*Delta*materialProperties.t;

%% 2. Check the rigid body modes of the element stiffness matrix

% Compute rank and size of the element stiffness matrix
[sizeK,~] = size(K_element);
rankK = rank(K_element);

rigid_body_modes = sizeK - rankK;

% For plane analysis the rigid body modes must be 3: 2 translatoric and one
% rotational modes. Otherwise print warning message
if rigid_body_modes~=3
   fprintf('Warning: Element stiffness matrix has %d rigid body modes',rigid_body_modes); 
end

end

