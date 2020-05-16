function K = computeStiffnessMatrixPlateInMembraneActionLinear ...
    (mesh, propParameters, propAnalysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Computes the stiffness matrix needed for linear analysis of a plate in
% membrane action for constant strain triangle (CST)
%
%              Input :
%               mesh : The mesh of the structure
%     propParameters : The material properties of the structure
%       propAnalysis : Analysis type (plane stress or plane strain)
%
%             Output :
%                  K : Master stiffness matrix
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all elements
%
% 2. Check the rigid body modes of the master stiffness matrix
%
%% Function main body

%% 0. Read input

% Number of nodes for the global system
no_nodes_global = length(mesh.nodes);

% Number of degrees of freedom for the global system
% for 2D displacement based analysis
no_dofs_global = 2*no_nodes_global;

% Number of DoFs at the element level (depends on the element type)
no_nodes_element = 3;
no_dofs_element = no_nodes_element*2;

% Initialize master stiffness matrix
K = zeros(no_dofs_global,no_dofs_global);

%% 1. Loop over all elements
for i=1:length(mesh.elements)
    
    % Get the current element in the mesh
    element = mesh.elements(i,2:no_nodes_element+1);
    
    % Get the nodes of the triangle in a counterclockwise fashion
    nodes = mesh.nodes(element,2:end);
    
    % Compute element stiffness matrix for the CST
    K_element = computeElementStiffnessMatrixPlateInMembraneActionLinearCST...
                                       (nodes,propParameters,propAnalysis);
    
    % Assemble to the global stiffness matrix via element freedom tables
    % Element freedom table
    EFT = zeros(1,no_dofs_element);
    
    % Assign the entried of the element freedom table recursively
    for j=1:no_nodes_element
        EFT(1,2*j-1) = 2*element(j)-1;
        EFT(1,2*j) = 2*element(j);
    end
    
    K(EFT,EFT) = K(EFT,EFT) + K_element;  
end

%% 2. Check the rigid body modes of the master stiffness matrix

% Compute rank and size of the master stiffness matrix
[sizeK,~] = size(K);
rankK = rank(K);

rigid_body_modes = sizeK - rankK;

% For plane analysis the rigid body modes must be 3: 2 translatoric and one
% rotational modes. Otherwise print warning message
if rigid_body_modes~=3
   fprintf('Warning: Master stiffness matrix has %d rigid body modes',rigid_body_modes); 
end

end

