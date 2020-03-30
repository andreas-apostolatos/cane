function [contactLength, contactForce, maxContactPressure] = ...
    computeContactResultants...
    (mesh, parameters, dHat, lambdaHat, nodeIDs_active)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
% Date : 27.02.2020
%
%% Function documentation
%
% Computes the length of contact, its reaction force and the maximal
% pressure. Used in unit test for the comparison with analytical values
% according to Hertz contact theory. Cannot detect an element that has all
% three nodes active.
%
%              Input :
%             strMsh : Nodes and elements in the mesh
%         parameters : Problem specific technical parameters
%               dHat : Nodal displacement field
%          lambdaHat : Nodal Lagrange Multipliers field
%     nodeIDs_active : The IDs of the active (contact) nodes
%
%             Output :
%      contactLength : Length of contact surface
%       contactForce : The total reaction force on the contact 
% maxContactPressure : The maximum compressive pressure on the contact
%
%% Function layout :
%
% 0. Read input
%
% 1. Find the mesh elements that have active nodes
% ->
%    1i. Find how many active nodes belong to each element
%
%   1ii. Check if element has 2 active nodes and store the mesh element
% <-
%
% 2. Keep only the non-zero entries of active mesh elements
%
% 3. Compute the coordinates of the displaced nodes
%
% 4. Compute contact length and contact pressure
% ->
%    4i. Get active nodes from active elements
%
%   4ii. Compute the coordinates of the first vertex of the FE edge
%
%  4iii. Compute the coordinates of the second vertex of the FE edge
%
%   4iv. Find the index of Lagrage Multipliers for each FE edge
%
%    4v. Compute the reaction force at the vertices of the segment
%
%   4vi. Compute the length of the current contact FE edge
%
%  4vii. Compute the contact pressure of the current contact FE edge
% <-
%
% 5. Get the maximum pressure from all FE edges
%
% 6. Compute the complete contact force using the nodal Lagrange Multipliers
%
%% Function main body

%% 0. Read input

% Number of total mesh elements
noMeshElements = length(mesh.elements);

% Initialize active mesh elements
activeMeshElements = zeros(noMeshElements,4);
activeElementNodes = zeros(noMeshElements,4);

% Initialize element counter
counter = 1;

%% 1. Find the mesh elements that have active nodes
for m = 1:noMeshElements
    %% 1i. Find how many active nodes belong to each element
    elementNodes = ismember(mesh.elements(m,:),nodeIDs_active);
    
    %% 1ii. Check if element has 2 active nodes and store the mesh element
    if(sum(elementNodes) == 2)
        
        % Store active elements and corresponding nodes
        activeElementNodes(counter,:) = elementNodes;
        activeMeshElements(counter,:) = mesh.elements(m,:);
        
        % Update counter
        counter = counter + 1;
    end
end
noActiveElements = counter-1;

%% 2. Keep only the non-zero entries of active mesh elements
activeMeshElements = activeMeshElements((1:noActiveElements),:);
activeElementNodes = activeElementNodes((1:noActiveElements),:);

%% 3. Compute the coordinates of the displaced nodes
nodesDisplaced = zeros(length(mesh.nodes),3);

% Loop over all mesh nodes and add displacements
for i = 1:length(mesh.nodes)
    nodesDisplaced(i,1) = mesh.nodes(i,1) + dHat(2*i - 1);
    nodesDisplaced(i,2) = mesh.nodes(i,2) + dHat(2*i);
end

%% 4. Compute contact length and contact pressure

% Initialize contact pressure and contact length
contactPressure = zeros(noActiveElements, 1);
contactLength = 0;

% Loop over all active mesh elements
for m = 1:noActiveElements
    %% 4i. Get active nodes from active elements
    
    % Get previously stored values 
    tmp_elementNodes = logical(activeElementNodes(m,:));
    tmp_meshElement = activeMeshElements(m,:);
    
    % Get only active nodes from the mesh element
    tmp_nodes = tmp_meshElement(tmp_elementNodes);
    
    %% 4ii. Compute the coordinates of the first vertex of the FE edge
    x0 = nodesDisplaced(tmp_nodes(1),1);
    y0 = nodesDisplaced(tmp_nodes(1),2);
    
    %% 4iii. Compute the coordinates of the second vertex of the FE edge
    x1 = nodesDisplaced(tmp_nodes(2),1);
    y1 = nodesDisplaced(tmp_nodes(2),2);
    
    %% 4iv. Find the index of Lagrage Multipliers for each FE edge
    lm1 = nodeIDs_active == tmp_nodes(1);
    lm2 = nodeIDs_active == tmp_nodes(2);
    
    %% 4v. Compute the reaction force at the vertices of the segment
    force0 = -lambdaHat(lm1);
    force1 = -lambdaHat(lm2);
    
    %% 4vi. Compute the length of the current contact FE edge
    contactLengthEl =  norm([x1 - x0; y1 - y0]);
    contactLength = contactLength + contactLengthEl;
    
    %% 4vii. Compute the contact pressure of the current contact FE edge
    contactPressure(m) = (force0 + force1)/2/parameters.t/contactLengthEl;
end

%% 5. Get the maximum pressure from all FE edges
maxContactPressure = max(contactPressure);

%% 6. Compute the complete contact force using the nodal Lagrange Multipliers
contactForce = sum(lambdaHat);

end