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
% according to Hertz contact theory. 
%
%              Input :
%             strMsh : Nodes and elements in the mesh
%         parameters : Problem specific technical parameters
%               dHat : Nodal displacement field
%          lambdaHat : Nodal Lagrange Multipliers field
%      nodeIDs_active: The IDs of the active (contact) nodes
%
%             Output :
%      contactLength : Length of contact surface
%       contactForce : The total reaction force on the contact 
% maxContactPressure : The maximum compressive pressure on the contact
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the coordinates of the displaced nodes
%
% 2. Compute the complete contact force using the nodal Lagrange Multipliers
%
% 3. Loop over all contact finite element edges
% ->
%    3i. Compute the coordinates of the first vertex of the FE edge
%
%   3ii. Compute the coordinates of the second vertex of the FE edge
%
%  3iii. Compute the reaction force at the vertices of the segment
%
%   3iv. Compute the length of the current contact FE edge
%
%    3v. Compute the contact pressure of the current contact FE edge
% <-
%
% 4. Get the maximum pressure from all FE edges
%
%% Function main body
warning('This function assumes that the contact nodes are sorted in a directional manner with respect to the boundary of the deformable body');

%% 0. Read input

% Initialize contact length
contactLength = 0;

% Number of contact finite element edges
noContactFEEdges = size(nodeIDs_active,1) - 1;

% Initialize contact pressure which is computed element-wise
contactPressure = zeros(noContactFEEdges, 1);

%% 1. Compute the coordinates of the displaced nodes
nodesDisplaced = zeros(length(mesh.nodes),3);
for i = 1:length(mesh.nodes)
    nodesDisplaced(i,1) = mesh.nodes(i,1) + dHat(2*i - 1);
    nodesDisplaced(i,2) = mesh.nodes(i,2) + dHat(2*i);
end

%% 2. Compute the complete contact force using the nodal Lagrange Multipliers
contactForce = sum(lambdaHat);

%% 3. Loop over all contact finite element edges
for iCE = 1:noContactFEEdges
    %% 3i. Compute the coordinates of the first vertex of the FE edge
    x0 = nodesDisplaced(nodeIDs_active(iCE, 1),1);
    y0 = nodesDisplaced(nodeIDs_active(iCE, 1),2);
    
    %% 3ii. Compute the coordinates of the second vertex of the FE edge
    x1 = nodesDisplaced(nodeIDs_active(iCE + 1, 1),1);
    y1 = nodesDisplaced(nodeIDs_active(iCE + 1, 1),2);
    
    %% 3iii. Compute the reaction force at the vertices of the segment
    force0 = - lambdaHat(iCE, 1);
    force1 = - lambdaHat(iCE + 1, 1);
    
    %% 3iv. Compute the length of the current contact FE edge
    contactLengthEl =  norm([x1 - x0; y1 - y0]);
    contactLength = contactLength + contactLengthEl;
    
    %% 3v. Compute the contact pressure of the current contact FE edge
    contactPressure(iCE, 1) = (force0 + force1)/2/parameters.t/contactLengthEl;
end

%% 4. Get the maximum pressure from all FE edges
maxContactPressure = max(contactPressure);

end