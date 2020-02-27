function [contactLength,contactForce,maxContactPressure] = computeContactResultants...
    (mesh,displacement,lagrange,parameters)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
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
%               mesh : Elements and nodes of the mesh
%       displacement : The resulting displacement field
%           lagrange : .multipliers  : values of the Lagrange multipliers
%                    : .active_nodes : node numbers of the active nodes
%         parameters : The material properties of the structure
%
%             Output :
%      contactLength : Length of contact area
%       contactForce : The total reaction force on the contact 
% maxContactPressure : The maximal compressive pressure on the contact 
%
%% Function main body

% Initialize the array of the displaced nodes
nodesDisplaced = zeros(length(mesh.nodes),3);

% Initialize pseudocounter
n=1;

for i=1:length(mesh.nodes)
    
    % Add the x and y components of the displacement field
    nodesDisplaced(i,1) = mesh.nodes(i,1) + displacement(2*n-1);
    nodesDisplaced(i,2) = mesh.nodes(i,2) + displacement(2*n);
    
    % Update counter
    n=n+1;
end

% Compute contact force, which is just the sum of lagrange multipliers
contactForce = sum(lagrange.multipliers);

% Initialize contact length
contactLength = 0;
contactPressure = zeros(size(lagrange.active_nodes,1)-1,1);

% Loop over every lagrange multiplier
for j=2:size(lagrange.active_nodes,1)
    
    % Find x0 and y0 coordinate
    x0 = nodesDisplaced(lagrange.active_nodes(j-1),1);
    y0 = nodesDisplaced(lagrange.active_nodes(j-1),2);
    % Find x1 and y1 coordinate
    x1 = nodesDisplaced(lagrange.active_nodes(j),1);
    y1 = nodesDisplaced(lagrange.active_nodes(j),2);
    
    % Find force at node 0 and 1
    f0 = -lagrange.multipliers(j-1);
    f1 = -lagrange.multipliers(j);
    
    % Calculate the distance between the nodes of contact
    contactLength = contactLength+sqrt((x0-x1)^2+(y0-y1)^2);
    
    % Calculate the pressure in each segment
    contactPressure(j-1) = (f0+f1)/(2*parameters.t*(sqrt((x0-x1)^2+(y0-y1)^2)));
    
end

% Get the maximum pressure of all segments
maxContactPressure = max(contactPressure);

end