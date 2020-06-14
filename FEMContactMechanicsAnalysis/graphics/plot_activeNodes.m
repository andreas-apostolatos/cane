function plot_activeNodes(mesh, dHat, nodeIDs_active)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% Plots a red dot for every node and if a direction is specified bars for
% the values of the Lagrange multipliers
% 
%            Input :
%             mesh : Elements and nodes of the mesh
%             dHat : The full displacement vector
%   nodeIDs_active : Nodal IDs of the nodes which are found to come in
%                    contact with any of the rigid segments
%
%           Output :
%                    graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Add the displacement field to the nodal coordinates
%
% 2. Plot active nodes
%
%% Function main body

%% 0. Read input 

% Check if active nodes exist
if(isempty(nodeIDs_active))
    fprintf('No contact nodes !\n');
    return; 
end

%% 1. Add the displacement field to the nodal coordinates
nodes_X = mesh.nodes(nodeIDs_active, 2) + dHat(2*nodeIDs_active - 1);
nodes_Y = mesh.nodes(nodeIDs_active, 3) + dHat(2*nodeIDs_active);

%% 2. Plot active nodes
scatter(nodes_X,nodes_Y,'magentao','b','LineWidth',3);

end
