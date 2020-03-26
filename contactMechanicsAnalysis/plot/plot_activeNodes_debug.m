function [] = plot_activeNodes_debug(mesh,displacement,nodeIDs_active)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
% Date: 04.02.2020
%
%% Function documentation
%
% Plots a red dot for every node and if a direction is specified bars for
% the values of the Lagrange multipliers
% 
%            Input :
%             mesh : Elements and nodes of the mesh
%     displacement : The full displacement vector
%         lagrange : Vector with Lagrange multipliers for every active node
%
%           Output : []
%
%% Function main body

% Check if active nodes exist
if(isempty(nodeIDs_active))
    fprintf('No contact nodes !\n');
    return; 
end

% Add the x and y components of the displacement field
nodes_X = mesh.nodes(nodeIDs_active,1) + displacement(2*nodeIDs_active-1);
nodes_Y = mesh.nodes(nodeIDs_active,2) + displacement(2*nodeIDs_active);

% Plot active nodes
hold on;
scatter(nodes_X,nodes_Y,'bx','b','LineWidth',2);
hold off;

end