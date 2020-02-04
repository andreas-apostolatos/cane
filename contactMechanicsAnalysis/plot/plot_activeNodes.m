function [] = plot_activeNodes(mesh,displacement,lagrange)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
% Date : 04.02.2020
%
%% Function documentation
%
% Task: Plots a red dot for every node and if a direction is specified bars
%       for the values of the Lagrange multipliers
% 
%            Input :
%             mesh : Elements and nodes of the mesh
%     displacement : The full displacement field
%         lagrange : Vector with Lagrange multipliers for every active node
%
%           Output : []
%
%% Function main body

% assign variables
active_nodes = lagrange.active_nodes;

% Check if active nodes exist
if(isempty(active_nodes))
    fprintf('No contact nodes !\n');
    return; 
end

%% Plot active nodes

% Add the x and y components of the displacement field
nodes_displacedX = mesh.nodes(active_nodes,1) + displacement(2*active_nodes-1);
nodes_displacedY = mesh.nodes(active_nodes,2) + displacement(2*active_nodes);
hold on;
scatter(nodes_displacedX,nodes_displacedY,'ro','b','LineWidth',3);
hold off;

end

