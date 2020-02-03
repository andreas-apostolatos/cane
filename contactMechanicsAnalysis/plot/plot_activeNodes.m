function [] = plot_activeNodes...
    (mesh,displacement,lagrange)
%% Function documentation
%
% Task: Plots a red dot for every node and if a direction is specified bars
%       for the values of the Lagrange multipliers
% 
%                Input :
%                 mesh : Elements and nodes of the mesh
%         displacement : The displacement field
%         active_nodes : List with the node indices óf the active nodes
% lagrange_multipliers : Vector with Lagrange multipliers for every active node
%
%             Output :
%                  []
%
%% 
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

