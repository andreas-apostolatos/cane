function [] = plot_lagrangeMultipliers...
    (mesh,displacement,lagrange,direction)
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
%            direction : Can take the values empty '' for no bar, 'h' for horizontal 
%                        drawing of the bars and 'v' for vertical drawing
%
%             Output :
%                  []
%
%% Check if active nodes exist

% assign variables
active_nodes = lagrange.active_nodes;
lagrange_multipliers = lagrange.multipliers;

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

%% Add the Lagrange multipliers bars

if(direction=='h')

    h1=gca;
    h2 = axes('Position',get(h1,'Position'));
    barh(h2,mesh.nodes(active_nodes,2)+displacement(2*active_nodes),abs(lagrange_multipliers),0.1,'EdgeColor',[1 0.1 0.5]);
    set(h2,'XAxisLocation','top','Color','none','YTickLabel',[]);
    set(h2,'Box','off');
    set(h2, 'xdir','reverse');
    z=zoom;
    setAllowAxesZoom(z,h2,false); 
    linkaxes([h1 h2],'y');
elseif (direction=='v')
    h1=gca;
    h2 = axes('Position',get(h1,'Position'));
    bar(h2,mesh.nodes(active_nodes,1)+displacement(2*active_nodes-1),abs(lagrange_multipliers),0.02,'EdgeColor',[1 0.1 0.5]);
    set(h2,'YAxisLocation','right','Color','none','XTickLabel',[]);
    set(h2,'XLim',get(h1,'XLim'),'Layer','bot');
    set(h2,'Box','off');
    z=zoom;
    setAllowAxesZoom(z,h2,false); 
    linkaxes([h1 h2],'x');  
end

end

