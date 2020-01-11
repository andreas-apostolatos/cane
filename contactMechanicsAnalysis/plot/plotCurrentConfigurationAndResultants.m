%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          =============                                %
%                          || AUTHORS ||                                %
%                          =============                                %
%                                                                       %
%            Efstathios E. Theotokoglou, Professor NTUA                 %
%                   (stathis@central.ntua.gr)                           %
%                                                                       %
%            Andreas Apostolatos, Research Associate TUM                %
%                 (andreas.apostolatos@tum.de)                          %
%                                                                       %
%           ===============================================             %
%           || plot_current_configuration_and_resultants ||             %
%           ===============================================             %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = plot_current_configuration_and_resultants(mesh,rb,displacement,graph)
%% Function documentation
%
% Plots the current configuration of a plate in membrane action given the
% displacement field and visualizes the displacement/strain or stress field
% over the initial configuration of the plate
%
%        input :
%         mesh : Elements and nodes of the mesh
%           rb : Vector of the Dirichlet boundary conditions with their
%                global numbering
% displacement : The displacement field sorted in a vector according to its
%                global numbering
%        graph : Structure on the graphics (indices etc.)
%
%       output :
%        index : The index of the current graph
%
%% Function main body

%% 1. Compute the new loactions for the vertices of the triangles in the mesh

% Initialize the array of the displaced nodes
nodes_displaced = zeros(length(mesh.nodes),3);

% Initialize pseudocounter
counter = 1;

for i=1:length(mesh.nodes)
    % Add the x and y components of the displacement field
    nodes_displaced(i,1) = mesh.nodes(i,1) + displacement(2*counter-1);
    nodes_displaced(i,2) = mesh.nodes(i,2) + displacement(2*counter);
    
    % Update counter
    counter = counter + 1;
end

%% 2. Visualize the displaced elements on the mesh

figure(graph.index)

% Reference configuration
if strcmp(graph.visualization.geometry,'reference_and_current')||strcmp(graph.visualization.geometry,'current');
    patch('faces',mesh.elements,'vertices',nodes_displaced(:,1:2),'facecolor','g','edgecolor','black');
    hold on;
end
% Current configuration
if strcmp(graph.visualization.geometry,'reference_and_current')||strcmp(graph.visualization.geometry,'reference');
    patch('faces',mesh.elements,'vertices',mesh.nodes(:,1:2),'facecolor','none','edgecolor','black');
end
axis equal off;
axis on;
title('The current configuration of the mesh');

%% 3. Visualize the Dirichlet boundary conditions on the mesh

% Create the supports
[xs,ys,zs] = createSupports(nodes_displaced,rb);

%supports
hold on;
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end
hold off;

% Update index
index = graph.index + 1;

end

