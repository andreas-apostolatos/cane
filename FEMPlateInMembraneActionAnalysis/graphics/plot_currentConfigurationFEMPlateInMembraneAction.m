function index = plot_currentConfigurationFEMPlateInMembraneAction...
    (mesh,homDBC,contactSegments,displacement,graph)
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the current configuration of a plate in membrane action given the
% resulting displacement field 
%
%           input :
%            mesh : Elements and nodes of the mesh
%          homDBC : Vector of the Dirichlet boundary conditions with their
%                   global numbering
% contactSegments : Containts the coordinates of the vertices of the 
%                   straight segments which form the rigid body's boundary:
%                       .points : Array containing the coordinates of the 
%                                 vertices of each line segment X0 - X1 in 
%                                 the form [x0, y0 ; x1 , y1]
%    displacement : The displacement field sorted in a vector according to 
%                   its global numbering
%           graph : On the graphics
%                       .index : The index of the current figure
%
%       output :
%        index : The index of the current figure
%
% Function layout :
%
% 1. Compute the new loactions for the vertices of the triangles in the mesh
%
% 2. Visualize the displaced elements on the mesh
%
% 3. Visualize the Dirichlet boundary conditions on the mesh
%
% 4. Plot the contact segments defining the boundary of the rigid body
%
% 5. Assign figure properties
%
% 6. Update figure index
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
if strcmp(graph.visualization.geometry,'reference_and_current')||...
   strcmp(graph.visualization.geometry,'current')
    patch('faces',mesh.elements,'vertices',nodes_displaced(:,1:2),'facecolor','g','edgecolor','black');
    hold on;
end
% Current configuration
if strcmp(graph.visualization.geometry,'reference_and_current')||...
   strcmp(graph.visualization.geometry,'reference')
    patch('faces',mesh.elements,'vertices',mesh.nodes(:,1:2),'facecolor','none','edgecolor','black');
end

%% 3. Visualize the Dirichlet boundary conditions on the mesh

% Create the supports
[xs,ys,zs] = createSupports(nodes_displaced,homDBC);

% Supports
hold on;
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end

%% 4. Plot the contact segments defining the boundary of the rigid body
if ~isempty(contactSegments)
    if isfield(contactSegments,'points')
        if ~isempty(contactSegments.points)
            for iPt = 1:size(contactSegments.points,3)
                % plot a line form start to end coordiante x0 to x1
                startPoint = contactSegments.points(:,1,iPt);
                endPoint = contactSegments.points(:,2,iPt);
                plot(startPoint,endPoint,'-','Linewidth',2,'Color','black');

                % create contour
                contourWall = createWall(contactSegments.points(:,:,iPt));

                % plot contours below the line to indicate orientation
                plot(contourWall(:,1),contourWall(:,2),':','LineWidth',2,'Color','black');
            end
        end
    end
end

%% 5. Assign figure properties
axis equal;
axis on;
grid on;
title('The current configuration of the mesh');
hold off;

%% 6. Update figure index
index = graph.index + 1;

end