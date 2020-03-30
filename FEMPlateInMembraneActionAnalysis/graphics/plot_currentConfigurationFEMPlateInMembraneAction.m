function index = plot_currentConfigurationFEMPlateInMembraneAction...
    (mesh, homDBC, segmentsContact, displacement, graph)
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
% segmentsContact : Containts the coordinates of the vertices of the 
%                   straight segments which form the rigid body's boundary:
%                  .numSegments : Number of rigid segments
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
% 4. Loop over all rigid contact segments to plot the boundary of the rigid body
% ->
%    4i. Plot the rigid segment
%
%   4ii. Plot a dashed segment indicating the side where the rigid body lies
% <-
%
% 5. Assign figure properties
%
% 6. Update figure index
%
%% Function main body

%% 0. Read input

% Color for the faces of the geometry
colorFace = [217 218 219]/255;

%% 1. Compute the new loactions for the vertices of the triangles in the mesh

% Initialize the array of the displaced nodes
nodes_displaced = zeros(length(mesh.nodes), 3);

% Initialize pseudocounter
counter = 1;

for i = 1:length(mesh.nodes)
    % Add the x and y components of the displacement field
    nodes_displaced(i, 1) = mesh.nodes(i, 1) + displacement(2*counter - 1);
    nodes_displaced(i, 2) = mesh.nodes(i, 2) + displacement(2*counter);
    
    % Update counter
    counter = counter + 1;
end

%% 2. Visualize the displaced elements on the mesh
figure(graph.index)

% Reference configuration
if strcmp(graph.visualization.geometry, 'reference_and_current') || ...
   strcmp(graph.visualization.geometry, 'current')
    patch('faces', mesh.elements, 'vertices', nodes_displaced(:,1:2), 'facecolor', ...
        colorFace, 'edgecolor', 'black');
    hold on;
end
% Current configuration
if strcmp(graph.visualization.geometry, 'reference_and_current') || ...
   strcmp(graph.visualization.geometry, 'reference')
    patch('faces', mesh.elements, 'vertices', mesh.nodes(:,1:2), ... 
        colorFace, 'none', 'edgecolor', 'black');
end

%% 3. Visualize the Dirichlet boundary conditions on the mesh

% Create the supports
[xs, ys, zs] = createSupports(nodes_displaced, homDBC);

% Supports
hold on;
for k = 1:length(xs(:, 1))
    plot3(xs(k, :), ys(k, :), zs(k, :), 'Linewidth', 2, 'Color', 'black');
end

%% 4. Loop over all rigid contact segments to plot the boundary of the rigid body
if ~isempty(segmentsContact)
    if isfield(segmentsContact,'points')
        if ~isempty(segmentsContact.points)
            for iSeg = 1:segmentsContact.numSegments
                %% 4i. Plot the rigid segment
                xCoord = segmentsContact.points(iSeg, [1 3]);
                yCoord = segmentsContact.points(iSeg, [2 4]);
                plot(xCoord, yCoord, '-', 'Linewidth', 2, 'Color', 'black');

                %% 4ii. Plot a dashed segment indicating the side where the rigid body lies
                segmentOffset = createSegmentOffset(segmentsContact.points(iSeg, :), ...
                    segmentsContact.normals(iSeg, :));
                plot(segmentOffset(:, 1),segmentOffset(:, 2), ':', 'LineWidth', ...
                    2, 'Color', 'black');
            end
        end
    end
end

%% 5. Assign figure properties
title('The current configuration of the mesh');
axis equal;
axis on;
grid on;
hold off;

%% 6. Update figure index
index = graph.index + 1;

end