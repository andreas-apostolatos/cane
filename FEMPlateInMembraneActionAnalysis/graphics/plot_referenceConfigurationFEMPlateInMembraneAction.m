function index = plot_referenceConfigurationFEMPlateInMembraneAction...
    (strMsh,analysis,F,homDBC,contactSegments,graph,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the reference configuration together with the supports and the
% loads on a 2D plate in membrane action problem discretized with classic
% finite elements.
%   
%           Input :
%          strMsh : Nodes and elements in the mesh
%        analysis : .type : The analysis type
%               F : The global load vector
%          homDBC : The global numbering of the nodes where homogeneous
%                   Dirichlet boundary conditions are applied
% contactSegments : Containts the coordinates of the vertices of the 
%                   straight segments which form the rigid body's boundary:
%                       .points : Array containing the coordinates of the 
%                                 vertices of each line segment X0 - X1 in 
%                                 the form [x0, y0 ; x1 , y1]
%           graph : On the graphics
%                       .index : The index of the current graph
%          outMsg : On outputting information
%
%          Output :
%           index : The index of the current graph
%
% Function layout :
%
% 0. Read input
%
% 1. Plot the geometry together with the loads and the boundary conditions
%
% 2. Plot the contact segments defining the boundary of the rigid body
%
% 3. Assign figure properties
%
% 4. Update the graph index
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________\n');
    fprintf('###########################################################\n');
    fprintf('Plotting the reference configuration of a plate in membrane\n');
    fprintf('action has been initiated\n');
    fprintf('___________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Assign the color values
edgeColor = 'black';
faceColor = 'green';

%% 1. Plot the geometry together with the loads and the boundary conditions
plot(2);
hold on;
trimesh(strMsh.elements,strMsh.nodes(:,1),strMsh.nodes(:,2),strMsh.nodes(:,3),'edgecolor',edgeColor,'facecolor',faceColor);
plot_boundaryConditionsOnMesh(strMsh,homDBC,F);

%% 2. Plot the contact segments defining the boundary of the rigid body
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

%% 3. Assign figure properties
axis equal;
grid on;
title('The initial mesh of the reference configuration');
axis on;
hold off;

%% 4. Update the graph index
index = graph.index + 1;

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nPlotting the reference configuration took %.2d seconds \n\n',computationalTime);
    fprintf('_______________________Parsing Ended_______________________\n');
    fprintf('###########################################################\n\n\n');
end

end

