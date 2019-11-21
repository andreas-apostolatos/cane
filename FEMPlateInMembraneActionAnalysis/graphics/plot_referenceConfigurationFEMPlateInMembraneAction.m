function index = plot_referenceConfigurationFEMPlateInMembraneAction...
    (strMsh,analysis,F,homDBC,graph,outMsg)
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
% loads on a 2D plate in membrane action problem discretized with classical
% finite elements.
%   
%           Input :
%          strMsh : Nodes and elements in the mesh
%        analysis : .type : The analysis type
%               F : The global load vector
%          homDBC : The global numbering of the nodes where homogeneous
%                   Dirichlet boundary conditions are applied
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
% 2. Appendix
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
axis equal;
trimesh(strMsh.elements,strMsh.nodes(:,1),strMsh.nodes(:,2),strMsh.nodes(:,3),'edgecolor',edgeColor,'facecolor',faceColor);
plot_boundaryConditionsOnMesh(strMsh,homDBC,F);
grid on;
title('The initial mesh of the reference configuration');
axis on;
hold off;

% Update the graph index
index = graph.index + 1;

%% 2. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nPlotting the reference configuration took %.2d seconds \n\n',computationalTime);
    fprintf('_______________________Parsing Ended_______________________\n');
    fprintf('###########################################################\n\n\n');
end

end

