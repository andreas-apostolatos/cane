function index = plot_analyticalSolutionAtTimeInstance ... 
    (mesh, t, propPostproc, propGraph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documentation
%
% Plots the analytical solution specified in the function handle for 2D
% heat transfer analysis
%
%             Input :
%              mesh : Nodes and elements for the fluid mesh
%                 t : Real time at which we want to visualize the results
%      propPostproc : Structure containing information on postprocessing
%                      .computeAnalytical : Function handle to computation 
%                                           of the analtical resultant
%         propGraph : Structure containing information on the graphics
%            outMsg : Enables outputting information in the command window 
%                     when selected as 'outputEnabled'
%
%% Function layout :
%
% 0. Read input
%
% 1. Initialize the arrays
%
% 2. Loop over all the evaluation points
%
%    2i. Get the current Cartesian location
%
%    2ii. Compute the resultant on the evaluation point
%
%    2iii. Update the Cartesian loaction along the x-direction
%
% 3. Visualize the selected resultant over the domain
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('____________________________________________________\n');
    fprintf('####################################################\n');
    fprintf('\n');
    fprintf('Analytical solution of the heated plate element for \n');
    fprintf('the heat transfer problem has been initiated\n');
    fprintf('\n');
    fprintf('____________________________________________________\n');
    fprintf('\n');
    tic;
end

%% 0. Read input

% Define the number of infinite series sum
propPostproc.k = 200;

% Define the grid sizes
numGridPts_x = 49;
numGridPts_y = 49;

% Lower left corner
x0(1) = min(mesh.nodes(:,2));
x0(2) = min(mesh.nodes(:,3));

% Lower right corner
x1(1) = max(mesh.nodes(:,2));
x1(2) = min(mesh.nodes(:,3));

% Upper left corner
x2(1) = min(mesh.nodes(:,2));
x2(2) = max(mesh.nodes(:,3));

% Length step at x- and y- direction
dx = (x1 - x0)/numGridPts_x;
dy = (x2 - x0)/numGridPts_y;

% Get height and width
propPostproc.height = abs(x2(2)-x0(2));
propPostproc.width = abs(x1(1)-x0(1));

%% 1. Initialize the arrays

% Initialize the field arrays
resultant = zeros(numGridPts_x, numGridPts_y);

% Initialize the Cartesian components
PCartesian = zeros(numGridPts_x, numGridPts_y, 2);

% Initialize the Coordinates on the Cartesian space
x = x0;

%% 2. Loop over all the evaluation points
for n = 1:numGridPts_x + 1
    for m = 1:numGridPts_y + 1
        %% 2i. Get the current Cartesian location
        PCartesian(n, m, 1:2) = x;

        %% 2ii. Compute the resultant on the evaluation 
        resultant(n, m) = propPostproc.computeAnalytical(x(1),x(2),t,propPostproc);
        
        %% 2iii. Update the Cartesian location along the x-direction
        x = x + dx;
    end
    x = x0;
    x = x + n*dy;
end

%% 3. Visualize the selected resultant over the domain

% Initialize figure handle
figure(propGraph.index);
set(figure(propGraph.index),'Name','Theoretical Solution');

% Draw the color map
pcolor(PCartesian(:,:,1), PCartesian(:,:,2), resultant(:,:));

% Graph properties
hold on;
colorbar('Location', 'EastOutside');
colormap('jet');
axis tight;
shading interp;
axis equal;
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);

% Graph title
title(['Temperature distribution at t = ',num2str(t),' s']);
hold off;

% Update the graph index
index = propGraph.index + 1;

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Plotting the analytical solution took %.2d seconds \n\n', computationalTime);
    fprintf('_________Plotting Analytical Solution Ended_________\n');
    fprintf('####################################################\n');
end

end