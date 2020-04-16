function index = plot_steadyStateBenchmarkProblemAnalyticalSolution ... 
    (mesh, valuesInhomDOFs, propGraph, outMsg)
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
% Plots the analytical solution of the Taylor-Green vortices for the 
% chosen time and chosen domain for a 2D incompressible flow.
%
%             Input :
%            fldMsh : Nodes and elements for the fluid mesh
%        parameters : The parameters of the flow (density, viscosity)
%                 t : Real time at which we want to visualize the results
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
numSeries = 200;

% Define the grid sizes
numGridPts_x = 49;
numGridPts_y = 49;

% Lower left corner
x0(1) = min(mesh.nodes(:,1));
x0(2) = min(mesh.nodes(:,2));

% Lower right corner
x1(1) = max(mesh.nodes(:,1));
x1(2) = min(mesh.nodes(:,2));

% Upper left corner
x2(1) = min(mesh.nodes(:,1));
x2(2) = max(mesh.nodes(:,2));

% Length step at x- and y- direction
dx = (x1 - x0)/numGridPts_x;
dy = (x2 - x0)/numGridPts_y;

% Get height and width
height = abs(x2(2)-x0(2));
width = abs(x1(1)-x0(1));

% Get temperatures
T1 = min(valuesInhomDOFs);
T2 = max(valuesInhomDOFs);

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
        coefficient = 0;
        
        % Loop over the number of infinite sum series
        for k = 1:numSeries
            coefficient = coefficient + ((( (-1)^(k+1) )+1)/k) * sin((k*pi*x(1))/width) * ...
                     ( sinh((k*pi*x(2))/width) / sinh((k*pi*height)/width) );
        end
        
        % Assign the computed temperature to the resultant
        resultant(n, m) = T1+(T2-T1)*(2/pi)*coefficient;
        
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
colorbar('Location','EastOutside');
axis tight;
shading interp;
axis equal;
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);

% Graph title
title('Temperature distribution');
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