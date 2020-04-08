function index = plot_transientTaylorGreenVortices2D ... 
    (fldMsh, parameters, t, propGraph, outMsg)
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
%   propFldDynamics : Transient analysis parameters : 
%                           T0 : Start time of the simulation
%                    	  TEnd : End time of the simulation
%                  noTimeSteps : Number of time steps
%                        	dt : Time step (automatically computed)
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
    fprintf('Analytical solution of the Taylor-Green vortices for\n');
    fprintf('the Navier-Stokes problem has been initiated\n');
    fprintf('\n');
    fprintf('____________________________________________________\n');
    fprintf('\n');
    tic;
end

%% 0. Read input

% Define the grid sizes
numGridPts_x = 49;
numGridPts_y = 49;

% Lower left corner
x0(1) = min(fldMsh.nodes(:,1));
x0(2) = min(fldMsh.nodes(:,2));

% Lower right corner
x1(1) = max(fldMsh.nodes(:,1));
x1(2) = min(fldMsh.nodes(:,2));

% Upper left corner
x2(1) = min(fldMsh.nodes(:,1));
x2(2) = max(fldMsh.nodes(:,2));

% Length step at x- and y- direction
dx = (x1 - x0)/numGridPts_x;
dy = (x2 - x0)/numGridPts_y;

%% 1. Initialize the arrays

% Initialize the field arrays
if strcmp(propGraph.postProcComponent,'velocityVectorPlot')
    resultant = zeros(numGridPts_x, numGridPts_y, 2);
else
    resultant = zeros(numGridPts_x, numGridPts_y);
end

% Initialize the Cartesian components
PCartesian = zeros(numGridPts_x, numGridPts_y, 2);

% Initialize the Coordinates on the Cartesian space
x = x0;

%% 2. Loop over all the evaluation points
for n = 1:numGridPts_x + 1
    for m = 1:numGridPts_y + 1
        %% 2i. Get the current Cartesian location
        PCartesian(n, m, 1:2) = x;

        %% 2ii. Compute the resultant on the evaluation point
        if strcmp(propGraph.postProcComponent, 'xVelocity')
            resultant(n, m) = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
        elseif strcmp(propGraph.postProcComponent, 'yVelocity')
            resultant(n, m) = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
        elseif strcmp(propGraph.postProcComponent, 'pressure')
            resultant(n, m) = -.25*(cos(2*x(1)) + cos(2*x(2)))*exp(-4*t*parameters.nue);
        elseif strcmp(propGraph.postProcComponent, '2normVelocity')
            uX = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
            uY = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
            resultant(n,m) = norm([uX uY]);
        elseif strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
            resultant(n, m, 1) = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
            resultant(n, m, 2) = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
        end

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

% Plot the resultant over the domain
if strcmp(propGraph.postProcComponent,'xVelocity') || ...
   strcmp(propGraph.postProcComponent,'yVelocity') || ...
   strcmp(propGraph.postProcComponent,'pressure')  || ...
   strcmp(propGraph.postProcComponent,'2normVelocity')

    % Draw the color map
    pcolor(PCartesian(:,:,1), PCartesian(:,:,2), resultant(:,:));

elseif strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
    quiver(PCartesian(:, :, 1), PCartesian(:, :, 2), resultant(:, :, 1), ...
        resultant(:, :, 2), 'autoscale', 'on');
end

% Define graph properties if it is not a velocity vector plot
if ~strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
    hold on;
    colorbar('Location','EastOutside');
    axis tight;
    shading interp;
end

hold on;
% Graph properties
axis equal;
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);

% Graph title
if strcmp(propGraph.postProcComponent, 'xVelocity')
    title(['Velocity component U at t = ',num2str(t)]);
elseif strcmp(propGraph.postProcComponent, 'yVelocity')
    title(['Velocity component V at t = ',num2str(t)]);
elseif strcmp(propGraph.postProcComponent, 'pressure')
    title(['Pressure distribution at t = ',num2str(t)]);
elseif strcmp(propGraph.postProcComponent, '2normVelocity')
    title(['Velocity magnitude ||u|| at t = ',num2str(t)]);
elseif strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
    title(['Contours of the velocity field u at t = ',num2str(t)]);
end
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