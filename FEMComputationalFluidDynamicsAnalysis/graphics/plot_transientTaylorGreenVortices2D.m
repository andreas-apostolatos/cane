function index = plot_transientTaylorGreenVortices2D ... 
    (fldMsh, parameters, propFldDynamics ,propGraph, outMsg)
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
% Plots the analytical solution of the Taylor-Green vortices over the 
% chosen time and chosen domain for a 2D incompressible flow.
%
%             Input :
%               p,q : Polynomial degrees of the 2D NURBS patch
%            Xi,Eta : The knot vectors of the NURBS patch
%                CP : The set of Control Point coordinates and weights for 
%                     the NURBS patch
%           isNURBS : Flag on whether the basis is a NURBS or a B-Spline
%        parameters : The parameters of the flow (density, viscosity)
%   propFldDynamics : Transient analysis parameters : 
%                       TStart : Start time of the simulation
%                    	  TEnd : End time of the simulation
%                  noTimeSteps : Number of time steps
%                        	dt : Time step (automatically computed)
%         propGraph : Structure containing information on the graphics
%            outMsg : Enables outputting information in the command window 
%                     when selected as 'outputEnabled'
%
% Function layout :
%
% 0. Read input
%
% 1. Get the vertices of the rectangular domain
%
% 2. Loop over all the time steps
%
%    2i. Preamble of the time stepping iterations
%
%        2iv.1. Get the current Cartesian location
%
%        2iv.2. Compute the resultant on the evaluation point
%
%        2iv.3. Update the Cartesian loaction along the x-direction
%
%    2v. Visualize the selected resultant over the domain
%
%   2vi. Update the time of the simulation
%
% 3. Appendix
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

% Initialize time
t = propFldDynamics.T0;

%% 1. Get the vertices of the rectangular domain

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

% Initialize string
reverseStr = '';

%% 2. Loop over all the time steps
for i = 1:propFldDynamics.noTimeSteps
    %% 2i. Preamble of the time stepping iterations
    
    % Print message on the progress of the load steps
    msg = sprintf('\t Time step %d/%d at real time %d seconds \n', i, ...
        propFldDynamics.noTimeSteps,t);
    fprintf([reverseStr msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    %% 2ii. Initialize the arrays
    
    % Initialize the field arrays
    if strcmp(propGraph.postProcComponent,'velocityVectorPlot')
        resultant = zeros(numGridPts_x, numGridPts_y, 2);
    else
        resultant = zeros(numGridPts_x, numGridPts_y);
    end
    
    % Initialize the Cartesian components
    PCartesian = zeros(numGridPts_x, numGridPts_y, 2);
    
    %% 2iii. Initialize the Coordinates on the Cartesian space
    x = x0;
    
    %% 2iv. Loop over all the evaluation points
    for n = 1:numGridPts_x + 1
        for m = 1:numGridPts_y + 1
            %% 2iv.1. Get the current Cartesian location
            PCartesian(n, m, 1:2) = x;
            
            %% 2iv.2. Compute the resultant on the evaluation point
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
            
            %% 2iv.3. Update the Cartesian location along the x-direction
            x = x + dx;
        end
        x = x0;
        x = x + n*dy;
    end
    
    %% 2v. Visualize the selected resultant over the domain
    
    % Initialize figure handle
    figure(propGraph.index);

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
    
    %% 1vi. Update the time of the simulation and graph index
    t = t + propFldDynamics.dt;
    propGraph.index = propGraph.index + 1;
    
end

% Assign output graph index
index = propGraph.index;
%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Plotting the analytical solution took %.2d seconds \n\n', computationalTime);
    fprintf('_________Plotting Analytical Solution Ended_________\n');
    fprintf('####################################################\n');
end

end