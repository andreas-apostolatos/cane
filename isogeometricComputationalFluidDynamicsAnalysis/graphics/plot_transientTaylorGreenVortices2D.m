function plot_transientTaylorGreenVortices2D ... 
    (p, Xi, q, Eta, CP, isNURBS, parameters, propFldDynamics, ...
    propGraph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
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

% Number of Control Points
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));

% Number of knots
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);

% Define the grid sizes
numGridPts_xi = 49;
numGridPts_eta = 49;

% Check input
checkInputForBSplineSurface ...
    (p, numKnots_xi, nxi, q, numKnots_eta, neta);

% Initialize time
t = propFldDynamics.TStart;

%% 1. Get the vertices of the rectangular domain

% Lower left corner
xi = 0;
eta = 0;
xiSpan = findKnotSpan(xi, Xi, nxi);
etaSpan = findKnotSpan(eta, Eta, neta);
R = computeIGABasisFunctionsAndDerivativesForSurface ...
    (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, 0);
x0 = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
    (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, R);

% Lower right corner
xi = 1;
eta = 0;
xiSpan = findKnotSpan(xi, Xi, nxi);
etaSpan = findKnotSpan(eta, Eta, neta);
R = computeIGABasisFunctionsAndDerivativesForSurface ...
    (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, 0);
x1 = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
    (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, R);

% Upper left corner
xi = 0;
eta = 1;
xiSpan = findKnotSpan(xi, Xi, nxi);
etaSpan = findKnotSpan(eta, Eta, neta);
R = computeIGABasisFunctionsAndDerivativesForSurface ...
    (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, 0);
x2 = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
    (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, R);

% Length step at x- and y- direction
dx = (x1 - x0)/numGridPts_xi;
dy = (x2 - x0)/numGridPts_eta;

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
        resultant = zeros(numGridPts_xi, numGridPts_eta, 2);
    else
        resultant = zeros(numGridPts_xi, numGridPts_eta);
    end
    
    % Initialize the Cartesian components
    PCartesian = zeros(numGridPts_xi, numGridPts_eta, 3);
    
    %% 2iii. Initialize the Coordinates on the Cartesian space
    x = x0;
    
    %% 2iv. Loop over all the evaluation points
    for j = 1:numGridPts_xi + 1
        for k = 1:numGridPts_eta + 1
            %% 2iv.1. Get the current Cartesian location
            PCartesian(j, k, 1:3) = x;
            
            %% 2iv.2. Compute the resultant on the evaluation point
            if strcmp(propGraph.postProcComponent, 'xVelocity')
                resultant(j, k) = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
            elseif strcmp(propGraph.postProcComponent, 'yVelocity')
                resultant(j, k) = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
            elseif strcmp(propGraph.postProcComponent, 'pressure')
                resultant(j, k) = -.25*(cos(2*x(1)) + cos(2*x(2)))*exp(-4*t*parameters.nue);
            elseif strcmp(propGraph.postProcComponent, '2normVelocity')
                uX = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
                uY = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
                resultant(j,k) = norm([uX uY]);
            elseif strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
                resultant(j, k, 1) = -cos(x(1))*sin(x(2))*exp(-2*t*parameters.nue);
                resultant(j, k, 2) = sin(x(1))*cos(x(2))*exp(-2*t*parameters.nue);
            end
            
            %% 2iv.3. Update the Cartesian loaction along the x-direction
            x = x + dx;
        end
        x = x0;
        x = x + j*dy;
    end
    
    %% 2v. Visualize the selected resultant over the domain
    
    % Initialize figure handle
    figure(propGraph.index);

    % Plot the resultant over the domain
    if strcmp(propGraph.postProcComponent,'xVelocity') || ...
            strcmp(propGraph.postProcComponent,'yVelocity') || ...
            strcmp(propGraph.postProcComponent,'pressure') || ...
            strcmp(propGraph.postProcComponent,'2normVelocity')
        surf(PCartesian(:, :, 1), PCartesian(:, :, 2), PCartesian(:, :, 3), ...
            resultant(:, :), 'EdgeColor', 'none');
    elseif strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
        quiver(PCartesian(:, :, 1), PCartesian(:, :, 2), resultant(:, :, 1), ...
            resultant(:, :, 2), 'autoscale', 'off');
    end
    
    if ~strcmp(propGraph.postProcComponent, 'velocityVectorPlot') && ...
            ~strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
        % Graphics options
        colormap('default');

        % On the color bar
        colorbar;
        hold on;
    end
    
    % Hold on the graph
    hold on;
    
    % Plot the edges of the elements
    if ~strcmp(propGraph.postProcComponent, 'velocityVectorPlot') && ...
            ~strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
        plot_knotsForBSplineSurfaceOnCartesianSpace ...
            (p, q, Xi, Eta, CP, isNURBS, 0, numGridPts_xi, numGridPts_eta);
    end

    % On the graph properties
    view(2);
    axis equal;
    xlabel('x', 'FontSize', 14);
    ylabel('y', 'FontSize', 14);
    if strcmp(propGraph.postProcComponent, 'xVelocity')
        title('Velocity component u_x');
    elseif strcmp(propGraph.postProcComponent, 'yVelocity')
        title('Velocity component u_y');
    elseif strcmp(propGraph.postProcComponent, 'pressure')
        title('Pressure distribution');
    elseif strcmp(propGraph.postProcComponent, '2normVelocity')
        title('Velocity magnitude ||u||_2');
    elseif strcmp(propGraph.postProcComponent, 'velocityVectorPlot')
        title('Contours of the velocity field u');
    end
    
    % Hold off the graph
    hold off;
    
    %% 1vi. Update the time of the simulation
    t = t + propFldDynamics.dt;
end
close(propGraph.index);

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Plotting the analytical solution took %.2d seconds \n\n', computationalTime);
    fprintf('_________Plotting Analytical Solution Ended_________\n');
    fprintf('####################################################\n');
end

end