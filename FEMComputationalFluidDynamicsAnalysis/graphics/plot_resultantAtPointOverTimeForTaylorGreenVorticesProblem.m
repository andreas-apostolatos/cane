function index = plot_resultantAtPointOverTimeForTaylorGreenVorticesProblem ...
    (x, y, fldMsh, parameters, upHistory, propFldDynamics, propGraph)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documenation
%
% Returns a 2D plot with the evolution of a resultant throught the
% transient analysis at a specific point for 2D incompressible flow problem
% using FEM analysis for the benchmark of Taylor-Green votices together with
% the analytical resultant at the chosen point.
%
%                     Input :
%                       x,y : Cartesian coordinates of evaluation point
%                    fldMsh : Nodes and elements for the fluid mesh
%                parameters : The parameters of the flow (density, 
%                             viscosity)
%                 upHistory : The history data of the transient analysis by
%                             solving the problem numerically
%           propFldDynamics : Transient analysis parameters : 
%                                  T0 : Start time of the simulation
%                                TEnd : End time of the simulation
%                         noTimeSteps : Number of time steps
%                                  dt : Time step (automatically computed)
%                 propGraph : Structure containing information on the 
%                             graphics and postprocessing values,
%                               .postProcComponent : Postprocessing 
%                                                    component to
%                                                    visualize :
%                                                    'xVelocity'
%                                                    'yVelocity'
%                                                    'pressure'
%                                                    '2normVelocity'
%
%                   Output :  graphics
%
%% Function layout :
%
% 0. Read input
%
% 1. Find the mesh element from cartesian coordinates
%
% 2. Compute the DOFs that corespond to element nodes
%
% 3. Loop over all the time steps
%
%    3i. Fill up the time step array accordingly
%
%   3ii. Get the current discrete solution vector
%
%  3iii. Get the actual resultants at the point and at the current time step
%
%   3iv. Assign the resultants from numerical and analytical solution
%
%    3v. Update the simulation time
%
% 4. Plot the resultant at the given point and over time
%
% 5. Update the figure handle index
%
%% Function main body

%% 0. Read input

% Initialize figure handle
figure(propGraph.index)

% Initialize the output arrays
resultantAtPointNumerical = zeros(propFldDynamics.noTimeSteps, 1);
resultantAtPointAnalytical = zeros(propFldDynamics.noTimeSteps, 1);
timeSteps = zeros(propFldDynamics.noTimeSteps, 1);

% Initialize time
t = propFldDynamics.dt;

%% 1. Find the mesh element from cartesian coordinates
for iElement = 1:size(fldMsh.elements,1)

    % Find the node IDs of an element
    node_IDs = fldMsh.elements(iElement,:);
    
    % Find the node coordinates of an element
    vertexI = fldMsh.nodes(node_IDs(1),1:2);
    vertexJ = fldMsh.nodes(node_IDs(2),1:2);
    vertexK = fldMsh.nodes(node_IDs(3),1:2);
    
    % Compute basis functions and check if the point is inside the element
    [N,~,isInside] = computeCST2DBasisFunctions(vertexI,vertexJ,vertexK,x,y);
    
    % Break the loop if the point is inside the element
    if isInside
       break
    end
end

%% 2. Compute the DOFs that corespond to element nodes
vertexI_DOFs = (node_IDs(1)*3-2):(node_IDs(1)*3);
vertexJ_DOFs = (node_IDs(2)*3-2):(node_IDs(2)*3);
vertexK_DOFs = (node_IDs(3)*3-2):(node_IDs(3)*3);

%% 3. Loop over all the time steps
for iTime = 1:propFldDynamics.noTimeSteps
    %% 3i. Fill up the time step array accordingly
    timeSteps(iTime) = t;
    
    %% 3ii. Get the current discrete solution vector
    upCurrent = upHistory(:,iTime+1);
    
    %% 3iii. Get the actual resultants at the point and at the current time step
    upVector = upCurrent(vertexI_DOFs)*N(1) + ...
               upCurrent(vertexJ_DOFs)*N(2) + ...
               upCurrent(vertexK_DOFs)*N(3);
    
    %% 3iv. Assign the resultants from numerical and analytical solution
    if strcmp(propGraph.postProcComponent, 'xVelocity')
        resultantAtPointNumerical(iTime, 1) = upVector(1);
        resultantAtPointAnalytical(iTime, 1) = -cos(x)*sin(y)*exp(-2*t*parameters.nue);
    elseif strcmp(propGraph.postProcComponent, 'yVelocity')
        resultantAtPointNumerical(iTime, 1) = upVector(2);
        resultantAtPointAnalytical(iTime, 1) = sin(x)*cos(y)*exp(-2*t*parameters.nue);
    elseif strcmp(propGraph.postProcComponent, 'pressure')
        resultantAtPointNumerical(iTime, 1) = upVector(3);
        resultantAtPointAnalytical(iTime, 1) = -.25*(cos(2*x) + cos(2*y))*exp(-4*t*parameters.nue);
    elseif strcmp(propGraph.postProcComponent, '2normVelocity')
        uXNumerical = upVector(1);
        uYNumerical = upVector(2);
        resultantAtPointNumerical(iTime, 1) = norm([uXNumerical uYNumerical]);
        uXAnalytical = -cos(x)*sin(y)*exp(-2*t*parameters.nue);
        uYAnalytical = sin(x)*cos(y)*exp(-2*t*parameters.nue);
        resultantAtPointAnalytical(iTime, 1) = norm([uXAnalytical uYAnalytical]);
    end
   
    %% 3v. Update the simulation time
    t = t + propFldDynamics.dt;
end

%% 4. Plot the resultant at the given point and over time
plot(timeSteps, resultantAtPointAnalytical, 'b',...
     timeSteps, resultantAtPointNumerical, 'r');
legend('Analytical', 'Navier-Stokes','Orientation','horizontal','Location','southoutside');
xlabel('time (seconds)');
if strcmp(propGraph.postProcComponent, 'xVelocity')
    yLabelString = 'x-velocity component u_x (m/s)';
elseif strcmp(propGraph.postProcComponent, 'yVelocity')
    yLabelString = 'y-velocity component u_y ';
elseif strcmp(propGraph.postProcComponent, 'pressure')
    yLabelString = 'pressure p (Pa)';
elseif strcmp(propGraph.postProcComponent,'2normVelocity')
    yLabelString = ' (m/s)';
end
ylabel(yLabelString);
title(sprintf('Evaluation point X = (%d, %d)', x, y));
    
%% 5. Update the figure handle index
index = propGraph.index + 1;

end