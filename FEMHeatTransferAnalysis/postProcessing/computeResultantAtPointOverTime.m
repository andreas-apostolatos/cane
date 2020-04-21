function [timeSpaceDiscrete, resultantNumerical, resultantAnalytical] = ...
    computeResultantAtPointOverTime ...
    (x, y, mesh, upHistory, propHeatDynamics, ...
    propPostproc, outMsg)
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
% Returns the evolution of the selected postprocessing resultant using a 
% numerical simulation and possibly the analytical solution throught the 
% transient analysis at a specific point for the 2D heat transfer problem 
% using FEM at the chosen Cartesian location.
%
%                     Input :
%                       x,y : Cartesian coordinates of the evaluation point
%                      mesh : Nodes and elements for the structure mesh
%                parameters : The parameters of the flow (density, 
%                             viscosity)
%                 upHistory : The history data of the transient analysis by
%                             solving the problem numerically
%          propHeatDynamics : Transient analysis parameters : 
%                                  T0 : Start time of the simulation
%                                TEnd : End time of the simulation
%                         noTimeSteps : Number of time steps
%                                  dt : Time step (automatically computed)
%              propPostproc : Structure containing information on the 
%                             postprocessing,
%                               .computeAnalytical : Function handle to the
%                                                    computation of the
%                                                    analtical resultant
%                   outMsg : On printing information during analysis in the
%                            command window
%
%                   Output :
%        timeSpaceDiscrete : Vector containing the discrete time space
%       resultantNumerical : Vector containing the values of the resultant
%                            over time from the numerical simulation
%      resultantAnalytical : Vector containing the values of the resultant
%                            over time from the analytical expression
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
% ->
%    3i. Fill up the time step array accordingly
%
%   3ii. Get the current discrete solution vector
%
%  3iii. Get the DOFs of the resultant at the element where (x, y) belongs to at the given time step
%
%   3iv. Compute the numerical resultant at (x, y) at the given time step
%
%    3v. Compute the analytical resultant at (x, y) at the given time step
%
%   3vi. Update the simulation time
% <-
%
% 4. Appendix
%
%% Function main body

isAnalytical = false;
if isa(propPostproc.computeAnalytical, 'function_handle')
    isAnalytical = true;
end

if strcmp(outMsg,'outputEnabled')
    fprintf('________________________________________________\n');
    fprintf('################################################\n');
    fprintf('\n');
    fprintf('Computation of Temperature evolution over time\n');
    fprintf('at point (%d, %d)\n', x, y);
    fprintf('has been initiated\n');
if isAnalytical
    fprintf('Analytical solution is computed using function\n');
    fprintf('%s\n', func2str(propPostproc.computeAnalytical));
end 
    fprintf('\n');
    fprintf('________________________________________________\n');
    fprintf('\n');
    tic;
end

%% 0. Read input
    
% Initialize the output arrays
resultantNumerical = zeros(propHeatDynamics.noTimeSteps + 1, 1);
if isAnalytical
    resultantAnalytical = zeros(propHeatDynamics.noTimeSteps + 1, 1);
else
    resultantAnalytical = 'undefined';
end
timeSpaceDiscrete = zeros(propHeatDynamics.noTimeSteps + 1, 1);

% Initialize time
t = propHeatDynamics.T0;

% Define the number of infinite series sum
propPostproc.k = 200;

% Lower left corner
x0(1) = min(mesh.nodes(:,1));
x0(2) = min(mesh.nodes(:,2));

% Lower right corner
x1(1) = max(mesh.nodes(:,1));
x1(2) = min(mesh.nodes(:,2));

% Upper left corner
x2(1) = min(mesh.nodes(:,1));
x2(2) = max(mesh.nodes(:,2));

% Get height and width
propPostproc.height = abs(x2(2)-x0(2));
propPostproc.width = abs(x1(1)-x0(1));

%% 1. Find the mesh element from cartesian coordinates
for iElement = 1:size(mesh.elements,1)

    % Find the node IDs of an element
    node_IDs = mesh.elements(iElement,:);
    
    % Find the node coordinates of an element
    vertexI = mesh.nodes(node_IDs(1),1:2);
    vertexJ = mesh.nodes(node_IDs(2),1:2);
    vertexK = mesh.nodes(node_IDs(3),1:2);
    
    % Compute basis functions and check if the point is inside the element
    [N, ~, isInside] = computeCST2DBasisFunctions ...
        (vertexI, vertexJ, vertexK, x, y);
    
    % Break the loop if the point is inside the element
    if isInside
       break
    end
end
if ~isInside
    error('Cartesian location (%d, %d) is outside the computational domain', x, y);
end

%% 2. Compute the DOFs that corespond to element nodes
DOFsI = node_IDs(1, 1);
DOFsJ = node_IDs(1, 2);
DOFsK = node_IDs(1, 3);

%% 3. Loop over all the time steps
for iTime = 1:propHeatDynamics.noTimeSteps + 1
    %% 3i. Fill up the time step array accordingly
    timeSpaceDiscrete(iTime, 1) = t;
    
    %% 3ii. Get the current discrete solution vector
    upCurrent = upHistory(:, iTime);
    
    %% 3iii. Get the DOFs of the resultant at the element where (x, y) belongs to at the given time step
    upVector = upCurrent(DOFsI, 1)*N(1) + ...
               upCurrent(DOFsJ, 1)*N(2) + ...
               upCurrent(DOFsK, 1)*N(3);
           
    %% 3iv. Compute the numerical resultant at (x, y) at the given time step
    resultantNumerical(iTime, 1) = upVector(1, 1);
    
    %% 3v. Compute the analytical resultant at (x, y) at the given time step
    if isAnalytical
        resultantAnalytical(iTime, 1) = ...
            propPostproc.computeAnalytical(x, y, t, propPostproc);
    end
   
    %% 3vi. Update the simulation time
    t = t + propHeatDynamics.dt;
end

%% 4. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Computation of Temperature field over time took %.2d seconds \n\n', computationalTime);
    fprintf('_______________Postprocessing Computation Ended________________\n');
    fprintf('###############################################################\n\n\n');
end

end