function [timeSpaceDiscrete, resultantNumerical, resultantAnalytical] = ...
    computeResultantAtPointOverTime ...
    (x, y, fldMsh, parameters, upHistory, propFldDynamics, ...
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
% transient analysis at a specific point for 2D the incompressible flow 
% problem using FEM at the chosen Cartesian location.
%
%                     Input :
%                       x,y : Cartesian coordinates of the evaluation point
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
%              propPostproc : Structure containing information on the 
%                             postprocessing,
%                               .postProcComponent : Postprocessing 
%                                                    component to
%                                                    visualize :
%                                                    'xVelocity'
%                                                    'yVelocity'
%                                                    'pressure'
%                                                    '2normVelocity'
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
if ~isempty(propPostproc)
    if ~ischar(propPostproc)
        if isfield(propPostproc, 'postProcComponent')
            if ischar(propPostproc.postProcComponent)
                if ~strcmp(propPostproc.postProcComponent, 'xVelocity') && ...
                        ~strcmp(propPostproc.postProcComponent, 'yVelocity') && ...
                        ~strcmp(propPostproc.postProcComponent, 'pressure') && ...
                        ~strcmp(propPostproc.postProcComponent, '2normVelocity')
                    error('propGraph.postProcComponent must be one of the following character arrays, xVelocity, yVelocity, pressure, 2normVelocity');
                end
            else
                error('propGraph.postProcComponent must be a character array');
            end
        else
            error('propGraph must define field postProcComponent');
        end
        if ~isfield(propPostproc, 'computeAnalytical')
            error('propGraph must define field computeAnalytical');
        end
    else
        error('propGraph must not be a string');
    end
else
    error('propGraph must not be an empty array');
end
isAnalytical = false;
if isa(propPostproc.computeAnalytical, 'function_handle')
    isAnalytical = true;
end
if strcmp(outMsg,'outputEnabled')
    fprintf('________________________________________________\n');
    fprintf('################################################\n');
    fprintf('\n');
    fprintf('Computation of he evolution of %s\n', propPostproc.postProcComponent);
    fprintf('over time at point (%d, %d)\n', x, y);
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
resultantNumerical = zeros(propFldDynamics.noTimeSteps + 1, 1);
if isAnalytical
    resultantAnalytical = zeros(propFldDynamics.noTimeSteps + 1, 1);
else
    resultantAnalytical = 'undefined';
end
timeSpaceDiscrete = zeros(propFldDynamics.noTimeSteps + 1, 1);

% Initialize time
t = propFldDynamics.T0;

%% 1. Find the mesh element from cartesian coordinates
for iElement = 1:size(fldMsh.elements,1)

    % Find the node IDs of an element
    node_IDs = fldMsh.elements(iElement,2:end);
    
    % Find the node coordinates of an element
    vertexI = fldMsh.nodes(node_IDs(1),2:3);
    vertexJ = fldMsh.nodes(node_IDs(2),2:3);
    vertexK = fldMsh.nodes(node_IDs(3),2:3);
    
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
DOFsI = (node_IDs(1, 1)*3 - 2):(node_IDs(1, 1)*3);
DOFsJ = (node_IDs(1, 2)*3 - 2):(node_IDs(1, 2)*3);
DOFsK = (node_IDs(1, 3)*3 - 2):(node_IDs(1, 3)*3);

%% 3. Loop over all the time steps
for iTime = 1:propFldDynamics.noTimeSteps + 1
    %% 3i. Fill up the time step array accordingly
    timeSpaceDiscrete(iTime, 1) = t;
    
    %% 3ii. Get the current discrete solution vector
    upCurrent = upHistory(:, iTime);
    
    %% 3iii. Get the DOFs of the resultant at the element where (x, y) belongs to at the given time step
    upVector = upCurrent(DOFsI, 1)*N(1) + ...
               upCurrent(DOFsJ, 1)*N(2) + ...
               upCurrent(DOFsK, 1)*N(3);
           
    %% 3iv. Compute the numerical resultant at (x, y) at the given time step
    if strcmp(propPostproc.postProcComponent, 'xVelocity')
        resultantNumerical(iTime, 1) = upVector(1, 1);
    elseif strcmp(propPostproc.postProcComponent, 'yVelocity')
        resultantNumerical(iTime, 1) = upVector(2, 1);
    elseif strcmp(propPostproc.postProcComponent, 'pressure')
        resultantNumerical(iTime, 1) = upVector(3, 1);
    elseif strcmp(propPostproc.postProcComponent, '2normVelocity')
        resultantNumerical(iTime, 1) = norm(upVector(1:2, 1));
    end
    
    %% 3v. Compute the analytical resultant at (x, y) at the given time step
    if isAnalytical
        resultantAnalytical(iTime, 1) = ...
            propPostproc.computeAnalytical(x, y, t, parameters);
    end
   
    %% 3vi. Update the simulation time
    t = t + propFldDynamics.dt;
end

%% 4. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Computation of %s over time took %.2d seconds \n\n', propPostproc.postProcComponent, computationalTime);
    fprintf('_______________Postprocessing Computation Ended________________\n');
    fprintf('###############################################################\n\n\n');
end

end