function index = plot_resultantAtPointOverTimeForTaylorGreenVorticesProblem ...
    (x, y, parameters, upHistory, propFldDynamics, propGraph)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documenation
%
% Returns a 2D plot with the evolution of a resultant throught the
% transient analysis at a specific point for two different solutions
% corresponding to a 2D incompressible flow problem using isogeometric 
% analysis for the benchmark of Taylor-Green votices together with the 
% analytical resultant at the chosen point.
%
%                     Input :
%                    xi,eta : The NURBS parametric coordinates of the
%                             points on which to plot the resultants over 
%                             the time
%                       p,q : The polynomial degrees of the NURBS patch
%                    Xi,Eta : The knot vectors of the NURBS patch
%                        CP : The set of Control Point coordinates and
%                             weights for the NURBS patch
%                   isNURBS : Function handle to whether the basis of the
%                             B-Spline patch is a NURBS or a B-Spline
%                parameters : The parameters of the flow (density, 
%                             viscosity)
%                upHistory1 : The history data of the transient analysis by
%                             solving the problem numerically with method 1
%                upHistory2 : The history data of the transient analysis by
%                             solving the problem numerically with method 2
%           propFldDynamics : Transient analysis parameters : 
%                              TStart : Start time of the simulation
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
%                   Output :
%                            graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Find the knot spans where the point lives in
%
% 2. Compute the Cartesian coordinates of the Point
%
% 3. Loop over all the time steps
%
%    3i. Fill up the time step array accordingly
%
%   3ii. Get the current discrete solution vector
%
%  3iii. Get the actual discrete solution vector
%
%   3iv. Get the actual resultants at the point and at the current time step
%
%    3v. Get the desirable resultant at the point and at the current time step
%
%   3vi. Update the simulation time
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
t = propFldDynamics.T0;

%% 3. Loop over all the time steps
for iTime = 1:propFldDynamics.noTimeSteps
    %% 3i. Fill up the time step array accordingly
    timeSteps(iTime, 1) = t;
    
    %% 3ii. Get the current discrete solution vector
    
    % Distribute the solution vector into the elememts for the current time
    % step :
    
    % element discrete solution vectors
  
    upEl = zeros(numKnots_xi - p - 1, numKnots_eta - q - 1, 3*(p + 1)*(q + 1));
    for j = (q + 1):(numKnots_eta - q - 1)
        for i = (p + 1):(numKnots_xi - p - 1)
            % Initialize the counter
            k = 1;
            for c = j - q - 1:j - 1 
                for b = i - p:i
                    
                    
                    % For solution 2
                    upEl(i, j, k) = upHistory(3*(c*numCPs_xi + b) - 2, iTime);
                    upEl(i, j, k + 1) = upHistory(3*(c*numCPs_xi + b) - 1, iTime);
                    upEl(i, j, k + 2) = upHistory(3*(c*numCPs_xi + b), iTime);

                    % Update the counter
                    k = k + 3;
                end
            end
        end
    end
    
    %% 3iii. Get the actual discrete solution vector
   
    upActual = upEl(xiSpan, etaSpan, :);
   
    upActualVector = zeros(3*(p + 1)*(q + 1), 1);
    for i = 1:3*(p + 1)*(q + 1)
        
        
        % For solution 2
        upActualVector(i) = upActual(1, 1, i);
    end
    
    %% 3iv. Get the actual resultants at the point and at the current time step
    
    
    
    % For method 2
    upVector = computeNodalVectorIncompressibleFlow2D ...
        (RMtx, p, q, upActualVector);
    
    %% 3v. Get the desirable resultant at the point and at the current time step
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
   
    %% 3vi. Update the simulation time
    t = t + propFldDynamics.dt;
end

%% 4. Plot the resultant at the given point and over time
plot(timeSteps, resultantAtPointAnalytical, 'b',...
     timeSteps, resultantAtPointNumerical, 'r');
legend('Analytical', 'Navier-Stokes');
xlabel('time (seconds)');
if strcmp(propGraph.postProcComponent, 'xVelocity')
    yLabelString = 'x-velocity component u_x (m/s)';
elseif strcmp(propGraph.postProcComponent, 'yVelocity')
    yLabelString = 'y-velocity component u_y (m/s)';
elseif strcmp(propGraph.postProcComponent, 'pressure')
    yLabelString = 'pressure p (Pa)';
elseif strcmp(propGraph.postProcComponent,'2normVelocity')
    yLabelString = 'velocity magnitude ||u||_2 (m/s)';
end
ylabel(yLabelString);
title(sprintf('Evaluation point X = (%d, %d)', xCartesian(1), xCartesian(2)));
    
%% 5. Update the figure handle index
index = propGraph.index + 1;

end