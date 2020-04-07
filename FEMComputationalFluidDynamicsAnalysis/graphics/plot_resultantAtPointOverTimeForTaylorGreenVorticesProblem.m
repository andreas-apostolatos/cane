function index = plot_resultantAtPointOverTimeForTaylorGreenVorticesProblem ...
    (xi, p, Xi, eta, q, Eta, CP, isNURBS, parameters, upHistory1, upHistory2, ...
    propFldDynamics, propGraph)
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

% Number of knots in u,v-direction
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);

% Number of Control Points in u,v-direction
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Check input
checkInputForBSplineSurface ...
    (p, numKnots_xi, numCPs_xi, q, numKnots_eta, numCPs_eta);

% Initialize figure handle
figure(propGraph.index)

% Initialize the output arrays
resultantAtPointNumerical1 = zeros(propFldDynamics.noTimeSteps, 1);
resultantAtPointNumerical2 = zeros(propFldDynamics.noTimeSteps, 1);
resultantAtPointAnalytical = zeros(propFldDynamics.noTimeSteps, 1);
timeSteps = zeros(propFldDynamics.noTimeSteps, 1);

% Initialize time
t = propFldDynamics.TStart;

%% 1. Find the knot spans where the point lives in
xiSpan = findKnotSpan(xi, Xi, numCPs_xi);
etaSpan = findKnotSpan(eta, Eta, numCPs_eta);

%% 2. Compute the Cartesian coordinates of the Point
RMtx = computeIGABasisFunctionsAndDerivativesForSurface ...
    (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, 0);
xCartesian = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
    (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, RMtx);
x = xCartesian(1);
y = xCartesian(2);

%% 3. Loop over all the time steps
for iTime = 1:propFldDynamics.noTimeSteps
    %% 3i. Fill up the time step array accordingly
    timeSteps(iTime, 1) = t;
    
    %% 3ii. Get the current discrete solution vector
    
    % Distribute the solution vector into the elememts for the current time
    % step :
    
    % element discrete solution vectors
    upEl1 = zeros(numKnots_xi - p - 1, numKnots_eta - q - 1, 3*(p + 1)*(q + 1));
    upEl2 = zeros(numKnots_xi - p - 1, numKnots_eta - q - 1, 3*(p + 1)*(q + 1));
    for j = (q + 1):(numKnots_eta - q - 1)
        for i = (p + 1):(numKnots_xi - p - 1)
            % Initialize the counter
            k = 1;
            for c = j - q - 1:j - 1 
                for b = i - p:i
                    % For solution 1
                    upEl1(i, j, k) = upHistory1(3*(c*numCPs_xi + b) - 2, iTime);
                    upEl1(i, j, k + 1) = upHistory1(3*(c*numCPs_xi + b) - 1, iTime);
                    upEl1(i, j, k + 2) = upHistory1(3*(c*numCPs_xi + b), iTime);
                    
                    % For solution 2
                    upEl2(i, j, k) = upHistory2(3*(c*numCPs_xi + b) - 2, iTime);
                    upEl2(i, j, k + 1) = upHistory2(3*(c*numCPs_xi + b) - 1, iTime);
                    upEl2(i, j, k + 2) = upHistory2(3*(c*numCPs_xi + b), iTime);

                    % Update the counter
                    k = k + 3;
                end
            end
        end
    end
    
    %% 3iii. Get the actual discrete solution vector
    upActual1 = upEl1(xiSpan, etaSpan, :);
    upActual2 = upEl2(xiSpan, etaSpan, :);
    upActualVector1 = zeros(3*(p + 1)*(q + 1), 1);
    upActualVector2 = zeros(3*(p + 1)*(q + 1), 1);
    for i = 1:3*(p + 1)*(q + 1)
        % For solution 1
        upActualVector1(i) = upActual1(1, 1, i);
        
        % For solution 2
        upActualVector2(i) = upActual2(1, 1, i);
    end
    
    %% 3iv. Get the actual resultants at the point and at the current time step
    
    % For method 1
    upVector1 = computeNodalVectorIncompressibleFlow2D ...
        (RMtx, p, q, upActualVector1);
    
    % For method 2
    upVector2 = computeNodalVectorIncompressibleFlow2D ...
        (RMtx, p, q, upActualVector2);
    
    %% 3v. Get the desirable resultant at the point and at the current time step
    if strcmp(propGraph.postProcComponent, 'xVelocity')
        resultantAtPointNumerical1(iTime, 1) = upVector1(1);
        resultantAtPointNumerical2(iTime, 1) = upVector2(1);
        resultantAtPointAnalytical(iTime, 1) = -cos(x)*sin(y)*exp(-2*t*parameters.nue);
    elseif strcmp(propGraph.postProcComponent, 'yVelocity')
        resultantAtPointNumerical1(iTime, 1) = upVector1(2);
        resultantAtPointNumerical2(iTime, 1) = upVector2(2);
        resultantAtPointAnalytical(iTime, 1) = sin(x)*cos(y)*exp(-2*t*parameters.nue);
    elseif strcmp(propGraph.postProcComponent, 'pressure')
        resultantAtPointNumerical1(iTime, 1) = upVector1(3);
        resultantAtPointNumerical2(iTime, 1) = upVector2(3);
        resultantAtPointAnalytical(iTime, 1) = -.25*(cos(2*x) + cos(2*y))*exp(-4*t*parameters.nue);
    elseif strcmp(propGraph.postProcComponent, '2normVelocity')
        uXNumerical1 = upVector1(1);
        uYNumerical1 = upVector1(2);
        resultantAtPointNumerical1(iTime, 1) = norm([uXNumerical1 uYNumerical1]);
        uXNumerical2 = upVector2(1);
        uYNumerical2 = upVector2(2);
        resultantAtPointNumerical2(iTime, 1) = norm([uXNumerical2 uYNumerical2]);
        uXAnalytical = -cos(x)*sin(y)*exp(-2*t*parameters.nue);
        uYAnalytical = sin(x)*cos(y)*exp(-2*t*parameters.nue);
        resultantAtPointAnalytical(iTime, 1) = norm([uXAnalytical uYAnalytical]);
    end
   
    %% 3vi. Update the simulation time
    t = t + propFldDynamics.dt;
end

%% 4. Plot the resultant at the given point and over time
plot(timeSteps, resultantAtPointAnalytical, 'b', timeSteps, ...
    resultantAtPointNumerical1, 'g', timeSteps, resultantAtPointNumerical2, 'r');
legend('Analytical', 'Solution 1', 'Solution  2');
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
title(sprintf('Evaluation point X = (%d, %d, %d)', xCartesian(1), ...
    xCartesian(2), xCartesian(3)));
    
%% 5. Update the figure handle index
index = propGraph.index + 1;

end