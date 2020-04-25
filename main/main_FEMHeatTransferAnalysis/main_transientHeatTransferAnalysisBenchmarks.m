%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Script documentation
%
% Task : Transient heat transfer analysis for two benchmark cases
%
% Date : 21.04.2020
%
%% Preamble
clear;
clc;
close all;
 
%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the low order basis functions
addpath('../../basisFunctions/');

% Add the equation system solvers
addpath('../../equationSystemSolvers/');

% Add all the efficient computation functions
addpath('../../efficientComputation/');

% Include functions related to transient analysis
addpath('../../transientAnalysis/');

% Add all functions related to body forces
addpath('../../FEMPlateInMembraneActionAnalysis/loads/');
 
% Add all functions related to heat transfer analysis
addpath('../../FEMHeatTransferAnalysis/solvers/',...
        '../../FEMHeatTransferAnalysis/solutionMatricesAndVectors/',...
        '../../FEMHeatTransferAnalysis/loads/',...
        '../../FEMHeatTransferAnalysis/graphics/',...
        '../../FEMHeatTransferAnalysis/output/',...
        '../../FEMHeatTransferAnalysis/initialConditions/',...
        '../../FEMHeatTransferAnalysis/postProcessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMHeatTransferAnalysis/';
caseName = 'transientSquareCavity';
% caseName = 'transientWallHeating';

% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    propParameters, propNLinearAnalysis, propHeatDynamics, propGaussInt] = ...
    parse_HeatModelFromGid(pathToCase, caseName, 'outputEnabled');

%% UI

% On the computation of the body forces
computeBodyForces = 'undefined';

% Equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Assign load
propNBC.tractionLoadVct = [0
                           0
                           0]; %computeConstantFlux

% On the writing the output function
propVTK.isOutput = false;
propVTK.writeOutputToFile = @writeOutputFEMHeatTransferAnalysisToVTK;
propVTK.VTKResultFile = 'undefined';

% On transient inhomogeneous Dirichlet boundary conditions
updateInhomDOFs = 'undefined';
propIDBC = [];

% Choose the appropriate matrix update computation corresponding to the
% chosen time integration scheme
if strcmp(propHeatDynamics.method,'IMPLICIT_EULER')
    propHeatDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsImplicitEuler;
    propHeatDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propHeatDynamics.method,'GALERKIN')
    propHeatDynamics.computeProblemMtrcsTransient =  ...
        @computeProblemMtrcsGalerkin;
    propHeatDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propHeatDynamics.method,'CRANK_NICOLSON')
    propHeatDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsCrankNicolson;
    propHeatDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
else
    error('Invalid time integration method selected in propStrDynamics.method as %s',propHeatDynamics.method);
end

% Initialize graphics index
propGraph.index = 1;

%% Define the initial condition function
computeInitialConditions = @computeInitCndsFEMHeatTransferAnalysis;
if strcmp(caseName, 'transientSquareCavity')
    propHeatDynamics.initialTemperature = 300;
elseif strcmp(caseName, 'transientWallHeating')
    propHeatDynamics.initialTemperature = 200;
end

%% Solve the transient heat transfer problem
[dHistory, minElSize] = solve_FEMHeatTransferTransient ...
    (strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateInhomDOFs, propNBC, @computeLoadVctFEMHeatTransferAnalysis, ...
    propParameters, computeBodyForces, propAnalysis, computeInitialConditions, ...
    @computeStiffMtxAndLoadVctFEMHeatTransferAnalysisCST,...
    propNLinearAnalysis, propIDBC, propHeatDynamics, solve_LinearSystem, ...
    @solve_FEMLinearSystem, propGaussInt, propVTK, caseName,'outputEnabled');

%% Define a function to compute analytical results 
if strcmp(caseName, 'transientSquareCavity')
    propPostproc.computeAnalytical = @(x,y,t,propPostproc) ...
        propPostproc.T1+(propPostproc.T2-propPostproc.T1) * (2/pi) * ...
        sum( ((( (-1).^( (1:propPostproc.k) +1) )+1)./ (1:propPostproc.k) ) .* ...
        sin(( (1:propPostproc.k) *pi*x)/propPostproc.width) .* ...
        ( sinh(( (1:propPostproc.k) *pi*y)/propPostproc.width) ./ ...
        sinh(( (1:propPostproc.k) *pi*propPostproc.height)/propPostproc.width) ) );   
    % NOTE -> this is NOT a time dependent function
    
    % Assign temperatures
    propPostproc.T1 = min(valuesInhomDOFs);
    propPostproc.T2 = max(valuesInhomDOFs);
    
elseif strcmp(caseName, 'transientWallHeating')
    propPostproc.computeAnalytical = @(x,y,t,propPostproc) ...
        propPostproc.T2 + (propPostproc.T1-propPostproc.T2)*(4/pi) * ...
        sum( (1./(2*(1:propPostproc.k)-1)) .* ...
        exp( -( ((2*(1:propPostproc.k)-1)*pi)/propPostproc.width ) *propPostproc.alpha*t ) .* ...
        sin( ((2*(1:propPostproc.k)-1)*pi*x)/propPostproc.width) );
    
    % Assign temperatures and alpha
    propPostproc.T1 = propHeatDynamics.initialTemperature;
    propPostproc.T2 = max(valuesInhomDOFs);
    propPostproc.alpha = propParameters.alpha;
end

%% Visualize analytical solution at one time instance
if strcmp(caseName, 'transientSquareCavity') || strcmp(caseName, 'transientWallHeating')
    t = propHeatDynamics.TEnd;
    propGraph.index = plot_analyticalSolutionAtTimeInstance...
        (strMsh,t,propPostproc,propGraph,'outputEnabled');
end

%% Compute the selected resultant at the chosen Cartesian location over time                                
x = 0.6;
y = 0.6;
[timeSpaceDiscrete, resultantNumerical, resultantAnalytical] = ...
    computeResultantAtPointOverTime...
    (x, y, strMsh, dHistory, ...
    propHeatDynamics, propPostproc, 'outputEnabled');

%% Plot the selected resultant at the chosen Cartesian location over time
figure(propGraph.index)
if ~ischar(resultantAnalytical)
    plot(timeSpaceDiscrete, resultantAnalytical, 'black',...
         timeSpaceDiscrete, resultantNumerical, 'blue');
    legend('Analytical', 'FEM', 'Orientation', 'horizontal', 'Location', 'southoutside');
else
    plot(timeSpaceDiscrete, resultantNumerical, 'blue');
end
xlabel('Time [seconds]');
ylabel('Temperature');
title(sprintf('Evolution of Temperature at point X = (%d, %d)', x, y));
propGraph.index = propGraph.index + 1;

%% END OF THE SCRIPT