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
% Date : 28.04.2020
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
addpath('../../FEMThermalConductionAnalysis/solvers/',...
        '../../FEMThermalConductionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMThermalConductionAnalysis/loads/',...
        '../../FEMThermalConductionAnalysis/graphics/',...
        '../../FEMThermalConductionAnalysis/output/',...
        '../../FEMThermalConductionAnalysis/initialConditions/',...
        '../../FEMThermalConductionAnalysis/postprocessing/', ...
        '../../FEMThermalConductionAnalysis/transientAnalysis/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMThermalConductionAnalysis/';
caseName = 'trapezoidalPlateHeatFlux';
% caseName = 'rectangularPlateHeatFlux';

% Parse the data from the GiD input file
[thermalMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    propParameters, propNLinearAnalysis, propThermalDynamics, propGaussInt] = ...
    parse_ThermalModelFromGid(pathToCase, caseName, 'outputEnabled');

%% UI

% On the computation of the body forces
computeBodyForces = 'undefined';

% Equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% On the writing the output function
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMThermalConductionAnalysisToVTK;
propVTK.VTKResultFile = 'undefined';

% Initialize graphics index
propGraph.index = 1;

% On transient inhomogeneous Dirichlet boundary conditions
updateInhomDOFs = 'undefined';
propIDBC = [];

% On post processing
propPostproc.computeAnalytical = 'undefined';

% Define the initial condition function 
computeInitialConditions = @computeInitCndsFEMThermalConductionAnalysis;

%% Define intial temperature and/or applied heat flux
if strcmp(caseName,'trapezoidalPlateHeatFlux') || strcmp(caseName,'rectangularPlateHeatFlux')
    propThermalDynamics.temperatureInit = 300;
    propNBC.flux = 100; 
end

%% Define the integration schemes function handles
computeTransientProblemMatrices = ...
    {@computeProblemMtrcsImplicitEulerThermalConduction ...
    @computeProblemMtrcsGalerkinThermalConduction ....
    @computeProblemMtrcsCrankNicolsonThermalConduction};
name = {'Implicit' 'Galerkin' 'Crank-Nicolson'};
color = {'blue', 'red', 'black'};
nameMethods = length(computeTransientProblemMatrices);

%% Initialize solution arrays
if ~propVTK.isOutput
    numNodes = length(thermalMsh.nodes);
    THistory = zeros(numNodes, propThermalDynamics.noTimeSteps + 1, nameMethods);
end

%% Loop over all integration schemes to compare the results
for iMet = 1:length(computeTransientProblemMatrices)
    %% Define the time integration rules
    propThermalDynamics.computeProblemMtrcsTransient = ...
        computeTransientProblemMatrices{iMet};
    propThermalDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
    
    %% Solve the transient heat transfer problem
    [THistory(:, :, iMet), WComplete, minElSize] = solve_FEMThermalConductionTransient ...
    (thermalMsh, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateInhomDOFs, propNBC, @computeLoadVctFEMThermalConductionAnalysis, ...
    propParameters, computeBodyForces, propAnalysis, computeInitialConditions, ...
    @computeStiffMtxAndLoadVctFEMThermalConductionAnalysisCST, ...
    propNLinearAnalysis, propIDBC, propThermalDynamics, solve_LinearSystem, ...
    @solve_FEMLinearSystem, propGaussInt, propVTK, caseName, 'outputEnabled');
end

%% Plot the evolution of the temperature over time at the selected Cartesian location
if ~propVTK.isOutput
    figure(propGraph.index)
    hold on;
    for iMet = 1:length(computeTransientProblemMatrices)
        x = -10;
        y = 0;
        [timeSpaceDiscrete, temperature_h, ~] = ...
            computeTemperatureAtPointOverTime ...
            (x, y, thermalMsh, THistory(:, :, iMet), propThermalDynamics, ...
            propPostproc, 'outputEnabled');
        plot(timeSpaceDiscrete, temperature_h, color{iMet},'DisplayName',name{iMet});
    end
    hold off;
    xlabel('Time [seconds]');
    ylabel('Temperature');
    title(sprintf('Evolution of Temperature at point X = (%d, %d)', x, y));
    legend('Orientation', 'horizontal', 'Location', 'southoutside');  
    hold off;
    propGraph.index = propGraph.index + 1;
end

%% END OF THE SCRIPT