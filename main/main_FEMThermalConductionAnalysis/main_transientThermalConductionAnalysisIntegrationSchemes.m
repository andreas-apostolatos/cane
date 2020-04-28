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
addpath('../../FEMThermalConductionAnalysis/solvers/',...
        '../../FEMThermalConductionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMThermalConductionAnalysis/loads/',...
        '../../FEMThermalConductionAnalysis/graphics/',...
        '../../FEMThermalConductionAnalysis/output/',...
        '../../FEMThermalConductionAnalysis/initialConditions/',...
        '../../FEMThermalConductionAnalysis/postProcessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMThermalConductionAnalysis/';
caseName = 'trapezoidalPlateHeatFlux';
% caseName = 'rectangularPlateHeatFlux';

% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    propParameters, propNLinearAnalysis, propThermalDynamics, propGaussInt] = ...
    parse_ThermalModelFromGid(pathToCase, caseName, 'outputEnabled');

%% UI

% On the computation of the body forces
computeBodyForces = 'undefined';

% Equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% On the writing the output function
propVTK.isOutput = false;
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
integrationSchemes = {@computeProblemMtrcsImplicitEulerThermalConduction
                      @computeProblemMtrcsGalerkinThermalConduction
                      @computeProblemMtrcsCrankNicolsonThermalConduction};
                  
% Define name and color for each scheme
name = {'Implicit';'Galerkin';'Crank-Nicolson'};
color = {'blue';'red';'green'};
       
%% Initialize figure properties
figure(propGraph.index)
hold on

%% Loop over all integration schemes to compare the results
for n = 1:length(integrationSchemes)
    
    %% Define the time integration rules
    propThermalDynamics.computeProblemMtrcsTransient = ...
        integrationSchemes{n};
    propThermalDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
    
    %% Solve the transient heat transfer problem
    [dHistory, minElSize] = solve_FEMThermalConductionTransient ...
    (strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateInhomDOFs, propNBC, @computeLoadVctFEMThermalConductionAnalysis, ...
    propParameters, computeBodyForces, propAnalysis, computeInitialConditions, ...
    @computeStiffMtxAndLoadVctFEMThermalConductionAnalysisCST, ...
    propNLinearAnalysis, propIDBC, propThermalDynamics, solve_LinearSystem, ...
    @solve_FEMLinearSystem, propGaussInt, propVTK, caseName,'outputEnabled');
    
    %% Compute the selected resultant at the chosen Cartesian location over time 
    x = -10;
    y = 0;
    [timeSpaceDiscrete, resultantNumerical, ~] = ...
        computeResultantAtPointOverTime...
        (x, y, strMsh, dHistory, ...
        propThermalDynamics, propPostproc, 'outputEnabled');

    %% Plot the selected resultant at the chosen Cartesian location over time
    plot(timeSpaceDiscrete, resultantNumerical, color{n},'DisplayName',name{n});
    
    temperature50s(n,1) = resultantNumerical(timeSpaceDiscrete == 50);
    temperature75s(n,1) = resultantNumerical(timeSpaceDiscrete == 75);
end

%% Update graph properties
xlabel('Time [seconds]');
ylabel('Temperature');
title(sprintf('Evolution of Temperature at point X = (%d, %d)', x, y));
legend('Orientation', 'horizontal', 'Location', 'southoutside');  
hold off
propGraph.index = propGraph.index + 1;

%% Print out the temperature at selected time instances
table(name,temperature50s,temperature75s)

%% END OF THE SCRIPT