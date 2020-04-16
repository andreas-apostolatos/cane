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
% Task : Steady state heat transfer analysis for different GiD input files
%
% Date : 14.04.2020
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

% Add all functions related to heat transfer analysis
addpath('../../FEMHeatTransferAnalysis/solvers/',...
        '../../FEMHeatTransferAnalysis/solutionMatricesAndVectors/',...
        '../../FEMHeatTransferAnalysis/loads/',...
        '../../FEMHeatTransferAnalysis/graphics/',...
        '../../FEMHeatTransferAnalysis/output/',...
        '../../FEMHeatTransferAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMHeatTransferAnalysis/';
%caseName = 'steadyStateBenchmark_1w1h';
%caseName = 'steadyStateBenchmark_2w1h';
caseName = 'plateWithTwoHoles';

% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, ~, propGaussInt] = ...
    parse_HeatModelFromGid(pathToCase, caseName, 'outputEnabled');

%% GUI

% On the body forces
computeBodyForces = @computeConstantVerticalStructureBodyForceVct;

% Choose solver for the linear equation system
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Output properties
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMHeatTransferAnalysisToVTK;
propVTK.VTKResultFile = 'undefined';

% Choose computation of the stiffness matrix
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMHeatTransferAnalysisCST;

% Linear analysis
propStrDynamics = 'undefined';

% Initialize graphics index
graph.index = 1;

% Assign load
%computeConstantFlux

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMHeatTransferAnalysis/';

%% Compute the flux(load) vector
t = 0;
F = computeLoadVctFEMHeatTransferAnalysis ...
    (strMsh, propAnalysis, propNBC, t, propGaussInt, 'outputEnabled');

%% Visualization of the configuration
% graph.index = plot_referenceConfigurationFEMPlateInMembraneAction ...
%     (strMsh, propAnalysis, F, homDOFs, [], graph, 'outputEnabled');

%% Initialize solution
numNodes = length(strMsh.nodes(:,1));
numDOFs = numNodes;
dHat = zeros(numDOFs,1);

%% Solve the plate in membrane action problem
[dHat, FComplete, minElSize] = solve_FEMHeatTransferSteadyState ...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    propNBC, computeBodyForces, parameters, computeStiffMtxLoadVct, ...
    solve_LinearSystem, propNLinearAnalysis, propGaussInt, propVTK, ...
    caseName, pathToOutput, 'outputEnabled');

%% Postprocessing

% Show analytical solution if we have a benchmark case
if strcmp(caseName, 'steadyStateBenchmark_1w1h') || strcmp(caseName, 'steadyStateBenchmark_2w1h')
    graph.index = plot_steadyStateBenchmarkProblemAnalyticalSolution(strMsh,valuesInhomDOFs,graph,'outputEnabled');
end

%% END OF THE SCRIPT
