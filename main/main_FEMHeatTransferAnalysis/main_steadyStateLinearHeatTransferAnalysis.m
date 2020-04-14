%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Heat transfer analysis for different GiD input files
%
% Date : 14.04.2020
%
%% Preamble

% Clear memory
clear;

% Clear the command window
clc;

% Close all generated windows
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

% Add all functions related to plate in membrane action analysis
addpath('../../FEMHeatTransferAnalysis/solvers/',...
        '../../FEMHeatTransferAnalysis/solutionMatricesAndVectors/',...
        '../../FEMHeatTransferAnalysis/loads/',...
        '../../FEMHeatTransferAnalysis/graphics/',...
        '../../FEMHeatTransferAnalysis/output/',...
        '../../FEMHeatTransferAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMHeatTransferAnalysis/';
caseName = 'test';

% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, propStrDynamics, propGaussInt] = ...
    parse_HeatModelFromGid(pathToCase, caseName, 'outputEnabled');

%% GUI

% On the body forces
computeBodyForces = @computeConstantVerticalStructureBodyForceVct;
% ADD A NEW FUNCTION

% Choose solver for the linear equation system
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Output properties
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMPlateInMembraneActionToVTK;
propVTK.VTKResultFile = 'undefined';

% Initialize graphics index
graph.index = 1;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Visualization of the configuration
% graph.index = plot_referenceConfigurationFEMPlateInMembraneAction(strMsh,analysis,F,homDBC,graph,'outputEnabled');

%% Initialize solution
numNodes = length(strMsh.nodes(:,1));
numDOFs = 2*numNodes;
dHat = zeros(numDOFs,1);

%% Solve the plate in membrane action problem
[dHat, FComplete, minElSize] = ...
    solve_FEMPlateInMembraneActionNLinear...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, ...
    computeBodyForces, parameters, solve_LinearSystem, propNLinearAnalysis, ...
    propGaussInt, propVTK, caseName, pathToOutput, 'outputEnabled');

%% Postprocessing
% graph.visualization.geometry = 'reference_and_current';
% resultant = 'stress';
% component = 'y';
% graph.index = plot_currentConfigurationAndResultants(strMsh,homDBC,dHat,parameters,analysis,resultant,component,graph);

%% END OF THE SCRIPT
