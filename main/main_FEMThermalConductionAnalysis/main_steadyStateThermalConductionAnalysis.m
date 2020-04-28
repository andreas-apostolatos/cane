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
% Task : Steady-state thermal conduction analysis
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

% Add all functions related to body forces
addpath('../../FEMPlateInMembraneActionAnalysis/loads/');

% Add all functions related to heat transfer analysis
addpath('../../FEMThermalConductionAnalysis/solvers/',...
        '../../FEMThermalConductionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMThermalConductionAnalysis/loads/',...
        '../../FEMThermalConductionAnalysis/graphics/',...
        '../../FEMThermalConductionAnalysis/output/',...
        '../../FEMThermalConductionAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMThermalConductionAnalysis/';
% caseName = 'steadyStateSquareCavity';
% caseName = 'steadyStateWallConduction';
caseName = 'rectangularPlateWithTwoHoles';
% caseName = 'rectangularPlateWithCenterHole';

% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, ~, propGaussInt] = ...
    parse_ThermalModelFromGid(pathToCase, caseName, 'outputEnabled');

%% GUI

% On the body force computation
computeBodyForces = 'undefined';

% Choose solver for the linear equation system
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% On the writing the output function
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMThermalConductionAnalysisToVTK;
propVTK.VTKResultFile = 'undefined';

% Choose computation of the stiffness matrix
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMThermalConductionAnalysisCST;

% Linear analysis
propHeatDynamics = 'undefined';

% Initialize graphics index
graph.index = 1;

% Assign load (computeConstantFlux)
propNBC.flux = 5e4;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMHeatTransferAnalysis/';

%% Initialize solution
numNodes = length(strMsh.nodes(:,1));
numDOFs = numNodes;
dHat = zeros(numDOFs,1);

%% Solve the plate in membrane action problem
[dHat, FComplete, minElSize] = solve_FEMThermalConductionSteadyState ...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    propNBC, computeBodyForces, parameters, computeStiffMtxLoadVct, ...
    solve_LinearSystem, propNLinearAnalysis, propGaussInt, propVTK, ...
    caseName, pathToOutput, 'outputEnabled');

%% Define a function to compute analytical results 
if strcmp(caseName, 'steadyStateSquareCavity')
    propPostproc.computeAnalytical = @(x,y,t,propPostproc) ...
        propPostproc.T1+(propPostproc.T2-propPostproc.T1) * (2/pi) * ...
        sum( ((( (-1).^( (1:propPostproc.k) +1) )+1)./ (1:propPostproc.k) ) .* ...
        sin(( (1:propPostproc.k) *pi*x)/propPostproc.width) .* ...
        ( sinh(( (1:propPostproc.k) *pi*y)/propPostproc.width) ./ ...
        sinh(( (1:propPostproc.k) *pi*propPostproc.height)/propPostproc.width) ) );
    
    % Assign temperatures
    propPostproc.T1 = min(valuesInhomDOFs);
    propPostproc.T2 = max(valuesInhomDOFs);
end

%% Visualize analytical solution
if strcmp(caseName, 'steadyStateSquareCavity')
    % At steady state time is infinite
    t = inf;
    graph.index = plot_analyticalSolutionAtTimeInstance...
        (strMsh,t,propPostproc,graph,'outputEnabled');
end

%% END OF THE SCRIPT