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
% Task : Transient heat transfer analysis for different GiD input files
%
% Date : 17.04.2020
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

% Add all functions related to heat transfer analysis
addpath('../../FEMHeatTransferAnalysis/solvers/',...
        '../../FEMHeatTransferAnalysis/solutionMatricesAndVectors/',...
        '../../FEMHeatTransferAnalysis/loads/',...
        '../../FEMHeatTransferAnalysis/graphics/',...
        '../../FEMHeatTransferAnalysis/output/',...
        '../../FEMHeatTransferAnalysis/initialConditions/',...
        '../../FEMHeatTransferAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMHeatTransferAnalysis/';
% caseName = 'transientBenchmark_2w1h';
caseName = 'transientBenchmark_test';


% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    propParameters, propNLinearAnalysis, propStrDynamics, propGaussInt] = ...
    parse_HeatModelFromGid(pathToCase, caseName, 'outputEnabled');

%% UI

% On the computation of the body forces
computeBodyForces = @computeConstantVerticalStructureBodyForceVct;

% Equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Function handle to the computation of the initial conditions
computeInitCnds = @computeInitCndsFEMHeatTransferAnalysis;
propStrDynamics.initialTemperature = 0;

% Assign load
propNBC.tractionLoadVct = [1e5; 0; 0];
%computeConstantFlux

% Assign the function handles for the computation of the stiffness matrix 
% and choose solver for the finite element system based on whether the 
% analysis is linear or nonlinear
isLinear = true;
if ~isempty(propNLinearAnalysis)
    if ~ischar(propNLinearAnalysis)
        if isfield(propNLinearAnalysis, 'method')
            if strcmp(propNLinearAnalysis.method, 'NEWTON_RAPHSON')
                isLinear = false;
            end
        else
            error('Structure propNLinearAnalysis should define member variable method')
        end
    end
end
if isLinear
    computeProblemMatricesSteadyState = @computeStiffMtxAndLoadVctFEMHeatTransferAnalysisCST;
    solve_FEMSystem = @solve_FEMLinearSystem;
else
    computeProblemMatricesSteadyState = @computeStiffMtxAndLoadVctFEMHeatTransferAnalysisCST;
    solve_FEMSystem = @solve_FEMNLinearSystem;
end

% On the writing the output function
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMHeatTransferAnalysisToVTK;
propVTK.VTKResultFile = 'undefined';

% On transient inhomogeneous Dirichlet boundary conditions
updateInhomDOFs = 'undefined';
propIDBC = [];

% Not a unit test case
isUnitTest = false;

% Choose the matric computation corresponding to the chosen time
% integration scheme
if strcmp(propStrDynamics.method,'EXPLICIT_EULER')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsExplicitEuler;
    propStrDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propStrDynamics.method,'BOSSAK')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossak;
    propStrDynamics.computeUpdatedVct = ...
        @computeBossakTransientUpdatedVctAccelerationField;
else
    error('Invalid time integration method selected in propStrDynamics.method as %s',propStrDynamics.method);
end

% Choose time integration scheme parameters
if strcmp(propStrDynamics.method,'BOSSAK')
    propStrDynamics.alphaB = -.1; % -.1
    propStrDynamics.betaB = .5; % .5
    propStrDynamics.gammaB = .6; % .6
end

% Initialize graphics index
graph.index = 1;

%% Solve the transient heat transfer problem
[dHistory, minElSize] = solve_FEMHeatTransferTransient ...
    (strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateInhomDOFs, propNBC, @computeLoadVctFEMHeatTransferAnalysis, ...
    propParameters, computeBodyForces, propAnalysis, computeInitCnds, ...
    computeProblemMatricesSteadyState, propNLinearAnalysis, propIDBC, ...
    propStrDynamics, solve_LinearSystem, solve_FEMSystem, propGaussInt, ...
    propVTK, caseName, isUnitTest, 'outputEnabled');

%% Postprocessing

%% END OF THE SCRIPT
