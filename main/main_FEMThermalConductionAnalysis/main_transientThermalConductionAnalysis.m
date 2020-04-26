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
%caseName = 'transientWallHeating';


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

% On the writing the output function
propVTK.isOutput = true;
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
propHeatDynamics.initialTemperature = 300;

%% Define boundary flux (load) value (computeConstantFlux)
propNBC.tractionLoadVct = [0 
                           0
                           0];

%% Solve the transient heat transfer problem
[dHistory, minElSize] = solve_FEMHeatTransferTransient ...
    (strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateInhomDOFs, propNBC, @computeLoadVctFEMHeatTransferAnalysis, ...
    propParameters, computeBodyForces, propAnalysis, computeInitialConditions, ...
    @computeStiffMtxAndLoadVctFEMHeatTransferAnalysisCST,...
    propNLinearAnalysis, propIDBC, propHeatDynamics, solve_LinearSystem, ...
    @solve_FEMLinearSystem, propGaussInt, propVTK, caseName,'outputEnabled');

%% END OF THE SCRIPT