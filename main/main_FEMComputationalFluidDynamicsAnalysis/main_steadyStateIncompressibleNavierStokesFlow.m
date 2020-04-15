%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Solves the steady-state incompressible Navier-Stokes equations in 
%        2D.
%
% Date : 03.06.2014
%
%% Preamble
clear;
clc;
close all;

%% Includes

% Add functions related to equation system solvers
addpath('../../equationSystemSolvers/');

% Add general math functions
addpath('../../generalMath/');

% Add the classical finite element basis functions
addpath('../../basisFunctions/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/');

% Add all functions related to the Finite Element Methods for Computational
% Fluid Dynamics problems
addpath('../../FEMComputationalFluidDynamicsAnalysis/solutionMatricesAndVectors/',...
        '../../FEMComputationalFluidDynamicsAnalysis/initialConditions',...
        '../../FEMComputationalFluidDynamicsAnalysis/solvers/',...
        '../../FEMComputationalFluidDynamicsAnalysis/loads/',...
        '../../FEMComputationalFluidDynamicsAnalysis/output/',...
        '../../FEMComputationalFluidDynamicsAnalysis/ALEMotion/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');

%% Parse the data from the GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
% caseName = 'flowAroundCylinderAdaptiveSteadyStateALE';
% caseName = 'NACA2412_AoA5_CFD';
caseName = 'unitTest_flowAroundCylinderAdaptiveSteadyState';

% Parse the data from the GiD input file
[fldMsh, homDBC, inhomDBC, valuesInhomDBC, nodesALE, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, propFldDynamics, propGaussInt, ~] = ...
    parse_FluidModelFromGid...
    (pathToCase, caseName, 'outputEnabled');

%% UI

% On the body forces
computeBodyForces = @computeConstantVerticalFluidBodyForceVct;

% On the initial conditions
% computeInitialConditions = @computeInitialConditionsFromVTKFileFEM4NSE2D;
computeInitialConditions = @computeNullInitialConditionsFEM4NSE;

% Output properties
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;
propVTK.VTKResultFile = 'undefined';

%% Choose the equation system solver
if strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
end

%% Initialize solution
[up, ~, ~, numIterStep] = computeInitialConditions...
    (propAnalysis, fldMsh, 'undefined', 'undefined', 'undefined', ...
    'undefined', 'undefined', 'undefined');

%% Solve the CFD problem
[up, FComplete, isConverged, minElSize] = solve_FEMVMSStabSteadyStateNSE ...
    (fldMsh, up, homDBC, inhomDBC, valuesInhomDBC, nodesALE, parameters, ...
    computeBodyForces, propAnalysis, solve_LinearSystem, propFldDynamics, ...
    propNLinearAnalysis, numIterStep, propGaussInt, propVTK, caseName, ...
    'outputEnabled');

%% END OF THE SCRIPT
