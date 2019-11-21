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

%% GUI

% On the body forces
computeBodyForces = @computeConstantVerticalBodyForceVct;

% On the initial conditions
% computeInitialConditions = @computeInitialConditionsFromVTKFileFEM4NSE2D;
computeInitialConditions = @computeNullInitialConditionsFEM4NSE2D;

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
% caseName = 'flowAroundCylinderAdaptiveSteadyStateALE';
caseName = 'NACA2412_AoA5_CFD';
% caseName = 'unitTest_flowAroundCylinderAdaptiveSteadyState';

%% Parse the data from the GiD input file
[fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,NBC,analysis,parameters,...
    propNLinearAnalysis,propFldDynamics,gaussInt] = ...
    parse_FluidModelFromGid...
    (pathToCase,caseName,'outputEnabled');

%% GUI

% On the transient analysis properties
if strcmp(propFldDynamics.method,'bossak')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE2D;
end 

%% Choose the equation system solver
if strcmp(analysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
end

%% Define the name of the vtk file from where to resume the simulation
% VTKResultFile = '_contourPlots_75';
VTKResultFile = 'undefined';

%% Solve the CFD problem
[up,FComplete,hasConverged,minElSize] = solve_FEMVMSStabSteadyStateNSE2D...
    (fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,parameters,...
    computeBodyForces,analysis,computeInitialConditions,...
    VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
    gaussInt,caseName,'outputEnabled');

%% END OF THE SCRIPT
