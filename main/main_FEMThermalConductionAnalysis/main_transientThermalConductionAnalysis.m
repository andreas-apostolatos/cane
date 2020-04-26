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
% Task : Transient thermal conduction analysis
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
% caseName = 'transientSquareCavity';
caseName = 'transientWallHeating';
% caseName = 'thermalConduction_trapezoid';

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
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMThermalConductionAnalysisToVTK;
propVTK.VTKResultFile = 'undefined';

% On transient inhomogeneous Dirichlet boundary conditions
updateInhomDOFs = 'undefined';
propIDBC = [];

% Choose the appropriate matrix update computation corresponding to the
% chosen time integration scheme
if strcmp(propThermalDynamics.method,'IMPLICIT_EULER')
    propThermalDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsImplicitEuler;
    propThermalDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propThermalDynamics.method,'GALERKIN')
    propThermalDynamics.computeProblemMtrcsTransient =  ...
        @computeProblemMtrcsGalerkin;
    propThermalDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propThermalDynamics.method,'CRANK_NICOLSON')
    propThermalDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsCrankNicolson;
    propThermalDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
else
    error('Invalid time integration method selected in propStrDynamics.method as %s', ...
        propThermalDynamics.method);
end

% Initialize graphics index
propGraph.index = 1;

%% Define the initial condition function
computeInitialConditions = @computeInitCndsFEMThermalConductionAnalysis;
propThermalDynamics.temperatureInit = 300; % 300

%% Define boundary flux (load) value (computeConstantFlux)
propNBC.flux = 1e2; % 1e2

%% Solve the transient heat transfer problem
[dHistory, minElSize] = solve_FEMThermalConductionTransient ...
    (strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateInhomDOFs, propNBC, @computeLoadVctFEMThermalConductionAnalysis, ...
    propParameters, computeBodyForces, propAnalysis, computeInitialConditions, ...
    @computeStiffMtxAndLoadVctFEMThermalConductionAnalysisCST,...
    propNLinearAnalysis, propIDBC, propThermalDynamics, solve_LinearSystem, ...
    @solve_FEMLinearSystem, propGaussInt, propVTK, caseName,'outputEnabled');

%% END OF THE SCRIPT