%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Geometrically nonlinear transient plane stress analysis
%
% Date : 19.02.2014
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

% Equation system solvers
addpath('../../equationSystemSolvers/');

% Equation system solvers
addpath('../../efficientComputation/');

% Include functions related to transient analysis
addpath('../../transientAnalysis/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/initialConditions/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMPlateInMembraneActionAnalysis/';

caseName = 'curvedPlateTipShearPlaneStressTransient';
% caseName = 'cantileverBeamPlaneStressTransientNLinear';

% Parse the data from the GiD input file
[strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, propStrDynamics, propGaussInt] = ...
    parse_StructuralModelFromGid(pathToCase, caseName, 'outputEnabled');

%% UI

% On the computation of the body forces
computeBodyForces = @computeConstantVerticalStructureBodyForceVct;

% Equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Function handle to the computation of the initial conditions
computeInitCnds = @computeInitCndsFEMPlateInMembraneAction;

% Define the amplitude of the externally applied load and time duration
propNBC.tractionVector = [-1e5; 0; 0];
propNBC.endTime = 1;

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
%     computeProblemMatricesSteadyState = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed;
    computeProblemMatricesSteadyState = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST;
    solve_FEMSystem = @solve_FEMLinearSystem;
else
    computeProblemMatricesSteadyState = @computeTangentStiffMtxResVctFEMPlateInMembraneAction;
    solve_FEMSystem = @solve_FEMNLinearSystem;
end

% On the writing the output function
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMPlateInMembraneActionToVTK;
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

%% Visualization of the configuration
% graph.index = plot_referenceConfigurationFEMPlateInMembraneAction(strMsh,analysis,F,homDBC,graph,'outputEnabled');

%% Solve the plate in membrane action problem
[~, minElSize] = solve_FEMPlateInMembraneActionTransient ...
    (propAnalysis, strMsh, homDBC, inhomDBC, valuesInhomDBC, ...
    updateInhomDOFs, propNBC, @computeLoadVctFEMPlateInMembraneAction, ...
    parameters, computeBodyForces, computeInitCnds, ...
    computeProblemMatricesSteadyState, propNLinearAnalysis, propIDBC, ...
    propStrDynamics, solve_LinearSystem, solve_FEMSystem, propGaussInt, ...
    propVTK, caseName, 'outputEnabled');

%% Postprocessing
% graph.visualization.geometry = 'reference_and_current';
% resultant = 'stress';
% component = 'y';
% graph.index = plot_currentConfigurationAndResultants(strMsh,homDBC,dHat,parameters,analysis,resultant,component,graph);

%% END OF THE SCRIPT
