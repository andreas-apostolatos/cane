%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Performs Multilevel Monte-Carlo for the 2D incompressible steady-
%        state Navier-Stokes equations
%
% Date : 16.12.2019
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
        '../../FEMComputationalFluidDynamicsAnalysis/ALEMotion/',...
        '../../FEMComputationalFluidDynamicsAnalysis/postProcessing/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');

%% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_testFEM4NavierStokesSteadyStateFlowAroundCylinder2D';

%% Parse the data from the GiD input file
[fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,~,analysis,parameters,...
    propNLinearAnalysis,propFldDynamics,gaussInt,postProc] = ...
    parse_FluidModelFromGid...
    (pathToCase,caseName,'');

%% GUI
% On the body forces
computeBodyForces = @computeConstantVerticalBodyForceVct;

% On the initial conditions
% computeInitialConditions = @computeInitialConditionsFromVTKFileFEM4NSE2D;
computeInitialConditions = @computeNullInitialConditionsFEM4NSE2D;

% On the transient analysis properties
if strcmp(propFldDynamics.method,'bossak')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE2D;
end

% Sampling function selection
% 'generateRandomUniformDistribution', 'generateRandomNormalDistribution', 
% 'generateRandomLatinHypercubeDistribution',
% 'generateRandomQuasiMonteCarloDistribution'
computeRandomDistribution = @generateRandomNormalDistribution;

% Initialize graph index
graph.index = 1;

% define parameters used in reference paper and simualiton
Ubar = 0.2; % mid velocity 
D = 0.1;    % diameter of the body
rho = parameters.rho; % density

%% Choose the equation system solver
if strcmp(analysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
end

%% Define the name of the vtk file from where to resume the simulation
VTKResultFile = 'undefined';

%% Generate input samples

% Number of samples
no_samples = 100;

% Mean value of the samples
mean_value = 0.3;

% Standard deviation of the samples
standard_deviation = 0.015;

% Uncertainty in the inlet amplitude
Umax_random = computeRandomDistribution(mean_value, standard_deviation, no_samples);

% histogram(input,100)
figure(graph.index)
histfit(Umax_random)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(Umax_random))
title('Input histogram')
graph.index = graph.index + 1;

%% Initialize output statistics
drag_coefficient_random = zeros(no_samples,1);
lift_coefficient_random = zeros(no_samples,1);

%% Set up progress bar
fprintf(['\n' repmat('.',1,no_samples) '\n\n']);
tic

%% Perform a Monte-Carlo Simulation for all samples
parfor i = 1:no_samples
    %% Change the input velocity to match the reference paper - parabolic input
    Umax = Umax_random(i,1);
    valuesInhomDBCModified = computeInletVelocityParabolic_unitTest(fldMsh, inhomDBC, valuesInhomDBC, Umax);
    
    %% Solve the CFD problem
    [~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh,homDBC,inhomDBC,valuesInhomDBCModified,nodesALE,parameters,...
        computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
        gaussInt,caseName,'');
    
    %% Compute postprocessing
    postProc_random = computePostProc(FComplete, analysis, parameters, postProc);
    forcesOnDomain = postProc_random.valuePostProc{1};
    Fx = forcesOnDomain(1,1);
    Fy = forcesOnDomain(2,1);
    drag_coefficient_random(i,1) = (2 * Fx)/(rho * Ubar * Ubar * D);
    liftCoefficient_random(i,1) = abs((2 * Fy)/(rho * Ubar * Ubar * D));
    
    %% Update progress bar
    fprintf('\b|\n');
end
disp(['Elapsed time: ', num2str(toc)])

%% Perform statistics to the output variables

% Drag coefficient
figure(graph.index)
histfit(drag_coefficient_random)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(drag_coefficient_random))
title('Output histogram for the drag coefficient')
graph.index = graph.index + 1;

% Lift coefficient
figure(graph.index)
histfit(lift_coefficient_random)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(lift_coefficient_random))
title('Output histogram for the lift coefficient')
graph.index = graph.index + 1;

%% END OF THE SCRIPT