%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Matthew Keller
%
%% Script documentation
%
% Task : Analyzes sensitivity action of a cylindrical object in 2D, steady 
%        state, incompressible Navier Stokes, channel flow. The design 
%        parameter is the cylinder radius. The sensitivities are computed 
%        using Finite Differencing whereas a gradient-based descend 
%        optimization algorithm is employed.
%
% Date : 18.04.2020
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
        '../../FEMComputationalFluidDynamicsAnalysis/boundaryConditions/',...
        '../../FEMComputationalFluidDynamicsAnalysis/solvers/',...
        '../../FEMComputationalFluidDynamicsAnalysis/loads/',...
        '../../FEMComputationalFluidDynamicsAnalysis/ALEMotion/',...
        '../../FEMComputationalFluidDynamicsAnalysis/postProcessing/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');

%% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_Cylinder2D_SensitivityCheck';

%% Parse the data from the GiD input file
[fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propALE, ~, propAnalysis, ...
    propParameters, propNLinearAnalysis, propFldDynamics, propGaussInt, ...
    propPostproc] = parse_FluidModelFromGid(pathToCase, caseName, '');

%% UI

% On the body forces
computeBodyForces = @computeConstantVerticalBodyForceVct;

% On the transient analysis properties
if strcmp(propFldDynamics.method,'bossak')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE2D;
end

% Free DOFs of the system (actual DOFs over which the solution is computed)
numNodes = length(fldMsh.nodes(:, 1));
numDOFs = propAnalysis.noFields*numNodes;
freeDOFs = 1:numDOFs;
prescribedDOFs = mergesorted(homDOFs, inhomDOFs);
prescribedDOFs = unique(prescribedDOFs);
freeDOFs(ismember(freeDOFs, prescribedDOFs)) = [];

% Define the name of the vtk file from where to resume the simulation
propVTK_true.isOutput = true;
propVTK_true.writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;
propVTK_false.isOutput = false;
propVTK_false.writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;

% Maximum optimization iterations and drag tolerance
maxIter = 30;
tolDrag = 1e-4;

% Perturbation size and design update scaling
epsilonTilde = 1e-2;
alphaTilde = 1e-2;

% Function handle for the computation of the Hessian
djdp1_ = @(djdp1_, p1) djdp1_;

% Dummy parameters
nodesSaved = 'undefined';
uMeshALE = 'undefined';

%% Choose the equation system solver
if strcmp(propAnalysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
end

%% Input parameters
% max input velocity defined in the reference paper
Umax = .1;

%% Change input velocity to have the parabolic distribution for each randomized input
% valuesInhomDBCModified = computeInletVelocityPowerLaw(fldMsh, inhomDOFs, valuesInhomDOFs, Umax);
valuesInhomDBCModified = valuesInhomDOFs;

%% Get the center of the cylinder in flow
propALE.propUser.x0 = mean(fldMsh.nodes(propPostproc.nodesDomain{1}(:, 1), 1));
propALE.propUser.y0 = mean(fldMsh.nodes(propPostproc.nodesDomain{1}(:, 1), 2));
radiusInit = max(fldMsh.nodes(propPostproc.nodesDomain{1}(:, 1), 1)) - propALE.propUser.x0;

%% Initializations

% Initialize the solution of the state equations
up = computeNullInitialConditionsFEM4NSE ...
    (propAnalysis, fldMsh, 'undefined', 'undefined', 'undefined', ... 
    'undefined', 'undefined', 'undefined');

% Initialize objective and design history
optHistory = zeros(maxIter, 2);

% Initialize the design update
dr = 0;

% Initialize the design
radius = radiusInit + dr;

% Initialize the optimization counter
counterOpt = 1;

% Print progress dots
fprintf(['\n' repmat('.',1,maxIter) '\n\n']);

% Start time count
tic;

%% Loop over all the optimization iterations
while counterOpt <= maxIter
    %% Update the mesh for the nominal state according to the updated design
    propALE.propUser.dr = dr;
    radius = radius + dr;
    [fldMsh, ~, ~, ~, ~, ~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D ...
        (fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, freeDOFs, ...
        nodesSaved, propALE, solve_LinearSystem, propFldDynamics, counterOpt);
    fldMsh.initialNodes = fldMsh.nodes;
    
    %% Solve the nominal steady-state CFD problem
    [up, FComplete, ~, ~] = solve_FEMVMSStabSteadyStateNSE ...
        (fldMsh, up, homDOFs, inhomDOFs, valuesInhomDBCModified, ...
        uMeshALE, propParameters, computeBodyForces, propAnalysis, ...
        solve_LinearSystem, propFldDynamics, propNLinearAnalysis, ...
        counterOpt, propGaussInt, propVTK_true, caseName, '');
    
    %% Compute the nominal drag force
    postProc_update = computePostProc ...
        (FComplete, propAnalysis, propParameters, propPostproc);
    forcesOnDomain = postProc_update.valuePostProc{1};
    drag_nom = forcesOnDomain(1, 1);
    
    %% Check if the objective function is below a threshold
    optHistory(counterOpt, 1) = drag_nom;
    optHistory(counterOpt, 2) = radius;
    if abs(drag_nom) < tolDrag
        break;
    end
 
    %% Move the mesh according to the prescribed perturbation
    propALE.propUser.dr = epsilonTilde;
    [fldMsh_p1, ~, ~, ~, ~, ~] = ...
        computeUpdatedMeshAndVelocitiesPseudoStrALE2D ...
        (fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, freeDOFs, ...
        nodesSaved, propALE, solve_LinearSystem, propFldDynamics, counterOpt);
    fldMsh_p1.initialNodes = fldMsh_p1.nodes;
    
    %% Solve the perturbed steady-state CFD problem
    [~, FComplete, ~, ~] = solve_FEMVMSStabSteadyStateNSE ...
        (fldMsh_p1, up, homDOFs, inhomDOFs, valuesInhomDBCModified, ...
        uMeshALE, propParameters, computeBodyForces, propAnalysis, ...
        solve_LinearSystem, propFldDynamics, propNLinearAnalysis, ...
        counterOpt, propGaussInt, propVTK_false, caseName, '');
    
    %% Compute the perturbed drag force
    postProc_update = computePostProc ...
        (FComplete, propAnalysis, propParameters, propPostproc);
    forcesOnDomain = postProc_update.valuePostProc{1};
    drag = forcesOnDomain(1, 1);

    %% Compute sensitivity via finite differencing
    drag_dp1 = (drag - drag_nom)/epsilonTilde;
    
    %% Compute Hessian of the system
    djd1 = djdp1_(drag_dp1);
    
    %% Compute the design update
    dr = -alphaTilde*djd1;
    

    %% Increment iteration counter
    counterOpt = counterOpt + 1;
        
    %% Update progress bar
    fprintf('\b|\n');
end
disp(['Elapsed time: ', num2str(toc)]);

%% Postprocessing

% Clean the optimization history from zero values
optHistory(counterOpt + 1:end,:) = [];

% Plot the history of the objective and the design
yyaxis left
plot(1:counterOpt - 1, optHistory(:, 1));
ylabel('Objective');
yyaxis right
plot(1:counterOpt - 1, optHistory(:, 2));
ylabel('Design');
axis on;
grid on;
xlabel('Optimization iteration');

%% END OF THE SCRIPT