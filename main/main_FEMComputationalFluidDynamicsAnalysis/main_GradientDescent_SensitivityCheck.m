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
% Date : 11.3.2020
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
        '../../FEMComputationalFluidDynamicsAnalysis/output/',...
        '../../FEMComputationalFluidDynamicsAnalysis/ALEMotion/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');
addpath('../../FEMComputationalFluidDynamicsAnalysis/postProcessing/');

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_Cylinder2D_SensitivityCheck';

%% Parse the data from the GiD input file
[fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,~,analysis,parameters,...
    propNLinearAnalysis,propFldDynamics,gaussInt,postProc] = ...
    parse_FluidModelFromGid...
    (pathToCase,caseName,'');
propNLinearAnalysis.maxIter = 10;
propNLinearAnalysis.eps = 1e-5;

%% GUI
% On the body forces
computeBodyForces = @computeConstantVerticalBodyForceVct;

% On the initial conditions
computeInitialConditions = @computeNullInitialConditionsFEM4NSE2D;

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
VTKResultFile = 'undefined';
propVTK_true.isOutput = true;
propVTK_false.isOutput = false;

%% Input parameters
% max input velocity defined in the reference paper
Umax = .1;

%% Change input velocity to have the parabolic distribution for each randomized input
valuesInhomDBCModified = computeInletVelocityPowerLaw(fldMsh, inhomDOFs, valuesInhomDOFs, Umax);

%% Variable initialization
i = 1; % Counter initialization for iteration search
iterationLimit = 15; % Limit search iterations
design_penalization = 1e4; % Design penalty used in gradient calculations

% Initialize abstract function containers
djdp1 = [];
p1_hist = [];

% Initialize accuracy parameters
djd1 = 1; % Initial accuracy of parameter 1
propALE.propUser.delta_p1 = .01; % Perturbation of parameter 1
propALE.propUser.iterate_p1 = 0;
gamma = 1e-4;

%Initialize state values
[x_Min,y_Min,p1_0] = computeStructureBoundary_Cylinder(fldMsh,propALE); % Minimum cylinder center coordinates and initial radius

% Initialize parameter states
p1 = p1_0;
propALE.propUser.x_Mid = x_Min + p1_0;
propALE.propUser.y_Mid = y_Min + p1_0;

% Initialization of the solution
noNodes = length(fldMsh.nodes(:,1));
noDOFs = 3*noNodes;
up = zeros(noDOFs,1);

% Parameters (I/O)
PlotFlag = 'False';

%% Define abstract base functions and logging process
djdp1_ = @(djdp1_, p1) djdp1_;

% Set up progress bar
fprintf(['\n' repmat('.',1,iterationLimit) '\n\n']);

%Start time count
tic

%% Main loop to solve CFD problem for each Monte Carlo random sampling and optimization processes
while (abs(djd1) > 1e-4 && i <= iterationLimit)
    %% Update internal variables with updated values from previous iteration
    propALE.propUser.p1 = p1;
     
    %% Solve the CFD problem in nominal state   
    [up,FComplete,~,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh,up,homDOFs,inhomDOFs,valuesInhomDBCModified,'undefined',parameters,...
        computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
        i,propVTK_true,gaussInt,caseName,'outputEnabled');
    
    % Calculate drag and lift force from the nodal forces
    postProc_update = computePostProc(FComplete, analysis, parameters, postProc);

    % Retrieve Fx and Fy from post processing
    forcesOnDomain = postProc_update.valuePostProc{1};
    Fx = forcesOnDomain(1,1);
    Fy = forcesOnDomain(2,1);

    % Set drag and lift
    lift = Fy;
    drag = Fx;

    % Write nominal output vector
    output_Nom(i,:) = [i, p1, lift, drag];  
    referenceDrag_Nom = drag; % Define reference drag of nominal state for accuracy calculations
 
    %% Solve the CFD problem with perturbed p1   
    p1 = propALE.propUser.p1 + propALE.propUser.delta_p1; % Adjust user-defined parameter 1
    propALE.propUser.Perturb_Flag = 'perturb';
   
    % Update the mesh for perturbed state 1
    [fldMsh_p1,~,~,~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D...
        (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,...
        solve_LinearSystem,propFldDynamics, i);
    
    fldMsh_p1.initialNodes = fldMsh_p1.nodes;
    
    [~,FComplete,~,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh_p1,up,homDOFs,inhomDOFs,valuesInhomDBCModified,propALE,parameters,...
        computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
        i,propVTK_false,gaussInt,caseName,'');

    postProc_update = computePostProc(FComplete, analysis, parameters, postProc);

    forcesOnDomain = postProc_update.valuePostProc{1};
    Fx = forcesOnDomain(1,1);
    Fy = forcesOnDomain(2,1);

    lift = Fy;
    drag = Fx;

    % Write output vector
    output_p1(i,:) = [i, p1, lift, drag];  

    % Compute sensitivities via finite differencing
    drag_dp1 = (drag - referenceDrag_Nom) / propALE.propUser.delta_p1;
              
    %% Compute gradient
    djd1 = djdp1_(drag_dp1,propALE.propUser.p1); % Compute gradient for parameter 1

    propALE.propUser.iterate_p1 = djd1*propALE.propUser.delta_p1;
    p1 = propALE.propUser.p1 - propALE.propUser.iterate_p1;
    
    propALE.propUser.Perturb_Flag = 'set_final';

    %% Update the mesh for the adjusted nominal state  
    [fldMsh,~,~,~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D...
        (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,...
        solve_LinearSystem,propFldDynamics, i);
    
    fldMsh.initialNodes = fldMsh.nodes;

    %% Increment iteration counter
    i = i + 1;
        
    %% Update progress bar
    fprintf('\b|\n');
        
end

%End time count
disp(['Elapsed time: ', num2str(toc)])
%% END OF THE SCRIPT