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
% Task : Performs steady-state fluid optimization for a rectangular tower
%        in flow whose design parameters depend on the defined geometry. The
%        sensitivities are computed using Finite Differencing whereas a
%        gradient-based descend optimization algorithm is employed.
%
% Date : 22.03.2020
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
addpath('../../FEMComputationalFluidDynamicsAnalysis/SolutionStructure/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');
addpath('../../FEMComputationalFluidDynamicsAnalysis/postProcessing/');

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = '2D_Building_CFD_optimization_validation_test';

%% Parse the data from the GiD input file
[fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,~,analysis,parameters,...
    propNLinearAnalysis,propFldDynamics,gaussInt,postProc] = ...
    parse_FluidModelFromGid...
    (pathToCase,caseName,'');
propNLinearAnalysis.maxIter = 10;
propNLinearAnalysis.eps = 1e-5;

%% GUI - ONLY WORKS WITH HEIGHT
% Request user defined parameter input
analysis_type = 'simplified';

[samplingCall,Umax,iterationLimit,design_penalization,propALE.propUser.delta_p1,...
 propALE.propUser.delta_p2,propALE.propUser.delta_p3,convergenceLimit,...
 h_limit,w_limit,t_limit] = MonteCarloInputGUI(analysis_type);

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

%% Change input velocity to have the power law distribution for each randomized input
valuesInhomDBCModified = computeInletVelocityPowerLaw(fldMsh,inhomDOFs,valuesInhomDOFs,Umax);

%% Variable initialization
i = 1; % Counter initialization for iteration tree search

limit_flag = false; % Initialize limit container

% Initialize boundary values of structure
[propALE.propUser.x_Base_Min,propALE.propUser.x_Mid,propALE.propUser.x_Base_Width,...
 propALE.propUser.x_Top_Width,propALE.propUser.y_Max] = computeStructureBoundary(fldMsh,propALE);

% Initial accuracy containers
djd1 = 1; % Initialize accuracy of parameter 1
djdp1 = [];
p1_hist = [];

% Initialize iteration steps
propALE.propUser.iterate_p1 = 0;

switch samplingCall
    case 'Height'
        p1_0 = propALE.propUser.y_Max; % Initialize state values
        p1 = p1_0; % Initialize parameter states
        p2 = '';
        p3 = '';
        p_IterateCall = [1,0,0]; % Define perturbed parameters: [height,width,taper]
        u_flag = 1; % Define mesh adjustment call
end

% Initialization of the solution
noNodes = length(fldMsh.nodes(:,1));
noDOFs = 3*noNodes;
up = zeros(noDOFs,1);

% Set up progress bar
fprintf(['\n' repmat('.',1,iterationLimit) '\n\n']);

%Start time count
tic

%% Main loop to solve CFD problem for each Monte Carlo random sampling and optimization processes
while (abs(djd1) > convergenceLimit && i <= iterationLimit && limit_flag == false)
    %% Check limit conditions
    limit_flag = limitCalculation(p1,p2,p3,h_limit,w_limit,t_limit,u_flag);
    
    %% Solve the CFD problem at the nominal (updated) condition
    [propALE,lift,drag] = nominalSolution(fldMsh,up,homDOFs,inhomDOFs,valuesInhomDBCModified,...
        parameters,computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,propALE,solve_LinearSystem,propFldDynamics,...
        propNLinearAnalysis,i,propVTK_true,gaussInt,caseName,p1,p2,p3,postProc);
    
    % Write nominal output vector
    output_Nom(i,:) = [i,p1,lift,drag];
    
    % Define reference drag of nominal state for accuracy calculations
    referenceDrag_Nom = drag;
    
    %% Solve the CFD problem with perturbed p1
    if (p_IterateCall(1) == 1)
        
        % Set caculation flag
        p_flag = 1;
        
        % Solve perturbed solution
        [p1,lift,drag] = perturbedSolution(fldMsh,up,homDOFs,inhomDOFs,valuesInhomDBCModified,...
            parameters,computeBodyForces,analysis,computeInitialConditions,...
            VTKResultFile,valuesInhomDOFs,propALE,solve_LinearSystem,propFldDynamics,...
            propNLinearAnalysis,i,propVTK_false,gaussInt,caseName,p_flag,postProc);
        
        % Write output vector
        output_p1(i,:) = [i,p1,lift,drag];
        
        % Compute sensitivities via finite differencing
        drag_dp1 = (drag - referenceDrag_Nom) / propALE.propUser.delta_p1;       
    end
    
    % Basic Gradient Descent - https://builtin.com/data-science/gradient-descent
    djd1 = drag_dp1;
    propALE.propUser.iterate_p1 = -(learning_rate*djd1);
    p1 = propALE.propUser.p1 + propALE.propUser.iterate_p1;
    
    % Update parameter history
    p1_hist(i,:) = propALE.propUser.p1; 
    djdp1(i,:) = djd1;
    
    %% Update the mesh for the adjusted nominal state    
    fldMsh = correctionSolution(fldMsh,homDOFs,inhomDOFs,...
                valuesInhomDOFs, propALE, solve_LinearSystem,...
                propFldDynamics,i,u_flag);
               
    %% Increment iteration counter
    i = i+1;
            
    %% Update progress bar
    fprintf('\b|\n');
        
end

%End time count
disp(['Elapsed time: ', num2str(toc)])

%% END OF THE SCRIPT