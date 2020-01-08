%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos, Matthew Keller
%
%% Script documentation
%
% Task : Performs Multilevel Monte-Carlo for the 2D Navier-Stokes equations
%        in 2D
%
% Date : 29.12.2019
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
addpath('../../FEMComputationalFluidDynamicsAnalysis/postProcessing/');

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_2D_building';

%% Parse the data from the GiD input file
[fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,nodesALE,~,analysis,parameters,...
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

%% Input parameters
% max input velocity defined in the reference paper
Umax = 0.3;

% define parameters used in reference paper and simualiton
D = 0.1;    % diameter of the body
Ubar = 0.2; % mid velocity 

%% Change input velocity to have the parabolic distribution for each randomized input
valuesInhomDBCModified = computeInletVelocityParabolic_unitTest(fldMsh, inhomDOFs, valuesInhomDOFs, Umax);

%% Variable initialization
i = 1; % Counter initialization for iteration tree search
iterationLimit = 10; % Limit search iterations
design_penalization = 1e4; % Design penalty used in J calculations

% Initialize abstract function containers
J = zeros(iterationLimit+1);
J(1) = 0;
dJdp1 = zeros(iterationLimit+1 , 1);
dJdp1(1) = 0;
dJdp2 = zeros(iterationLimit+1 , 1);
dJdp2(1) = 0;
p1_hist = zeros(iterationLimit+1 , 1);
p1_hist(1) = 0;
p2_hist = zeros(iterationLimit+1 , 1);
p2_hist(1) = 0;

% Initialize accuracy parameters
djd1 = 1; % Initial accuracy of parameter 1
djd2 = 1; % Initial accuracy of parameter 2
delta_p1 = 0.01; % Perturbation of parameter 1
delta_p2 = 0; % Perturbation of parameter 2 = 0 to analyze p1 only

%Initialize initial state values
p1_0 = 10; % Initial height of 2d building from GID input file
p2_0 = 2.5; % Initial width of 2d building from GID input file

p1 = p1_0; % Initialize initial states
p2 = p2_0;

% Parameters (I/O)
PlotFlag = 'False';

%% Define abstract base functions and logging process
J_  = @(J_, p1, p2) J_ + design_penalization*0.5*((p1 - p1_0)^2.0 + (p2 - p2_0)^2.0);
dJdp1_ = @(dJdp1_, p1) dJdp1_ + design_penalization*(p1 - p1_0);
dJdp2_ = @(dJdp2_, p2) dJdp2_ + design_penalization*(p2 - p2_0);

% Set up progress bar
fprintf(['\n' repmat('.',1,iterationLimit) '\n\n']);

%Start time count
tic

%% Main loop to solve CFD problem for each Monte Carlo random sampling and optimization processes
while (max(abs(djd1),abs(djd2)) > 1e-4 && i <= iterationLimit)    
    %% Initialize temp variables with updated values from previous iteration
    temp_p1 = p1;
    temp_p2 = p2;
    
    %% Solve the CFD problem in nominal state   
    [~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh,homDOFs,inhomDOFs,valuesInhomDBCModified,nodesALE,parameters,...
        computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
        gaussInt,caseName,'');

    %% Calculate drag and lift force from the nodal forces
    postProc_update = computePostProc(FComplete, analysis, parameters, postProc);

    %% Calculate drag and lift coefficient based on drag and lift force
    rho = parameters.rho; % density

    % get Fx and Fy from post processing
    forcesOnDomain = postProc_update.valuePostProc{1};
    Fx = forcesOnDomain(1,1);
    Fy = forcesOnDomain(2,1);

    % calculate drag and lift coefficiet
    dragCoefficient = (2 * Fx)/(rho * Ubar * Ubar * D);
    liftCoefficient = (2 * Fy)/(rho * Ubar * Ubar * D);

    % find absolute value so we don't get negative coefficients
    liftCoefficient = abs(liftCoefficient);

    %% Write nominal output vector
    output_Nom(i,:) = [i, liftCoefficient, dragCoefficient];  
    referenceDrag_Nom = dragCoefficient;
 
    %% Solve the CFD problem with perturbed p1
    p2 = temp_p2; % Ensure parameter 2 is nominal
    p1 = p1 + delta_p1; % Adjust parameter 1
        
    [~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh,homDOFs,inhomDOFs,valuesInhomDBCModified,nodesALE,parameters,...
        computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
        gaussInt,caseName,'');

    %% Calculate drag and lift force from the nodal forces
    postProc_update = computePostProc(FComplete, analysis, parameters, postProc);

    %% Calculate drag and lift coefficient based on drag and lift force
    rho = parameters.rho; % density

    % get Fx and Fy from post processing
    forcesOnDomain = postProc_update.valuePostProc{1};
    Fx = forcesOnDomain(1,1);
    Fy = forcesOnDomain(2,1);

    % calculate drag and lift coefficiet
    dragCoefficient = (2 * Fx)/(rho * Ubar * Ubar * D);
    liftCoefficient = (2 * Fy)/(rho * Ubar * Ubar * D);

    % find absolute value so we don't get negative coefficients
    liftCoefficient = abs(liftCoefficient);

    % Write output vector
    output_p1(i,:) = [p1, liftCoefficient, dragCoefficient];

    % Compute sensitivities via finite differences
    drag_dp1 = (dragCoefficient - referenceDrag_Nom) / delta_p1;
        
%     %% Solve the CFD problem with perturbed p2  
%     p1 = temp_p1; % Ensure parameter 1 is nominal
%     p2 = p2 + delta_p2; % Adjust parameter 2
%       
%     [~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
%         (fldMsh,homDBC,inhomDBC,valuesInhomDBCModified,nodesALE,parameters,...
%         computeBodyForces,analysis,computeInitialConditions,...
%         VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
%         gaussInt,caseName,'');
% 
%     %% Calculate drag and lift force from the nodal forces
%     postProc_update = computePostProc(FComplete, analysis, parameters, postProc);
% 
%     %% Calculate drag and lift coefficient based on drag and lift force
%     rho = parameters.rho; % density
% 
%     % get Fx and Fy from post processing
%     forcesOnDomain = postProc_update.valuePostProc{1};
%     Fx = forcesOnDomain(1,1);
%     Fy = forcesOnDomain(2,1);
% 
%     % calculate drag and lift coefficient
%     dragCoefficient = (2 * Fx)/(rho * Ubar * Ubar * D);
%     liftCoefficient = (2 * Fy)/(rho * Ubar * Ubar * D);
% 
%     % find absolute value so we don't get negative coefficients
%     liftCoefficient = abs(liftCoefficient);
% 
%     % Write output vector
%     output_p2(i,:) = [p2, liftCoefficient, dragCoefficient];
%  
%     % Compute sensitivities via finite differences
%     drag_dp2 = (dragCoefficient - referenceDrag_Nom) / delta_p2;
     
    %% Compute objective function and gradients
    J(i+1) = J_(referenceDrag_Nom,temp_p1,temp_p2); % Compute J
    
    dJdp1(i+1) = dJdp1_(drag_dp1,temp_p1); % Compute gradient for parameter 1
    djd1 = dJdp1(i+1); % Populate accuracy parameter 1
    p1_hist(i+1) = temp_p1; % Record p1 values
    
    dJdp2(i+1) = 0; % To anlayze p1 only set djdp2 = 0   
%     dJdp2(i+1) = dJdp2_(drag_dp2,temp_p2); % Compute gradient for parameter 2
%     djd2 = dJdp2(i+1); % Populate accuracy parameter 2
    p2_hist(i+1) = temp_p2; % Record p2 values
    
    %% Update step size - Barzilai-Borwein step length for gradient descent
    denom = (dJdp1(i+1) - dJdp1(i))^2 + (dJdp2(i+1) - dJdp2(i))^2;
    gamma = ((p1_hist(i+1) - p1_hist(i)) * (dJdp1(i+1) - dJdp1(i)) + (p2_hist(i+1) - p2_hist(i)) * (dJdp2(i+1) - dJdp2(i)))/denom;
    
    %% Update parameter values for next iteration
    p1 = temp_p1 - sign(dJdp1(i+1))*min(abs(gamma*dJdp1(i+1)), 0.1); 
%     p2 = temp_p2 - sign(dJdp2(i+1))*min(abs(gamma*dJdp2(i+1)), 0.1);

    %% Update the mesh
    %
    % Comments: The function below returns the updated mesh given the
    %           boundary conditions defined in nodesALE via the function
    %           handles for each of the nodes.
    % 
    [fldMsh,~,~,~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D...
        (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,nodesALE,...
        solve_LinearSystem,propFldDynamics,0);
        
    %% Increment iteration counter
    i = i+1;
        
    %% Update progress bar
    fprintf('\b|\n');
        
end

%End time count
disp(['Elapsed time: ', num2str(toc)])

%% Plot output results
% Drag coefficient
if strcmp(PlotFlag, 'True')
    figure(graph.index)
    histfit(output(:,2))
    yt = get(gca, 'YTick');
    set(gca, 'YTick', yt, 'YTickLabel', yt/numel(output(:,2)))
    title('Output histogram for the drag coefficient')
    graph.index = graph.index + 1;

    % Lift coefficient
    figure(graph.index)
    histfit(output(:,3))
    yt = get(gca, 'YTick');
    set(gca, 'YTick', yt, 'YTickLabel', yt/numel(output(:,3)))
    title('Output histogram for the lift coefficient')
    graph.index = graph.index + 1;
end
%% END OF THE SCRIPT