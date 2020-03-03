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
%        in flow whose design parameters are its width and height. The
%        sensitivities are computed using Finite Differencing whereas a
%        gradient-based descend optimization algorithm is employed.
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
caseName = '2D_Building_CFD_optimization_validation_test';

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
propVTK_true.isOutput = true;
propVTK_false.isOutput = false;

%% Input parameters

% max input velocity defined in the reference paper
Umax = .1;

% define parameters used in reference paper and simualiton
D = 0.1;    % diameter of the body
Ubar = 0.2; % mid velocity
rho = parameters.rho; % density

%% Change input velocity to have the parabolic distribution for each randomized input
% valuesInhomDBCModified = computeInletVelocityParabolic_unitTest(fldMsh, inhomDOFs, valuesInhomDOFs, Umax);
valuesInhomDBCModified = computeInletVelocityPowerLaw(fldMsh, inhomDOFs, valuesInhomDOFs, Umax);

%% Variable initialization

i = 1; % Counter initialization for iteration tree search
iterationLimit = 20; % Limit search iterations
design_penalization = 1e4; % Design penalty used in J calculations

% Initialize abstract function containers
djdp1 = [];
djdp2 = [];
p1_hist = [];
p2_hist = [];

% Initialize accuracy parameters
djd1 = 1; % Initial accuracy of parameter 1
djd2 = 1; % Initial accuracy of parameter 2
propALE.propUser.delta_p1 = .001; %.1, 0.005; % Perturbation of parameter 1
propALE.propUser.delta_p2 = .05;%,.05, 0.001; % Perturbation of parameter 2 = 0 to analyze p1 only
propALE.propUser.iterate_p1 = 0;
propALE.propUser.iterate_p2 = 0;
gamma = 1e-4;

%Initialize state values
p1_0 = 0.1; % Initial height of 2d building from GID input file
p2_0 = 0.02; % Initial width of 2d building from GID input file
x_0 = computeStructureBoundary(fldMsh,propALE); % Initial x location of 2d building from GID input file

% Initialize parameter states
p1 = p1_0;
p2 = p2_0;
propALE.propUser.x_Mid = x_0 + (p2*0.5); % Locate center of building for dx motion

% Initialization of the solution
noNodes = length(fldMsh.nodes(:,1));
noDOFs = 3*noNodes;
up = zeros(noDOFs,1);

% Parameters (I/O)
PlotFlag = 'False';

%% Define abstract base functions and logging process
%J_  = @(J_, p1, p2) J_ + design_penalization*0.5*((p1 - p1_0)^2.0 + (p2 - p2_0)^2.0);
% djdp1_ = @(djdp1_, p1) djdp1_ + design_penalization*(p1 - p1_0);
djdp1_ = @(djdp1_, p1) djdp1_;
% djdp2_ = @(djdp2_, p2) djdp2_ + design_penalization*(p2 - p2_0);

% Set up progress bar
fprintf(['\n' repmat('.',1,iterationLimit) '\n\n']);

%Start time count
tic

%% Main loop to solve CFD problem for each Monte Carlo random sampling and optimization processes
while (max(abs(djd1),abs(djd2)) > 1e-4 && i <= iterationLimit)    
    %% Update internal variables with updated values from previous iteration
    propALE.propUser.p1 = p1;
    propALE.propUser.p2 = p2;
     
    %% Solve the CFD problem in nominal state   
    [up,FComplete,~,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh,up,homDOFs,inhomDOFs,valuesInhomDBCModified,'undefined',parameters,...
        computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
        i,propVTK_true,gaussInt,caseName,'outputEnabled');
    
    % Calculate drag and lift force from the nodal forces
    postProc_update = computePostProc(FComplete, analysis, parameters, postProc);

    % Calculate drag and lift coefficient based on drag and lift force
    % Retrieve Fx and Fy from post processing
    forcesOnDomain = postProc_update.valuePostProc{1};
    Fx = forcesOnDomain(1,1);
    Fy = forcesOnDomain(2,1);

    % Calculate drag and lift coefficients
%     dragCoefficient = (2 * Fx)/(rho * Ubar * Ubar * D);
%     liftCoefficient = (2 * Fy)/(rho * Ubar * Ubar * D);
    lift = Fy;
    drag = Fx;

    % Remove negative coefficients
%     liftCoefficient = abs(liftCoefficient);
    liftCoefficient = abs(lift);

    % Write nominal output vector
%     output_Nom(i,:) = [i, liftCoefficient, dragCoefficient];  
    output_Nom(i,:) = [i, lift, drag];  
%     referenceDrag_Nom = dragCoefficient; % Define reference drag of nominal state for accuracy calculations
    referenceDrag_Nom = drag; % Define reference drag of nominal state for accuracy calculations
 
    %% Solve the CFD problem with perturbed p1   
    p1 = propALE.propUser.p1 + propALE.propUser.delta_p1; % Adjust user-defined parameter 1
    propALE.propUser.Perturb_Flag = 'dy';
    
    % Update the mesh for perturbed state 1
    [fldMsh_p1,~,~,~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D...
        (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,...
        solve_LinearSystem,propFldDynamics, i);

    [~,FComplete,~,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh_p1,up,homDOFs,inhomDOFs,valuesInhomDBCModified,propALE,parameters,...
        computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
        i,propVTK_false,gaussInt,caseName,'');

    postProc_update = computePostProc(FComplete, analysis, parameters, postProc);

    forcesOnDomain = postProc_update.valuePostProc{1};
    Fx = forcesOnDomain(1,1);
    Fy = forcesOnDomain(2,1);

%     dragCoefficient = (2 * Fx)/(rho * Ubar * Ubar * D);
%     liftCoefficient = (2 * Fy)/(rho * Ubar * Ubar * D);
    lift = Fy;
    drag = Fx;

%     liftCoefficient = abs(liftCoefficient);
    liftCoefficient = abs(lift);

    % Write output vector
%     output_p1(i,:) = [p1, liftCoefficient, dragCoefficient];
    output_p1(i,:) = [p1, lift, drag];  

    % Compute sensitivities via finite differencing
%     drag_dp1 = (dragCoefficient - referenceDrag_Nom) / propALE.propUser.delta_p1;
    drag_dp1 = (drag - referenceDrag_Nom) / propALE.propUser.delta_p1;
          
    %% Solve the CFD problem with perturbed p2  
%     p2 = propALE.propUser.p2 + propALE.propUser.delta_p2; % Adjust user-defined parameter 1
%     propALE.propUser.Perturb_Flag = 'dx';
%     
%     % Update the mesh for perturbed state 2
%     [fldMsh_p2,~,~,~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D...
%         (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,...
%         solve_LinearSystem,propFldDynamics, i);
%     
%     [~,FComplete,~,~] = solve_FEMVMSStabSteadyStateNSE2D...
%         (fldMsh_p2,up,homDOFs,inhomDOFs,valuesInhomDBCModified,propALE,parameters,...
%         computeBodyForces,analysis,computeInitialConditions,...
%         VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
%         i,propVTK_false,gaussInt,caseName,'');
% 
%     postProc_update = computePostProc(FComplete, analysis, parameters, postProc);
% 
%     forcesOnDomain = postProc_update.valuePostProc{1};
%     Fx = forcesOnDomain(1,1);
%     Fy = forcesOnDomain(2,1);
% 
% %     dragCoefficient = (2 * Fx)/(rho * Ubar * Ubar * D);
% %     liftCoefficient = (2 * Fy)/(rho * Ubar * Ubar * D);
%     lift = Fy;
%     drag = Fx;
% 
% %     liftCoefficient = abs(liftCoefficient);
%     liftCoefficient = abs(lift);
% 
% %     output_p2(i,:) = [p2, liftCoefficient, dragCoefficient];
%     output_p2(i,:) = [p2, lift, drag];
%  
% %     drag_dp2 = (dragCoefficient - referenceDrag_Nom) / propALE.propUser.delta_p2;
%     drag_dp2 = (drag - referenceDrag_Nom) / propALE.propUser.delta_p2;
     
    %% Compute objective function and gradients
    if i > 1
        djd1 = djdp1_(drag_dp1,propALE.propUser.p1); % Compute gradient for parameter 1
%         djd2 = djdp2_(drag_dp2,propALE.propUser.p2); % Compute gradient for parameter 2
    end
%    djd2 = 0; % To anlayze p1 only set djd2 = 0      
    
    %% Update step size - Barzilai-Borwein step length for gradient descent
%     if isempty(p1_hist) == 0      
%         denom = (djd1 - djdp1(end))^2 + (djd2 - djdp2(end))^2;
%         gamma = ((propALE.propUser.p1 - p1_hist(end)) * (djd1 - djdp1(end)) + (propALE.propUser.p2 - p2_hist(end)) * (djd2 - djdp2(end)))/denom;
%     end
    
    %% Update parameter values for next iteration
%     propALE.propUser.iterate_p1 = sign(djd1)*min(abs(gamma*djd1), (0.01*p1_0)); % Calculate gradient descent
%     propALE.propUser.iterate_p2 = sign(djd2)*min(abs(gamma*djd2), (0.01*p2_0)); % Calculate gradient descent
%     p1 = propALE.propUser.p1 + propALE.propUser.iterate_p1; % Gradient descent update
%     p2 = propALE.propUser.p2 +  propALE.propUser.iterate_p2; % Gradient descent update
%     p2 = max(p2, (0.5*p2_0)); % Width constraint
    
    propALE.propUser.iterate_p1 = djd1*propALE.propUser.delta_p1;
    p1 = propALE.propUser.p1 + propALE.propUser.iterate_p1;
    
    % Update parameter history
    p1_hist = [p1_hist, propALE.propUser.p1]; 
%     p2_hist = [p2_hist, propALE.propUser.p2];
    djdp1 = [djdp1, djd1];
%     djdp2 = [djdp2, djd2];   
    
    %% Update the mesh for the adjusted nominal state
    propALE.propUser.Perturb_Flag = 'dxdy';
    
    [fldMsh,~,~,~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D...
        (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,...
        solve_LinearSystem,propFldDynamics, i);
    
    %% Increment iteration counter
    i = i + 1;
        
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