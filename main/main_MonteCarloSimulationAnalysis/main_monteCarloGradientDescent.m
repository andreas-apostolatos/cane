%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos, Matthew Keller
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
caseName = '2D_Building_validation_test';

%% Parse the data from the GiD input file
[fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,propALE,~,analysis,parameters,...
    propNLinearAnalysis,propFldDynamics,gaussInt,postProc] = ...
    parse_FluidModelFromGid...
    (pathToCase,caseName,'');
propNLinearAnalysis.maxIter = 10;
propNLinearAnalysis.eps = 1e-5;

%% GUI
% Request user defined parameter input
[samplingCall,Umax,iterationLimit,design_penalization,learning_rate,propALE.propUser.delta_p1, ...
    propALE.propUser.delta_p2,propALE.propUser.delta_p3,h_limit,w_limit] = MonteCarloInputGUI();

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
valuesInhomDBCModified = computeInletVelocityPowerLaw(fldMsh, inhomDOFs, valuesInhomDOFs, Umax);

%% Variable initialization
i = 1; % Counter initialization for iteration tree search

% Initialize boundary values of structure
[propALE.propUser.x_Base_Min,propALE.propUser.x_Mid,propALE.propUser.x_Base_Width,...
  propALE.propUser.x_Top_Width,propALE.propUser.y_Max] = computeStructureBoundary(fldMsh,propALE);

% Initialize gamma
gamma = 1e-4;

% Initial accuracy containers
djd1 = 1; % Initial accuracy of parameter 1
djd2 = 1; % Initial accuracy of parameter 2
djd3 = 1; % Initial accuracy of parameter 3
djdp1 = [];
djdp2 = [];
djdp3 = [];
p1_hist = [];
p2_hist = [];
p3_hist = [];

% Initialize iteration steps
propALE.propUser.iterate_p1 = 0;
propALE.propUser.iterate_p2 = 0;
propALE.propUser.iterate_p3 = 0;

switch samplingCall
    case 'Height'
        p1_0 = propALE.propUser.y_Max; % Initialize state values
        p1 = p1_0; % Initialize parameter states
        p2 = '';
        p3 = '';
        p_IterateCall = [1,0,0]; % Define perturbed parameters: [height,width,taper]
        u_flag = 1; % Define mesh adjustment call
    case 'Width'
        p2_0 = propALE.propUser.x_Base_Width;
        p2 = p2_0;
        p1 = '';
        p3 = '';
        p_IterateCall = [0,1,0];
        u_flag = 2;
    case 'Height and Width'
        p1_0 = propALE.propUser.y_Max;
        p2_0 = propALE.propUser.x_Base_Width;
        p1 = p1_0;
        p2 = p2_0;
        p3 = '';
        p_IterateCall = [1,1,0];
        u_flag = 3;
    case 'Taper'       
        p_Base_0 = propALE.propUser.x_Base_Width;
        p_Top_0 = propALE.propUser.x_Top_Width;
        p3_0 = p_Base_0 - p_Top_0;
        p3 = p3_0;
        p1 = '';
        p2 = '';
        p_IterateCall = [0,0,1];
        u_flag = 4;
    case 'Height and Taper'
        p1_0 = propALE.propUser.y_Max;
        p_Base_0 = propALE.propUser.x_Base_Width;
        p_Top_0 = propALE.propUser.x_Top_Width;
        p3_0 = p_Base_0 - p_Top_0;
        p1 = p1_0;
        p3 = p3_0;
        p2 = '';
        p_IterateCall = [1,0,1];
        u_flag = 5;
end

% Initialization of the solution
noNodes = length(fldMsh.nodes(:,1));
noDOFs = 3*noNodes;
up = zeros(noDOFs,1);

%% Define abstract base functions and logging process
djdp1_ = @(djdp1_, p1) djdp1_ + design_penalization*(p1 - p1_0);
djdp2_ = @(djdp2_, p2) djdp2_ + design_penalization*(p2 - p2_0);
djdp3_ = @(djdp3_, p3) djdp3_ + design_penalization*(p3 - p3_0);

% Set up progress bar
fprintf(['\n' repmat('.',1,iterationLimit) '\n\n']);

%Start time count
tic

%% Main loop to solve CFD problem for each Monte Carlo random sampling and optimization processes
while (max([abs(djd1), abs(djd2), abs(djd3)]) > 1e-4 && i <= iterationLimit &&  p1 >= h_limit &&  p2 >= w_limit)    
    
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
        
        % Compute gradient for parameter 1
        djd1 = djdp1_(drag_dp1,propALE.propUser.p1);      
    end
    
    %% Solve the CFD problem with perturbed p2
    if (p_IterateCall(2) == 1)
        
        p_flag = 2;
        
        [p2,lift,drag] = perturbedSolution(fldMsh,up,homDOFs,inhomDOFs,valuesInhomDBCModified,...
            parameters,computeBodyForces,analysis,computeInitialConditions,...
            VTKResultFile,valuesInhomDOFs,propALE,solve_LinearSystem,propFldDynamics,...
            propNLinearAnalysis,i,propVTK_false,gaussInt,caseName,p_flag,postProc);

        output_p2(i,:) = [i,p2,lift,drag];
        
        drag_dp2 = (drag - referenceDrag_Nom) / propALE.propUser.delta_p2;
        
        djd2 = djdp2_(drag_dp2,propALE.propUser.p2);
    end
    
    %% Solve the CFD problem with perturbed p3
    if (p_IterateCall(3) == 1)

        p_flag = 3;
        
        [p3,lift,drag] = perturbedSolution(fldMsh,up,homDOFs,inhomDOFs,valuesInhomDBCModified,...
            parameters,computeBodyForces,analysis,computeInitialConditions,...
            VTKResultFile,valuesInhomDOFs,propALE,solve_LinearSystem,propFldDynamics,...
            propNLinearAnalysis,i,propVTK_false,gaussInt,caseName,p_flag,postProc);

        output_p3(i,:) = [i,p3,lift,drag];
        
        drag_dp3 = (drag - referenceDrag_Nom) / propALE.propUser.delta_p3;
        
        djd3 = djdp3_(drag_dp3,propALE.propUser.p3);
    end
    
    %% Update step size - Barzilai-Borwein step length for gradient descent   
    if isempty(p1_hist) == 0
        [denom,gamma] = gradientDescent(propALE,p1_hist,p2_hist,p3_hist,...
                        djd1,djd2,djd3,djdp1,djdp2,djdp3,u_flag);
    end
    
    %% Update parameter values for next iteration
    if (p_IterateCall(1) == 1)
        propALE.propUser.iterate_p1 = sign(djd1)*min(abs(gamma*djd1), (0.01*p1_0)); % Calculate gradient descent
        p1_hist(i,:) = propALE.propUser.p1; % Update parameter history
        djdp1(i,:) = djd1;
    end
    if (p_IterateCall(2) == 1)
        propALE.propUser.iterate_p2 = sign(djd2)*min(abs(gamma*djd2), (0.01*p2_0));
        p2_hist(i,:) = propALE.propUser.p2;
        djdp2(i,:) = djd2; 
    end
    if (p_IterateCall(3) == 1)
        propALE.propUser.iterate_p3 = sign(djd3)*min(abs(gamma*djd3), (0.01*p3_0));
        p3_hist(i,:) = propALE.propUser.p3;
        djdp3(i,:) = djd3; 
    end
    
    % Gradient descent update
    p1 = propALE.propUser.p1 - propALE.propUser.iterate_p1;
    p2 = propALE.propUser.p2 -  propALE.propUser.iterate_p2;
    p3 = propALE.propUser.p3 -  propALE.propUser.iterate_p3;
    
    %% Update the mesh for the adjusted nominal state    
    fldMsh = correctionSolution(fldMsh,homDOFs,inhomDOFs,...
                valuesInhomDOFs,propALE,solve_LinearSystem,...
                propFldDynamics,i,u_flag);  
    
    %% Increment iteration counter
    i = i+1;
        
    %% Update progress bar
    fprintf('\b|\n');
        
end

%End time count
disp(['Elapsed time: ', num2str(toc)])

%% END OF THE SCRIPT