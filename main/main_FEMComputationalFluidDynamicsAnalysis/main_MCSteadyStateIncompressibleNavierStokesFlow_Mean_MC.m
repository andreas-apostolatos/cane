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

%Add parabolic input function to parallel pool
poolobj = gcp;
poolobj.addAttachedFiles(fullfile(pwd, '../../FEMComputationalFluidDynamicsAnalysis/postProcessing/'));

% Define the path to the case
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

%% Input distribution definitions 
% Request user sample call selection
sample_list = {'randomUniform','randomNormal','latinHyperCube', 'quasiMonteCarloHalton'};
[indx,tf] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',sample_list);
switch indx
    case 1
        samplingCall = 'randomUniform';   
    case 2
        samplingCall = 'randomNormal';   
    case 3
        samplingCall = 'latinHyperCube';   
    case 4
        samplingCall = 'quasiMonteCarloHalton';
end

% User definition of distribution parameters for Monte Carlo Selection
prompt = {'Enter P1 Mean: ','Enter P1 Standard Deviation: ', 'Enter P2 Mean: ','Enter P2 Standard Deviation: ', 'Enter number of samples: '};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0.1','0.02','0.2','0.04','10'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
mean_value_p1 = str2double(answer(1));
standard_deviation_p1 = str2double(answer(2));
mean_value_p2 = str2double(answer(3));
standard_deviation_p2 = str2double(answer(4));
nSample = str2double(answer(5));

%% Generate Monte Carlo input vector
% Set input vector depending on distribution type and number of samples
switch samplingCall
    case 'randomUniform'
        % Uniform distribution with bounds
        distributionFunction_p1 = @(varargin) 2*standard_deviation_p1*rand(nSample,1)+mean_value_p1-standard_deviation_p1;
        distributionFunction_p2 = @(varargin) 2*standard_deviation_p2*rand(nSample,1)+mean_value_p2-standard_deviation_p2;
%     case 'randomNormal'
%         % Random normal distribution sampling
%         distributionFunction_p1 = @(varargin) normrnd(mean_value_p1, standard_deviation_p1, nSample, 1);    
%         distributionFunction_p2 = @(varargin) normrnd(mean_value_p2, standard_deviation_p2, nSample, 1); 
%     case 'latinHyperCube'
%         % Latin hypercube sampling
%         distributionFunction_p1 =  @(varargin) lhsnorm(mean_value_p1, (standard_deviation_p1*standard_deviation_p1), nSample); 
%         distributionFunction_p2 =  @(varargin) lhsnorm(mean_value_p2, (standard_deviation_p2*standard_deviation_p2), nSample);
%     case 'quasiMonteCarloHalton'
%         % Quasi Monte Carlo Halton Sequence sampling
%         p = haltonset(1,'Skip',1e3,'Leap',1e2);
%         p = scramble(p,'RR2');
%         haltonvector = net(p,nSample); % Halton Sequence sampling
%         distributionFunction_p1 = @(varargin) norminv(haltonvector, mean_value_p1, standard_deviation_p1); % Invert uniform Halton sequence to normal distribution
%         p = haltonset(1,'Skip',1e3,'Leap',1e2);
%         p = scramble(p,'RR2');
%         haltonvector = net(p,nSample); % Halton Sequence sampling
%         distributionFunction_p2 = @(varargin) norminv(haltonvector, mean_value_p2, standard_deviation_p2); % Invert uniform Halton sequence to normal distribution
%     otherwise
%         % Error
%         error('Invalid input sampling method')
end

% Generate initial inputs
input_p1 = distributionFunction_p1();
input_p2 = distributionFunction_p2();

% Initialize graph index
graph.index = 1;

% Initialize optimization parameters
design_penalization = 1e4;
% dp1 = 1e-3; % perturbation in parameter 1
% dp2 = 1e-3; % perturbation in parameter 2

% Parameters (I/O)
PlotFlag = 'False';

% Generate histrogram of inputs
if strcmp(PlotFlag, 'True')
    figure(graph.index)
    histfit(input)
    yt = get(gca, 'YTick');
    set(gca, 'YTick', yt, 'YTickLabel', yt/numel(input))
    title('Input histogram')
    graph.index = graph.index + 1;
end

%% Define abstract base functions
J_  = @(J_, p1, p2) J_ + design_penalization*0.5*((p1 - p1_0)*2.0 + (p2 - p2_0)*2.0);
dJdp1_ = @(dJdp1_, p1) dJdp1_ + design_penalization*(p1 - p1_0);
dJdp2_ = @(dJdp2_, p2) dJdp2_ + design_penalization*(p2 - p2_0);

%% Logging processes
% % Set up parallel processing - No while loop applied
if max(size(gcp)) == 0                  % parallel pool needed
    parpool                             % create the parallel pool
end

%% Input parameters
% max input velocity defined in the reference paper
Umax = 0.3;
% define parameters used in reference paper and simualiton
D_0 = 0.1;    % diameter of the body
Ubar_0 = 0.2; % mid velocity 


%% Change input velocity to have the parabolic distribution for each randomized input
valuesInhomDBCModified = computeInletVelocityParabolic_unitTest(fldMsh, inhomDBC, valuesInhomDBC, Umax);

%% Variable initialization
i = 1; % Counter initialization for iteration tree search
iterationLimit = 3; %Limit to search depth of monte carlo tree search
dragCoefficient_hist = zeros(nSample);

% Initialize abstract function containers
J = zeros(iterationLimit+1);
J(1) = 0;
dJdp1 = zeros(iterationLimit+1);
dJdp1(1) = 0;
dJdp2 = zeros(iterationLimit+1);
dJdp2(1) = 0;
p1_hist = zeros(iterationLimit+1);
p1_hist(1) = 0;
p2_hist = zeros(iterationLimit+1);
p2_hist(1) = 0;

% Initialize accuracy parameters
djd1 = 1;
djd2 = 1;

%Initialize temp parameters
p1_bestCase = mean_value_p1;
p2_bestCase = mean_value_p2;
temp_d_drag_dp1 = 0;
temp_d_drag_dp2 = 0;
temp_p1 = D_0;
temp_p2 = Ubar_0;

%Initialize initial values
p1_0 = D_0;
p2_0 = Ubar_0;
D = D_0;
Ubar = Ubar_0;

% Set up progress bar
fprintf(['\n' repmat('.',1,iterationLimit) '\n\n']);

%Start time count
tic

%% Main loop to solve CFD problem for each Monte Carlo random sampling and optimization processes
while (max(abs(djd1),abs(djd2)) > 1e-4 && i <= iterationLimit)    
    %% Solve the CFD problem in nominal state
   
    [~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh,homDBC,inhomDBC,valuesInhomDBCModified,nodesALE,parameters,...
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

    %% Write output
    output(i,:) = [i, liftCoefficient, dragCoefficient];       
    referenceDrag_Nom = dragCoefficient;

    
    %% Solve the CFD problem with shifted rho
    Ubar = temp_p2;
    
    parfor k=1:nSample
        %Build new struct for parameters use inside parallel loop
        D = input_p1(k);
        
        [~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
            (fldMsh,homDBC,inhomDBC,valuesInhomDBCModified,nodesALE,parameters,...
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
        
        output_p1(k,:) = [input_p1(k), liftCoefficient, dragCoefficient];
        
    end
    
    %Determine location of best case drag from MC search
    [d_drag_dp1,index] = min(output_p1(:,3) - referenceDrag_Nom);
    referenceDrag_p1 = output_p1(index,3);
    
    %Update best case p1 if drag is reduced
    if (d_drag_dp1 < 0)       
        p1_bestCase = output_p1(index,1);
    end
    
    % Compute sensitivities via finite differences
    dp1 = abs(referenceDrag_p1-temp_d_drag_dp1);
    drag_dp1 = (referenceDrag_p1 - referenceDrag_Nom) / dp1;
    temp_d_drag_dp1 = referenceDrag_p1;
        
    %% Solve the CFD problem with shifted nue   
    D = temp_p1;
    
    parfor k=1:nSample
        %Build new struct for parameters use inside parallel loop
        Ubar = input_p2(k);
        
        [~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
            (fldMsh,homDBC,inhomDBC,valuesInhomDBCModified,nodesALE,parameters,...
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
        
        output_p2(k,:) = [input_p2(k), liftCoefficient, dragCoefficient];
        
    end
    
    %Determine location of best case drag from MC search
    [d_drag_dp2,index] = min(output_p2(:,3) - referenceDrag_Nom);
    referenceDrag_p2 = output_p2(index,3);
    
    %Update best case p1 if drag is reduced
    if (d_drag_dp2 < 0)       
        p2_bestCase = output_p2(index,1);
    end
    
    % Compute sensitivities via finite differences
    dp2 = abs(referenceDrag_p2-temp_d_drag_dp2);
    drag_dp2 = (referenceDrag_p2 - referenceDrag_Nom) / dp2;
    temp_d_drag_dp2 = referenceDrag_p2;
    
    %% Check the CFD problem in combined state
    D = p1_bestCase;
    Ubar = p2_bestCase;
    
    [~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
        (fldMsh,homDBC,inhomDBC,valuesInhomDBCModified,nodesALE,parameters,...
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
    
    %% Update distribution depending on best case
    if (dragCoefficient < referenceDrag_Nom && dragCoefficient < referenceDrag_p1 &&  dragCoefficient < referenceDrag_p2)
        %Refine distribution region for both parameters
            standard_deviation_p1 = (((mean_value_p1+standard_deviation_p1) - (mean_value_p1-standard_deviation_p1)) - (2.0 * abs(mean_value_p1 - p1_bestCase)))/2.0;
            mean_value_p1 = p1_bestCase;

            standard_deviation_p2 = (((mean_value_p2+standard_deviation_p2) - (mean_value_p2-standard_deviation_p2)) - (2.0 * abs(mean_value_p2 - p2_bestCase)))/2.0;
            mean_value_p2 = p2_bestCase;
            
        %Update both parameters
%             D = p1_bestCase;
%             Ubar = p2_bestCase; 
            
    elseif (dragCoefficient < referenceDrag_Nom && dragCoefficient > referenceDrag_p1 &&  dragCoefficient < referenceDrag_p2)
        %Refine distribution region for parameter 1   
            standard_deviation_p1 = (((mean_value_p1+standard_deviation_p1) - (mean_value_p1-standard_deviation_p1)) - (2.0 * abs(mean_value_p1 - p1_bestCase)))/2.0;
            mean_value_p1 = p1_bestCase;
        
        %Update parameter 1
%             D = p1_bestCase;
%             Ubar = temp_p2;
        
    elseif (dragCoefficient < referenceDrag_Nom && dragCoefficient < referenceDrag_p1 &&  dragCoefficient > referenceDrag_p2)
        %Refine distribution region for parameter 2   
            standard_deviation_p2 = (((mean_value_p2+standard_deviation_p2) - (mean_value_p2-standard_deviation_p2)) - (2.0 * abs(mean_value_p2 - p2_bestCase)))/2.0;
            mean_value_p2 = p2_bestCase;
        
        %Update parameter 2
%             D = temp_p1;
%             Ubar = p2_bestCase;      
    end
    
    % Generate updated inputs
    input_p1 = distributionFunction_p1();
    input_p2 = distributionFunction_p2();
    
    %% Compute objective function and gradients
    J(i+1) = J_(referenceDrag_Nom,p1_bestCase,p2_bestCase);
    dJdp1(i+1) = dJdp1_(drag_dp1,p1_bestCase);  
    dJdp2(i+1) = dJdp2_(drag_dp2,p2_bestCase);  
    p1_hist(i+1) = p1_bestCase;
    p2_hist(i+1) = p2_bestCase;
    
    %% Update step size - Barzilai-Borwein step length for gradient descent
    denom = (dJdp1(i+1) - dJdp1(i))*2 + (dJdp2(i+1) - dJdp2(i))*2;
    gamma = ((p1_hist(i+1) - p1_hist(i)) * (dJdp1(i+1) - dJdp1(i)) + (p2_hist(i+1) - p2_hist(i)) * (dJdp2(i+1) - dJdp2(i)))/denom;

    D = temp_p1 - sign(dJdp1(i+1))*min(abs(gamma*dJdp1(i+1)), 0.1);
    Ubar = temp_p2 - sign(dJdp2(i+1))*min(abs(gamma*dJdp2(i+1)), 0.1);
    
    %% Increment iteration counter
        i = i+1;
        
    %% Update progress bar
        fprintf('\b|\n');
        
end

%End time count
disp(['Elapsed time: ', num2str(toc)])

%% Perform statistics to the output variables
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