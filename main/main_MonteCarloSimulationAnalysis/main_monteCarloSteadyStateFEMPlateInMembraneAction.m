%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Performs Multilevel Monte-Carlo for the steady-state linear
%        elastic plane stress problem
%
% Date : 24.12.2019
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

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');

%% Define the path to the case
pathToCase = '../../inputGiD/FEMPlateInMembraneActionAnalysis/';
caseName = 'cantileverBeamPlaneStress';

%% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,~] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

%% GUI

% On the body forces
bodyForces = @computeConstantVerticalBodyForceVct;

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Choose computation of the stiffness matrix
% computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST;
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed;

% Quadrature for the stiffness matrix and the load vector of the problem
% 'default', 'user'
intLoad.type = 'default';
intDomain.type = 'default';
intLoad.noGP = 1;
intDomain.noGP = 1;

% Quadrature for the L2-norm of the error
intError.type = 'user';
intError.noGP = 4;

% Linear analysis
propStrDynamics = 'undefined';

% On whether the case is a unit test
isUnitTest = true;

% Initialize graphics index
graph.index = 1;

%% Define the name of the vtk file from where to resume the simulation
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Input statistics

% Number of samples
no_samples = 1000;

% Mean value of the samples
mean = -1e3;

% Standard deviation of the samples
standard_deviation = 1e1;

% Sampling function selection
sampling_method_list = {'randomUniform','randomNormal','latinHyperCube', 'quasiMonteCarloHalton'};
sampling_method = sampling_method_list{2};
switch sampling_method
    case 'randomUniform'
        % Uniform distribution with bounds
        distributionFunction = @(varargin) 2*standard_deviation*rand(no_samples,1) + mean - standard_deviation;
    case 'randomNormal'
        % Random normal distribution sampling
        distributionFunction = @(varargin) normrnd(mean, standard_deviation, no_samples, 1);
    case 'latinHyperCube'
        % Latin hypercube sampling
        distributionFunction =  @(varargin) lhsnorm(mean, (standard_deviation*standard_deviation), no_samples);
    case 'quasiMonteCarloHalton'
        % Quasi Monte Carlo Halton Sequence sampling
        p = haltonset(1,'Skip',1e3,'Leap',1e2);
        p = scramble(p,'RR2');
        haltonvector = net(p,no_samples); % Halton Sequence sampling
        distributionFunction = @(varargin) norminv(haltonvector, mean, standard_deviation); % Invert uniform Halton sequence to normal distribution
    otherwise
        % Error
        error('Invalid input sampling method')
end

%% Generate input samples

% Uncertainty in the load amplitude
amplitude = distributionFunction();

% Re-initialize the function handle array in the NBC struct
NBC = rmfield(NBC,'fctHandle');
no_NBC = length(NBC.loadType(:,1));
NBC.fctHandle = cell(no_NBC,1);
NBC_random = repmat(NBC,no_samples);

% Plot the probability distribution function for the random input

% histogram(input,100,'Normalization','probability')
% histogram(input,100)
figure(graph.index)
histfit(amplitude)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(amplitude))
title('Input histogram')
graph.index = graph.index + 1;

%% Initialize output statistics

% Support forces in the x-direction
support_forces_x = zeros(no_samples,1);

% Support forces in the y-direction
support_forces_y = zeros(no_samples,1);

%% Set up progress bar
fprintf(['\n' repmat('.',1,no_samples) '\n\n']);
tic

%% Perform a Monte-Carlo Simulation for all samples
parfor i = 1:no_samples
    %% Compute the force vector subject to a load with random amplitude
    for j = 1:no_NBC
        NBC_random(i).fctHandle{j} = @(x,y,z,t) [0; amplitude(i,1); 0];
    end
    F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC_random(i),0,intLoad,'');
    
    %% Solve the CFD problem
    [~,FComplete,minElSize] = solve_FEMPlateInMembraneAction...
        (analysis,strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC_random(i),bodyForces,...
        parameters,computeStiffMtxLoadVct,solve_LinearSystem,...
        propNLinearAnalysis,propStrDynamics,intDomain,caseName,pathToOutput,...
        isUnitTest,'');
    
    %% Compute the support reaction forces
    support_forces_x(i,1) = sum(FComplete(homDBC(1:2:end),1));
    support_forces_y(i,1) = sum(FComplete(homDBC(2:2:end),1));
    
    %% Update progress bar
    fprintf('\b|\n');
end
disp(['Elapsed time: ', num2str(toc)])

%% Perform statistics to the output variables

% Drag coefficient
figure(graph.index)
histfit(support_forces_x)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(support_forces_x))
title('Output histogram for the support forces in x-direction')
graph.index = graph.index + 1;

% Lift coefficient
figure(graph.index)
histfit(support_forces_y)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(support_forces_y))
title('Output histogram for the support forces in y-direction')
graph.index = graph.index + 1;

%% END OF THE SCRIPT