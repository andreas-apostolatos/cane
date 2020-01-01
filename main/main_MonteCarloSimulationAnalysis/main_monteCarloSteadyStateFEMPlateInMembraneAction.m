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
    
% Add all functions related to the Monte Carlo Simulation
addpath('../../MonteCarloSimulationAnalysis/probabilityDistributionFunctions/',...
        '../../MonteCarloSimulationAnalysis/solvers/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');

%% Define the path to the case
pathToCase = '../../inputGiD/FEMPlateInMembraneActionAnalysis/';
caseName = 'cantileverBeamPlaneStress';
% caseName = 'PlateWithMultipleHolesPlaneStress';

%% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,~] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

% Enforce geometrically linear analysis
propNLinearAnalysis.method = 'UNDEFINED';
propNLinearAnalysis.noLoadSteps = NaN;
propNLinearAnalysis.eps = NaN;
propNLinearAnalysis.maxIter = NaN;

%% GUI

% On the body forces
bodyForces = @computeConstantVerticalBodyForceVct;

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Choose computation of the stiffness matrix
% computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST;
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed;

% Sampling function selection
% 'generateRandomUniformDistribution', 'generateRandomNormalDistribution', 
% 'generateRandomLatinHypercubeDistribution',
% 'generateRandomQuasiMonteCarloDistribution'
computeRandomDistribution = @generateRandomNormalDistribution;

% Quadrature for the stiffness matrix and the load vector of the problem
% 'default', 'user'
propIntLoad.type = 'default';
propIntDomain.type = 'default';
propIntLoad.noGP = 1;
propIntDomain.noGP = 1;

% Quadrature for the L2-norm of the error
intError.type = 'user';
intError.noGP = 4;

% Linear analysis
propStrDynamics = 'undefined';

% On whether the case is a unit test
isUnitTest = false;

% Initialize graphics index
graph.index = 1;

%% Define the name of the vtk file from where to resume the simulation
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Generate input samples

% Number of samples
no_samples = 1000;

% Mean value of the samples
mean_value = -1e3;
% mean_value = -1e2;

% Standard deviation of the samples
% standard_deviation = 1e1;
standard_deviation = 1e0;

% Uncertainty in the load amplitude
amplitude = computeRandomDistribution(mean_value, standard_deviation, no_samples);

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

% Compute mean value and variance of the input data
mean_value_amplitude = mean(amplitude);
standard_deviation_amplitude = std(amplitude);
fprintf('Input data possess mean value %d and standard deviation %d\n\n',mean_value_amplitude, standard_deviation_amplitude);

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
%         NBC_random(i).fctHandle{j} = @(x,y,z,t) [0; amplitude(i,1); 0];
        NBC_random(i).fctHandle{j} = @(x,y,z,t) [amplitude(i,1); 0; 0];
    end
    F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC_random(i),0,propIntLoad,'');

    %% Solve the CFD problem
    [~,FComplete,minElSize] = solve_FEMPlateInMembraneAction...
        (analysis,strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC_random(i),bodyForces,...
        parameters,computeStiffMtxLoadVct,solve_LinearSystem,...
        propNLinearAnalysis,propStrDynamics,propIntDomain,caseName,pathToOutput,...
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
% figure(graph.index)
% histfit(support_forces_x)
% yt = get(gca, 'YTick');
% set(gca, 'YTick', yt, 'YTickLabel', yt/numel(support_forces_x))
% title('Output histogram for the support forces in x-direction')
% graph.index = graph.index + 1;

% Lift coefficient
figure(graph.index)
histfit(support_forces_y)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(support_forces_y))
title('Output histogram for the support forces in y-direction')
graph.index = graph.index + 1;

% Compute mean value and variance of the input data
mean_value_support_forces_y = mean(support_forces_y);
standard_deviation_support_forces_y = std(support_forces_y);
fprintf('\n\nOutput data possess mean value %d and standard deviation %d\n\n',mean_value_support_forces_y, standard_deviation_support_forces_y);

%% END OF THE SCRIPT