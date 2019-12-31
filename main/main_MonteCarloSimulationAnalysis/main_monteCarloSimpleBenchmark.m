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

% Add all functions related to the Monte Carlo Simulation
addpath('../../MonteCarloSimulationAnalysis/probabilityDistributionFunctions/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');

%% GUI

% Initialize graphics index
graph.index = 1;

% Sampling function selection
% 'generateRandomUniformDistribution', 'generateRandomNormalDistribution', 
% 'generateRandomLatinHypercubeDistribution',
% 'generateRandomQuasiMonteCarloDistribution'
computeRandomDistribution = @generateRandomNormalDistribution;

%% Generate input samples

% Number of samples
no_samples = 1000;

% Mean value of the samples
mean_value_U = 0.0;

% Standard deviation of the samples
standard_deviation = 1.0;

% Uncertainty in the load amplitude
U_var = computeRandomDistribution(mean_value_U, standard_deviation, no_samples);

% Plot the probability distribution function for the random input

% histogram(input,100,'Normalization','probability')
% histogram(input,100)
figure(graph.index)
histfit(U_var)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(U_var))
title('Input histogram for U')
graph.index = graph.index + 1;

% Compute mean value and variance of the input data
mean_value_U = mean(U_var);
standard_deviation_U = std(U_var);
fprintf('Input data possess mean value %d and standard deviation %d\n\n',mean_value_U, standard_deviation_U);

%% Define model
model_y = @(U) exp(U);

%% Initialize output statistics

% Support forces in the x-direction
Y_var = zeros(no_samples,1);

%% Set up progress bar
fprintf(['\n' repmat('.',1,no_samples) '\n\n']);
tic

%% Perform a Monte-Carlo Simulation for all samples
parfor i = 1:no_samples    
    %% Compute the Y variable via the provided model
    Y_var(i,1) = model_y(U_var(i,1));
    
    %% Update progress bar
    fprintf('\b|\n');
end
disp(['Elapsed time: ', num2str(toc)]);

%% Perform statistics to the output variables
figure(graph.index)
histfit(Y_var)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt/numel(Y_var))
title('Output histogram for Y')
graph.index = graph.index + 1;

% Compute mean value and variance of the input data
mean_value_Y = mean(Y_var);
standard_deviation_Y = std(Y_var);
fprintf('\n\nOutput data possess mean value %d (%d) and standard deviation %d (%d)\n\n',mean_value_Y, sqrt(exp(1)), standard_deviation_Y, exp(1)*(exp(1)-1));

%% END OF THE SCRIPT