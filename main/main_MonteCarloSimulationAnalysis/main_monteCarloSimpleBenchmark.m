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

% Number of samples for the ultimate Monte Carlo Simulation
num_samples_ult = 1000;

% Number of Monte Carlo Simulations
num_mc = 1000;

% Number of samples
no_samples = zeros(num_mc, 1);

% Sample mean
mean_sample = zeros(num_mc, 1);

% Interval of confidence
alpha = 0.05;
interval_confidence_lower = zeros(num_mc, 1);
interval_confidence_upper = zeros(num_mc, 1);

% Mean value of the samples
mean_value_U = 0.0;

% Standard deviation of the samples
standard_deviation = 1.0;

% Plot the probability distribution function for the random input

% histogram(input,100,'Normalization','probability')
% histogram(input,100)
% figure(graph.index)
% histfit(U_var)
% yt = get(gca, 'YTick');
% set(gca, 'YTick', yt, 'YTickLabel', yt/numel(U_var))
% title('Input histogram for U')
% graph.index = graph.index + 1;

% Compute mean value and variance of the input data
% mean_value_U_est = mean(U_var);
% standard_deviation_U_est = std(U_var);
% fprintf('Input data possess mean value %d and standard deviation %d\n\n',mean_value_U_est, standard_deviation_U_est);

%% Define model
model_y = @(U) exp(U);

%% Loop over all the Monte Carlo Simulations
for iMC = 1:num_mc
    %% Assign the number of samples for the given MC Simulation
    no_samples(iMC,1) = iMC*num_samples_ult/num_mc;
    
    %% Initialize output
    Y_var = zeros(no_samples(iMC,1),1);
    
    %% Compute the samples for the load amplitude
    U_var = computeRandomDistribution(mean_value_U, standard_deviation, no_samples(iMC,1));
    
    %% Set up progress bar
    fprintf('Iteration : %d\n', iMC);
%     fprintf(['\n' repmat('.',1,no_samples(iMC,1)) '\n\n']);
%     tic

    %% Perform a Monte-Carlo Simulation for all samples
    parfor iSamples = 1:no_samples(iMC)  
        %% Compute the Y variable via the provided model
        Y_var(iSamples,1) = model_y(U_var(iSamples,1));

        %% Update progress bar
%         fprintf('\b|\n');
    end
%     disp(['Elapsed time: ', num2str(toc)]);
    
    %% Compute mean_samplethe mean value for the given MC Simulation
    mean_sample(iMC, 1) = mean(Y_var(:,1));
    standard_deviation_sample = std(Y_var(:,1));
    
    %% Compute the interval of confidence
    interval_confidence_lower(iMC, 1) = mean_sample(iMC, 1) - norminv(1 - alpha/2)*standard_deviation_sample/sqrt(no_samples(iMC,1));
    interval_confidence_upper(iMC, 1) = mean_sample(iMC, 1) + norminv(1 - alpha/2)*standard_deviation_sample/sqrt(no_samples(iMC,1));
    
    %% Clear variables
    clear U_var Y_var;
end

%% Perform statistics to the output variables
% figure(graph.index)
% histfit(Y_var)
% yt = get(gca, 'YTick');
% set(gca, 'YTick', yt, 'YTickLabel', yt/numel(Y_var))
% title('Output histogram for Y')
% graph.index = graph.index + 1;

% Compute mean value and variance of the input data
% mean_value_Y = mean(Y_var);
% standard_deviation_Y = std(Y_var);
% fprintf('\n\nOutput data possess mean value %d (%d) and standard deviation %d (%d)\n\n',mean_value_Y, sqrt(exp(1)), standard_deviation_Y, exp(1)*(exp(1)-1));

%% Plot the sample mean
figure(graph.index)
plot(no_samples, ones(length(no_samples), 1)*sqrt(exp(1)), '-black', no_samples, mean_sample, '-blue',...
    no_samples, interval_confidence_lower, '-.black', no_samples, interval_confidence_upper, '-.black');
legend('Exact solution', 'MC solution');
title('Sample mean');
graph.index = graph.index + 1;

%% END OF THE SCRIPT