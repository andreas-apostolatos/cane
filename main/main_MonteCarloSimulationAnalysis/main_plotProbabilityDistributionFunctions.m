%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Script documentation
%
% Task : Demonstrate different probability distribution functions
%
% Date : 02.01.2019
%
%% Preamble
clear;
clc;
close all;

%% Includes

% Add all functions related to the Monte Carlo Simulation
addpath('../../MonteCarloSimulationAnalysis/probabilityDistributionFunctions/');

%% Input data
meanValue = 10;
standardDeviation = 1;
noSamples = 10000;
% number of bins of histograms
noBins = 50;

%% Different sampling functions
randomUniformDistribution = generateRandomUniformDistribution...
                                 (meanValue, standardDeviation, noSamples);
randomNormalDistribution = generateRandomNormalDistribution...
                                 (meanValue, standardDeviation, noSamples);
quasiMonteCarloDistribution = generateRandomQuasiMonteCarloDistribution...
                                 (meanValue, standardDeviation, noSamples);
latinHypercubeDistribution = generateRandomLatinHypercubeDistribution...
                                 (meanValue, standardDeviation, noSamples);

%% Visualize the sampling methods via histograms
figure(1);
histfit(randomUniformDistribution, noBins);
title('Random Uniform Distribution');

figure(2);
histfit(randomNormalDistribution, noBins);
title('Random Normal Distribution');

figure(3);
histfit(quasiMonteCarloDistribution, noBins);
title('Quasi Monte Carlo Distribution');

figure(4);
histfit(latinHypercubeDistribution, noBins);
title('Latin Hypercube Distribution');