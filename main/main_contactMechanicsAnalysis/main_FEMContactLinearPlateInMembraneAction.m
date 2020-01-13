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
% Task : Analysis of Signorini contact problem
%
% Date : 12.01.2020
%
%% Preamble
clear;
clc;
close all;

%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the low order basis functions
addpath('../../basisFunctions/');

% Add all equation system solvers
addpath('../../equationSystemSolvers/');

% Add all the efficient computation functions
addpath('../../efficientComputation/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/',...
        '../../FEMPlateInMembraneActionAnalysis/errorComputation/');

% NOTE - check what is not neaded and rearange it
% Add all functions related to signorini frictionless contact problem
addpath('../../contactMechanicsAnalysis/plot',...
        '../../contactMechanicsAnalysis/planeStressAnalysis',...
        '../../contactMechanicsAnalysis/lagrangeMultipliers');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/';
caseName = 'example_01_noGap';
% caseName = ... NOTE - add all test cases

% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,propStrDynamics,gaussInt,contactNodes] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

%% GUI
% Analysis type -> get from parser, here for convenience
analysis.type = 'planeStress';

% On the graph
graph.index = 1;

% On the geometry visualization
graph.visualization.geometry = 'current';

%% rigid wall- line  [(x0,y0) ; (x1,y1)]
% define bottom contact line segment
%wall_1 = [-5, 0; 5,0];
wall_2 = [-5, -0.1; 5,-0.1];

% add a wall to the segments of points
segmentPoints(:,:,1) = wall_2;
%segmentPoints(:,:,2) = wall_1;

% Define the structrure array
% candidateNodes(1)=struct('indices',cn1);
% candidateNodes(2)=struct('indices',cn2);

%% Compute the load vector
time = 0;
F = computeLoadVctFEMPlateInMembraneAction(strMsh,NBC,time,gaussInt,'outputEnabled');

%% Output data to a VTK format
% pathToOutput = '../../outputVTK/FEMContactLinearPlateInMembraneAction/';
% not implemented yet

%% Visualization of the configuration
% NOTE - combine this function 
graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
    (strMsh,F,homDBC,graph,'outputEnabled');

% plot the wall segment
plot_segments(segmentPoints); 

%% Solve the system and get the displacement field
ts = cputime;

maxIteration = 100;

[displacement, lagrange] = solveSignoriniLagrange1(strMsh,homDBC,contactNodes,F,segmentPoints,parameters,analysis,maxIteration); % Use Algorithm 1 
%[displacement,lagrange] = solveSignoriniLagrange2(strMsh,homDBC,contactNodes,F,segmentsPoints,parameters,analysis,maxIteration); % or use Algorithm 2

% fix this so it will also work with only one line
%[displacement, lagrange] = multSolveSignoriniLagrange1(strMsh,homDBC,contactNodes,F,segmentPoints,parameters,analysis,maxIteration); % Use Algorithm 1 

fprintf('\t Time :   %4.2f \n',cputime-ts);
%% Postprocessing
graph.index = plot_currentConfigurationFEMPlateInMembraneAction(strMsh,homDBC,displacement,graph);
plot_segments(segmentPoints);
% To show vertical bars for the Lagrange multiplier values insert 'v' as last parameter
plot_lagrangeMultipliers(strMsh,displacement,lagrange,''); 

%% End of the script