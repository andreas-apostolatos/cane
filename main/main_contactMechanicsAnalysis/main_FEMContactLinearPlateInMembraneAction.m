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

% Add all functions related to signorini frictionless contact problem
addpath('../../contactMechanicsAnalysis/plot',...
        '../../contactMechanicsAnalysis/solvers',...
        '../../contactMechanicsAnalysis/supportFunctions');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/';
% caseName = 'example_01_bridge';
caseName = 'example_02_wedge';

% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,propStrDynamics,gaussInt,contactNodes] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

%% GUI - on the graph
graph.index = 1;

% On the geometry visualization
graph.visualization.geometry = 'current';

%% Rigid wall- line  [(x0,y0) ; (x1,y1)]

% different line segments for different cases
if strcmp(caseName,'example_01_bridge')
    % define bottom contact line segment
    wall_1 = [-0.5, -0.5; 2,-0.1];
    wall_2 = [2, -0.1; 4.5,-0.5];

    % add a wall to the segments of points
    segments.points(:,:,1) = wall_1;
    segments.points(:,:,2) = wall_2;
    
%     segments = createCircleSegments(2,-5.1,5,19);
    
elseif strcmp(caseName,'example_02_wedge')
    % define bottom contact line segment
    wall_1 = [-1, 3.75; -1, -2];
    wall_2 = [-0.33333, -2; 1.2, 3.75];
    
    % add a wall to the segments of points
    segments.points(:,:,1) = wall_1;
    segments.points(:,:,2) = wall_2;

end
% computeConstantVerticalLoad
%% Compute the load vector
time = 0;
F = computeLoadVctFEMPlateInMembraneAction(strMsh,NBC,time,gaussInt,'outputEnabled');

%% Output data to a VTK format
% pathToOutput = '../../outputVTK/FEMContactLinearPlateInMembraneAction/';
% not implemented yet

%% Visualization of the configuration
graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
    (strMsh,F,homDBC,graph,'outputEnabled');

% plot the wall segment
plot_segments(segments); 

%% Solve the system and get the displacement field
ts = cputime;

maxIteration = 30;

[displacement,lagrange] = solveSignoriniLagrange_1(strMsh,homDBC,contactNodes,F,segments,parameters,analysis,maxIteration); 
%[displacement,lagrange] = solveSignoriniLagrange_2(strMsh,homDBC,contactNodes,F,segments,parameters,analysis,maxIteration);

fprintf('\t Time: %4.2f \n',cputime-ts);

%% Postprocessing
graph.index = plot_currentConfigurationFEMPlateInMembraneAction(strMsh,homDBC,displacement,graph);
plot_segments(segments);
plot_activeNodes(strMsh,displacement,lagrange); 

%% End of the script