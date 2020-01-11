%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%% Script documentation
%
% Task : xxxxxx
%        xxxxxx
%
% Date : 10.01.2020
%
%% Preamble
clear;
clc;
close all;

%% Includes
addpath('../nurbs','../nurbs/nurbs_base_vectors','../nurbs/nurbs_basis_functions',...
        '../nurbs/nurbs_refinement','../nurbs/nurbs_geometry','../mesh',...
        '../mesh2D','../nurbs/nurbs_graphics','../nurbs/nurbs_surface',...
        '../algebra','../nurbs/nurbs_curve','../basisFunctions','../plot','../load',...
        '../supports','../stiffnessMatrices','../plane_stress_analysis',...
        '../boundaryConditions','../planeStressAnalysis','../lagrangeMultipliers');
    
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

% Include performance optimzed functions
addpath('../../efficientComputation/');




% must - checked by marko   
addpath('../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors')
addpath('../../FEMPlateInMembraneActionAnalysis/loads');
addpath('../../basisFunctions');


%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/';
caseName = 'example_01_noGap';
% caseName = 'example_02_withGap';

% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,propStrDynamics] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

%% Material constants

% Young's modulus
parameters.E = 1e5;

parameters.rho = 1.0;

% Poisson ratio
parameters.nue = 0.3;

% Thickness of the plate
parameters.t = 1;

%% GUI

% Analysis type
analysis.dimension = '2d';
analysis.dofs = 'displacements';
analysis.physics = 'plainStrain';
analysis.type = 'linear';

% On the graph
graph.index = 1;

figure(graph.index)

% On the geometry visualization
graph.visualization.geometry = 'current';

% Quadrature for the stiffness matrix and the load vector of the problem
% 'default', 'user'
intLoad.type = 'default';
intDomain.type = 'default';
intLoad.noGP = 1;
intDomain.noGP = 1;

%% rigid wall- line  [(x0,y0) ; (x1,y1)]

wall=[-5, 0; 5 ,0];

segmentsPoints(1,:,:)=wall;


%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMContactLinearPlateInMembraneAction/';

%% Compute the load vector
t = 0;
F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,intLoad,'outputEnabled');

%% Visualization of the configuration
graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
    (strMsh,analysis,F,homDBC,graph,'outputEnabled');


% On the Contact boundary (either the whole boundary is a contact area or
% with the function getContactNodes selected regions :

% get this from GiD

%candidateContactNodes=getContactNodesOnBoundary(p,U,q,V,CP,boundaryConditions,mesh);% Select certain area
candidateContactNodes=strMsh.boundaryNodes;% Define the whole boundary as possible contact area









%% Plot the initial configuration
hold on;
index = plotBoundaryConditionsOnMesh(strMsh,homDBC,F,graph);
plotSegments(strMsh,wall); % plot the wall segment
hold off;

%% Solve the system and get the displacement field
ts=cputime;

maxIterForLagrange=100;

[displacement, lagrange] = solveSignoriniLagrange1(strMsh,homDBC,candidateContactNodes,F,segmentsPoints,parameters,analysis,maxIterForLagrange); % Use Algorithm 1 
%[displacement,lagrange] = solveSignoriniLagrange2(mesh,homDBC,candidateContactNodes,F,segmentsPoints,materialProperties,analysis,maxIterForLagrange); % or use Algorithm 2

fprintf('\t Time :   %4.2f \n',cputime-ts);
%% Postprocessing
graph.index = plotCurrentConfigurationAndResultants(strMsh,homDBC,displacement,graph);
plotSegments(strMsh,segmentsPoints);
plotLagrangeMultipliers( strMsh, displacement,lagrange.active_nodes,lagrange.multipliers,'' );% To show vertical bars for the Lagrange multiplier values insert 'v' as last parameter

%% End of the script