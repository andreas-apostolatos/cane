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
caseName = 'example_01_bridge';
%caseName = 'example_02_wedge';
%caseName = 'example_03_hertz';


% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,propStrDynamics,gaussInt,propContact] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

%% GUI - on the graph
graph.index = 1;

% On the geometry visualization
graph.visualization.geometry = 'current';

%% Rigid wall- line  [(x0,y0) ; (x1,y1)]

% different line segments for different cases
if strcmp(caseName,'example_01_bridge')
   
    % Either define bottom contact line segment and add it to the segments
%     segments.points(:,:,1) = [-0.5, -0.5; 1.2,-0.1];
%     segments.points(:,:,2) = [1.2, -0.1; 2,-0.1];
%     segments.points(:,:,3) = [2, -0.1; 4.5,-0.5];
    
    % ...or define a circular contact boundary
    center = [2,-4.1];
    radius = 4;
    startAngle = 3*pi/4;
    endAngle = pi/4;
    nSegments = 20;
    
    % Create circular segments
    segments = createCircleSegments(center,radius,startAngle,endAngle,nSegments);
   
elseif strcmp(caseName,'example_02_wedge')
    
    % define bottom contact line segment and add it to the segments
    segments.points(:,:,1) = [-1, 3.75; -1, -2];
    segments.points(:,:,2) = [-0.33333, -2; 1.2, 3.75];
    segments.points(:,:,3) = [-1, -2; -0.33333, -2];

elseif strcmp(caseName,'example_03_hertz')
    % define bottom contact line segment
    wall_1 = [5, -1; 5, 5];
    
    % add a wall to the segments of points
    segments.points(:,:,1) = wall_1;    
    
end

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

[displacement,lagrange] = solveSignoriniLagrange_1(strMsh,homDBC,propContact,F,segments,parameters,analysis,maxIteration,'outputEnabled'); 
%[displacement,lagrange] = solveSignoriniLagrange_2(strMsh,homDBC,propContact,F,segments,parameters,analysis,maxIteration,'outputEnabled');

fprintf('\t Time: %4.2f \n',cputime-ts);

%% Postprocessing
graph.index = plot_currentConfigurationFEMPlateInMembraneAction(strMsh,homDBC,displacement,graph);
plot_segments(segments);
plot_activeNodes(strMsh,displacement,lagrange); 

% Get the length of the contact area and the reaction force on the contact
if strcmp(caseName,'example_03_hertz')
    [contactLength,contactForce,maxContactPressure] = ...
        computeContactResultants(strMsh,displacement,lagrange,parameters);
    
    radius = 5;
    force = sum(F);
    hertzContactLength = sqrt(4*(2*force)*radius*((1-parameters.nue^2)/parameters.E)/(pi*parameters.t));
    hertzPressure = 2*(2*force)/(parameters.t*pi*hertzContactLength);

    fprintf('\t The COMPUTED contact length is: %f \n',contactLength);
    fprintf('\t The (HERTZ)  contact length is: %f \n\n',hertzContactLength);

    fprintf('\t The maximal COMPUTED pressure is: %f \n',maxContactPressure);
    fprintf('\t The maximal (HERTZ)  pressure is: %f \n',hertzPressure);
end

%% End of the script