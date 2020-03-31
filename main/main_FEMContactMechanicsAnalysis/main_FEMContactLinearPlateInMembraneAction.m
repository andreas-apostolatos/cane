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
addpath('../../FEMContactMechanicsAnalysis/graphics/',...
        '../../FEMContactMechanicsAnalysis/solvers',...
        '../../FEMContactMechanicsAnalysis/auxiliary/',...
        '../../FEMContactMechanicsAnalysis/contactSegments/',...
        '../../FEMContactMechanicsAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/';
caseName = 'bridge';
% caseName = 'cantilever_beam';
% caseName = 'wedge';
% caseName = 'hertz';

% Parse the data from the GiD input file
[strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC, analysis, parameters, ...
    propNLinearAnalysis, propStrDynamics, propGaussInt, propContact] = ...
    parse_StructuralModelFromGid(pathToCase, caseName, 'outputEnabled');

%% UI

% On the body forces
bodyForces = @computeConstantVerticalStructureBodyForceVct;

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Choose computation of the stiffness matrix
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST;
% computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed;

% Maximum number of contact iterations
propContact.maxIter = 50;

% On whether the case is a unit test
isUnitTest = false;

% On the geometry visualization
graph.visualization.geometry = 'current';

% Initialize graphics index
graph.index = 1;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Define boundary segments of the rigid contact wall
if strcmp(caseName,'bridge')
    % Amplitude of the externally applied boundary traction
    propNBC.tractionLoadVct = [0; - 1e4; 0];
    
    % Define the contact segments
    typeSegment = 'tessellation_circle'; % 'one', 'two', 'three', 'tessellation_circle'
    if strcmp(typeSegment, 'one')
        contactSegments.numSegments = 1;
        contactSegments.points = zeros(contactSegments.numSegments, 4);
        contactSegments.points(1,:) = [0.5 -0.5 3.5 -0.5];
    elseif strcmp(typeSegment, 'two')
        contactSegments.numSegments = 2;
        contactSegments.points = zeros(contactSegments.numSegments, 4);
        contactSegments.points(1,:) = [0.5 -0.5 2 -0.1];
        contactSegments.points(2,:) = [2 -0.1 3.5 -0.5];
    elseif strcmp(typeSegment, 'three')
        contactSegments.numSegments = 3;
        contactSegments.points = zeros(contactSegments.numSegments, 4);
        contactSegments.points(1,:) = [-0.5 -0.5 1.2 -0.1];
        contactSegments.points(2,:) = [1.2 -0.1 2 -0.1];
        contactSegments.points(3,:) = [2 -0.1 4.5 -0.5];
    elseif strcmp(typeSegment, 'tessellation_circle')
        center = [2,-4.1];
        radius = 4;
        startAngle = 3*pi/4;
        endAngle = pi/4;
        numSegments = 30;
        contactSegments = createCircleSegments ...
            (center, radius, startAngle, endAngle, numSegments);
    else
        error('Segment type %s not existent');
    end
elseif strcmp(caseName,'cantilever_beam')
    % Amplitude of the externally applied boundary traction
    propNBC.tractionLoadVct = [0; - 1e3; 0];
    
    % Geometric parameters for extending the rigid walls and render the
    % contact algorithm stable
    x_translation = 1e-0;  
    
    % Define the contact segments
    typeSegment = 'tessellation_circle'; % 'one', 'two', 'three', 'tessellation_circle'
    if strcmp(typeSegment, 'one')
        contactSegments.numSegments = 1;
        contactSegments.points = zeros(contactSegments.numSegments, 4);
        contactSegments.points(1,:) = [1.0+x_translation 0.4 10+x_translation 0.4];
    elseif strcmp(typeSegment, 'two')
        contactSegments.numSegments = 2;
        contactSegments.points = zeros(contactSegments.numSegments, 4);
        contactSegments.points(1,:) = [0.5 -0.5 2 -0.1];
        contactSegments.points(2,:) = [2 -0.1 3.5 -0.5];
    elseif strcmp(typeSegment, 'three')
        contactSegments.numSegments = 3;
        contactSegments.points = zeros(contactSegments.numSegments, 4);
        contactSegments.points(1,:) = [-0.5 -0.5 1.2 -0.1];
        contactSegments.points(2,:) = [1.2 -0.1 2 -0.1];
        contactSegments.points(3,:) = [2 -0.1 4.5 -0.5];
    elseif strcmp(typeSegment, 'tessellation_circle')
        center = [2,-4.1];
        radius = 4;
        startAngle = 3*pi/4;
        endAngle = pi/4;
        numSegments = 30;
        contactSegments = createCircleSegments ...
            (center, radius, startAngle, endAngle, numSegments);
    else
        error('Segment type %s not existent');
    end
elseif strcmp(caseName,'wedge')
    % Amplitude of the externally applied boundary traction
    propNBC.tractionLoadVct = [0; - 1e5; 0];

    % Geometric parameters for extending the rigid walls and render the
    % contact algorithm stable
    x_translation = 1e-3;
    y_translation = .1;
    x_extension = 1;
    
    % Define the contact segments
    contactSegments.numSegments = 3;
    contactSegments.points = zeros(contactSegments.numSegments, 4);
    contactSegments.points(1,:) = [-1 3.75 -1 -3];
    contactSegments.points(2,:) = [-0.33333+x_translation -2 1.2+x_translation 3.75];
    contactSegments.points(3,:) = [-1-x_extension -2+y_translation -0.33333+x_extension -2+y_translation];
elseif strcmp(caseName,'hertz')
    % Amplitude of the externally applied boundary traction
    propNBC.tractionLoadVct = [1e4; 0; 0];
    
    % Define the contact segments
    contactSegments.numSegments = 1;
    contactSegments.points = zeros(contactSegments.numSegments, 4);
    contactSegments.points(1,:) = [5 -1 5 5];
end

%% Compute normals to segments
contactSegments = computeUnitNormalVctsToSegments(contactSegments);

%% Compute the load vector
time = 0;
F = computeLoadVctFEMPlateInMembraneAction...
    (strMsh, analysis, propNBC, time, propGaussInt,'');

%% Plot the reference configuration
graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
    (strMsh, analysis, F, homDBC, contactSegments, graph, 'outputEnabled');

%% Solve the system for the displacement field and the Lagrange Multipliers fields
[dHat, lambdaHat, nodeIDs_active, numIter, FComplete, minElSize] = ...
    solveSignoriniFrictionlessContact2D...
    (analysis, strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC, bodyForces, ...
    parameters, contactSegments, computeStiffMtxLoadVct, solve_LinearSystem, ...
    propNLinearAnalysis, propContact , propGaussInt, caseName, pathToOutput,...
    isUnitTest, 'outputEnabled');

%% Postprocessing

% Plot the current configuration
resultant = 'stress';
component = '2Principal';
[graph.index, ~, ~] = plot_currentConfigurationAndResultants...
    (analysis, strMsh, homDBC, dHat, nodeIDs_active, contactSegments, ...
    parameters, resultant, component, graph);

% Compute the contact length, the contact force and the maximum contact
% pressure
[contactLength, contactForce, maxContactPressure] = ...
    computePostprocResultantsSignoriniFrictionlessContact2D...
    (strMsh, parameters, dHat, lambdaHat, nodeIDs_active);

% Get the length of the contact area and the reaction force on the contact
if strcmp(caseName,'hertz')    
    radius = 5;
    force = sum(F);
    hertzContactLength = sqrt(4*(2*force)*radius*((1 - parameters.nue^2)/parameters.E)/(pi*parameters.t));
    hertzPressure = 2*(2*force)/(parameters.t*pi*hertzContactLength);

    fprintf('\t The numerical contact length is: %.15d \n', contactLength);
    fprintf('\t The analytical (Hertz) contact length is: %.15d \n\n', hertzContactLength);

    fprintf('\t The maximum numerical pressure is: %.15d \n', maxContactPressure);
    fprintf('\t The maximum analytical (Hertz)  pressure is: %.15d \n', hertzPressure);
end

%% End of the script
