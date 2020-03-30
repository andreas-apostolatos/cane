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
        '../../FEMContactMechanicsAnalysis/postprocessing/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/';
% caseName = 'example_01_bridge';
caseName = 'example_02_wedge';
% caseName = 'example_03_hertz';
% caseName = 'example_04_cantilever_beam_contact';
% caseName = 'example_05_rectangle_contact_debug';
% caseName = 'example_06_wedge_contact_debug';

% Parse the data from the GiD input file
[strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC, analysis, parameters, ...
    propNLinearAnalysis, propStrDynamics, propGaussInt, propContact] = ...
    parse_StructuralModelFromGid(pathToCase, caseName, 'outputEnabled');
%     computeConstantVerticalLoad
%     computeConstantHorizontalLoad

%% UI

% On the body forces
bodyForces = @computeConstantVerticalStructureBodyForceVct;

% Amplitude of the externally applied boundary traction
propNBC.tractionLoadVct = [0; - 1e5; 0];

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Choose computation of the stiffness matrix
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST;
% computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed;

% Update the contact properties
propContact.maxIter = 30;

% On whether the case is a unit test
isUnitTest = false;

% On the geometry visualization
graph.visualization.geometry = 'current';

% Initialize graphics index
graph.index = 1;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Rigid wall- line  [(x0,y0) ; (x1,y1)]

% different line segments for different cases
if strcmp(caseName,'example_01_bridge') || strcmp(caseName,'example_04_cantilever_beam_contact') || ...
        strcmp(caseName,'example_05_rectangle_contact_debug')
   
    % Either define bottom contact line segment and add it to the segments
%     contactSegments.points(:,:,1) = [-0.5, -0.5; 1.2,-0.1];
%     contactSegments.points(:,:,2) = [1.2, -0.1; 2,-0.1];
%     contactSegments.points(:,:,3) = [2, -0.1; 4.5,-0.5];
    
%     contactSegments.points(:,:,1) = [0.5, -0.5; 2,-0.1];
%     contactSegments.points(:,:,2) = [2, -0.1; 3.5,-0.5];
    
%     contactSegments.points(:,:,1) = [0.5, -0.5; 3.5,-0.5];

% For the cantilever case
%     contactSegments.points(:,:,1) = [1.0, 0.4; 10, 0.4];

% For the debug case
%     contactSegments.points(:,:,1) = [-10.0, 0; 10, 0];
    
    
    % ...or define a circular contact boundary
    center = [2,-4.1];
    radius = 4;
    startAngle = 3*pi/4;
    endAngle = pi/4;
    numSegments = 30; %17
    
    % Create circular segments
    contactSegments = createCircleSegments(center,radius,startAngle,endAngle,numSegments);
elseif strcmp(caseName,'example_02_wedge') || strcmp(caseName,'example_06_wedge_contact_debug')
    x_translation = 1e-3;
     
    y_translation = .1;
     
    x_extension = 1;
    
    % define bottom contact line segment and add it to the segments
    contactSegments.numSegments = 3;
    contactSegments.points = zeros(contactSegments.numSegments, 4);
    contactSegments.points(1,:) = [-1 3.75 -1 -3];
    contactSegments.points(2,:) = [-0.33333+x_translation -2 1.2+x_translation 3.75];
    contactSegments.points(3,:) = [-1-x_extension -2+y_translation -0.33333+x_extension -2+y_translation];

elseif strcmp(caseName,'example_03_hertz')
    % define bottom contact line segment
    wall_1 = [5, -1; 5, 5];
    
    % add a wall to the segments of points
    contactSegments.points(:,:,1) = wall_1;    

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

%% Solve the system and get the displacement field
[dHat, lambdaHat, nodeIDs_active] = solveSignoriniFrictionlessContactProblem2D...
    (analysis, strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC,bodyForces, ...
    parameters, contactSegments, computeStiffMtxLoadVct, solve_LinearSystem, ...
    propNLinearAnalysis, propContact , propGaussInt, caseName, pathToOutput,...
    isUnitTest, 'outputEnabled');

%% Postprocessing
graph.index = plot_currentConfigurationFEMPlateInMembraneAction...
    (strMsh, homDBC, contactSegments, dHat, graph);
hold on;
plot_activeNodes(strMsh, dHat, nodeIDs_active);
hold off;

% Get the length of the contact area and the reaction force on the contact
if strcmp(caseName,'example_03_hertz')
    [contactLength,contactForce,maxContactPressure] = ...
        computeContactResultants...
        (strMsh, parameters, dHat, lambdaHat, nodeIDs_active);
    
    radius = 5;
    force = sum(F);
    hertzContactLength = sqrt(4*(2*force)*radius*((1 - parameters.nue^2)/parameters.E)/(pi*parameters.t));
    hertzPressure = 2*(2*force)/(parameters.t*pi*hertzContactLength);

    fprintf('\t The COMPUTED contact length is: %f \n',contactLength);
    fprintf('\t The (HERTZ)  contact length is: %f \n\n',hertzContactLength);

    fprintf('\t The maximal COMPUTED pressure is: %f \n',maxContactPressure);
    fprintf('\t The maximal (HERTZ)  pressure is: %f \n',hertzPressure);
end

%% End of the script