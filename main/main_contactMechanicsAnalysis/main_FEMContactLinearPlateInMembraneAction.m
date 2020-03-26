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
% caseName = 'example_02_wedge';
%caseName = 'example_03_hertz';

% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,propStrDynamics,gaussInt,propContact] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

%% UI

% On the body forces
bodyForces = @computeConstantVerticalStructureBodyForceVct;

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
if strcmp(caseName,'example_01_bridge')
   
    % Either define bottom contact line segment and add it to the segments
%     contactSegments.points(:,:,1) = [-0.5, -0.5; 1.2,-0.1];
%     contactSegments.points(:,:,2) = [1.2, -0.1; 2,-0.1];
%     contactSegments.points(:,:,3) = [2, -0.1; 4.5,-0.5];
    
%     contactSegments.points(:,:,1) = [0.5, -0.5; 2,-0.1];
%     contactSegments.points(:,:,2) = [2, -0.1; 3.5,-0.5];
    
    contactSegments.points(:,:,1) = [0.5, -0.5; 3.5,-0.5];
    
    
    % ...or define a circular contact boundary
    center = [2,-4.1];
    radius = 4;
    startAngle = 3*pi/4;
    endAngle = pi/4;
    nSegments = 17; %17
    
    % Create circular segments
%     contactSegments = createCircleSegments(center,radius,startAngle,endAngle,nSegments);
%     computeConstantVerticalLoad
%     computeConstantHorizontalLoad
elseif strcmp(caseName,'example_02_wedge')
    
    % define bottom contact line segment and add it to the segments
    contactSegments.points(:,:,1) = [-1, 3.75; -1, -2];
    contactSegments.points(:,:,2) = [-0.33333, -2; 1.2, 3.75];
    contactSegments.points(:,:,3) = [-1, -2; -0.33333, -2];

elseif strcmp(caseName,'example_03_hertz')
    % define bottom contact line segment
    wall_1 = [5, -1; 5, 5];
    
    % add a wall to the segments of points
    contactSegments.points(:,:,1) = wall_1;    

end

%% Compute normals to segments
contactSegments = buildSegmentsData(contactSegments);

%% Compute the load vector
time = 0;
F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,time,gaussInt,'outputEnabled');

%% Plot the reference configuration
graph.index = plot_referenceConfigurationFEMPlateInMembraneAction...
    (strMsh,analysis,F,homDBC,contactSegments,graph,'outputEnabled');

%% Solve the system and get the displacement field
[dHat,lambdaHat,nodeIDs_active] = solveSignoriniLagrange_1_debug...
    (analysis,strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,bodyForces,parameters,...
    contactSegments,computeStiffMtxLoadVct,solve_LinearSystem,...
    propNLinearAnalysis,propContact,gaussInt,caseName,pathToOutput,...
    isUnitTest,'outputEnabled',graph);
% [dHat,lambdaHat,nodeIDs_active] = solveSignoriniLagrange_1...
%     (analysis,strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,bodyForces,parameters,...
%     contactSegments,computeStiffMtxLoadVct,solve_LinearSystem,...
%     propNLinearAnalysis,propContact,gaussInt,caseName,pathToOutput,...
%     isUnitTest,'outputEnabled');
% [dHat,lambdaHat,nodeIDs_active] = solveSignoriniLagrange_2...
%     (analysis,strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,bodyForces,parameters,...
%     contactSegments,computeStiffMtxLoadVct,solve_LinearSystem,...
%     propNLinearAnalysis,propContact,gaussInt,caseName,pathToOutput,...
%     isUnitTest,'outputEnabled'); 
%% Postprocessing
% graph.index = plot_currentConfigurationFEMPlateInMembraneAction(strMsh,homDBC,contactSegments,dHat,graph);
% plot_activeNodes(strMsh,dHat,nodeIDs_active);

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