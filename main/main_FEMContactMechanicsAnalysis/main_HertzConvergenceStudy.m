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
% Date : 01.04.2020
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
maxIter = 50;

% On whether the case is a unit test
isUnitTest = false;

% Enable output
outMsg = 'outputEnabled';
% outMsg = '';

% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/refinementStudyHertz/';

%% Initialize variables for convergence study

% Number of cases defined in GiD input folder
numCases = 4;
zeroVector = zeros(numCases,1);

% Array of Hertz (reference) contact lengths and pressures
hertzContactLength = zeroVector;
hertzPressure = zeroVector;
appliedForce = zeroVector;

% Intialize array of calculated post-processing resultants
contactLength = zeroVector;
contactForce = zeroVector;
maxContactPressure = zeroVector;

% Intialize array of number of elements
numberOfElements = zeroVector;
clear zeroVector;

%% Define boundary segments of the rigid contact wall
contactSegments.numSegments = 1;
contactSegments.points = zeros(contactSegments.numSegments, 4);
contactSegments.points(1,:) = [5 -1 5 5];

% Define the radius of Hertz benchmark problem setup
radius = 5;

%% Compute normals to segments
contactSegments = computeUnitNormalVctsToSegments(contactSegments);

%% Amplitude of the externally applied boundary traction
tractionLoadVct = [5e1
                   0
                   0];

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/refinementStudyHertz/';

% Collect all case names into an array of case names
% refined mesh with reduced element size at the tip
caseName = cellstr(num2str((numCases:-1:1)', 'refinementStudyHertz_%d'));

%% Preform a FEM calculation for each case
for iCases = 1:numCases
    %% Print message
    if strcmp(outMsg,'outputEnabled')
        fprintf('Case number :%d\n', iCases);
    end
    
    %% Parse data from GiD input file
    [strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC, analysis, parameters, ...
        propNLinearAnalysis, propStrDynamics, propGaussInt, propContact] = ...
        parse_StructuralModelFromGid(pathToCase, char(caseName(iCases)), outMsg);

    % Assign values to each iteration step
    propContact.maxIter = maxIter;
    propNBC.tractionLoadVct = tractionLoadVct;

    %% Compute the load vector
    time = 0;
    F = computeLoadVctFEMPlateInMembraneAction...
        (strMsh, analysis, propNBC, time, propGaussInt, '');

    %% Solve the system for the displacement field and the Lagrange Multipliers fields
    [dHat, lambdaHat, nodeIDs_active, numIter, FComplete, minElSize] = ...
        solveSignoriniFrictionlessContact2D...
        (analysis, strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC, bodyForces, ...
        parameters, contactSegments, computeStiffMtxLoadVct, solve_LinearSystem, ...
        propNLinearAnalysis, propContact , propGaussInt, char(caseName(iCases)), ...
        pathToOutput, isUnitTest, outMsg);

    %% Postprocessing
    
    % Compute the contact length, the contact force and the maximum contact
    % pressure
    [contactLength(iCases, 1), contactForce(iCases, 1), maxContactPressure(iCases, 1)] = ...
        computePostprocResultantsSignoriniFrictionlessContact2D...
        (strMsh, parameters, dHat, lambdaHat, nodeIDs_active);

    % Get the length of the contact area and the reaction force on the contact
    appliedForce(iCases, 1) = sum(F);
    hertzContactLength(iCases, 1) = sqrt(4*(2*appliedForce(iCases, 1))*radius*...
        ((1 - parameters.nue^2)/parameters.E)/(pi*parameters.t));
    hertzPressure(iCases, 1) = 2*(2*appliedForce(iCases, 1))/(parameters.t*pi*...
        hertzContactLength(iCases, 1));
    
    % Get the number of mesh elements
    numberOfElements(iCases, 1) = size(strMsh.elements, 1);
end

%% Convergence study

% Plot the maximum contact pressure versus the number of mesh elements
figure('Name','Max Contact Pressure')
hold on
grid on
plot(numberOfElements,maxContactPressure, 'r-o', 'LineWidth', 2);
plot(numberOfElements,hertzPressure, 'b--d','LineWidth', 2);
hold off
set(gca, 'xscale', 'log');
legend('FEM', 'reference', 'location', 'southeast');
xlabel('number of mesh elements');
ylabel('max contact pressure [Pa]');

% Plot the contact length versus the number of mesh elements
figure('Name', 'Contact Length');
hold on;
grid on;
plot(numberOfElements,contactLength, 'r-o', 'LineWidth', 2);
plot(numberOfElements,hertzContactLength, 'b--d', 'LineWidth', 2);
hold off;
set(gca,'xscale', 'log');
legend('FEM', 'reference', 'location', 'southeast');
xlabel('number of mesh elements');
ylabel('contact length [m]');

% Calculate the difference in applied force and resultant contact force
forceDifference = max(abs(appliedForce - contactForce));
fprintf('Maximum difference between applied and computed contact force is : %d\n', ...
    forceDifference);

%% End of the script