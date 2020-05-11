function testFrictionlessSignoriniContactWedge2D(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Script documentation
%
% Task : Signorini frictionless contact analysis for the case of 2D
%        wedge which is compressed inside a rigid support.
%
% Date : 30.03.2020
%
%% Function layout
%
% 0. Read input
%
% 1. Parse data from GiD input file
%
% 2. UI
%
% 3. Define boundary segments of the rigid contact wall
%
% 4. Compute normals to segments
%
% 5. Compute the numerical solution in terms of the contact length and maximum contact pressure
%
% 6. Define the expected solution
%
% 7. Verify the results
%
%% 0. Read input

% Define tolerances
absTol = 1e-15;
absTol0 = absTol*0;
absTol2 = absTol*1e3;
absTol7 = absTol2*1e5;

%% 1. Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/';
caseName = 'unitTest_wedge';

% Parse the data from the GiD input file
[strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC, analysis, parameters, ...
    propNLinearAnalysis, ~, propGaussInt, propContact] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'');

%% 2. UI

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

% Maximum number of contact iterations
propContact.maxIter = 50;

% On whether the case is a unit test
isUnitTest = true;

% Path to output
pathToOutput = 'undefined';

%% 3. Define boundary segments of the rigid contact wall

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

%% 4. Compute normals to segments
contactSegments = computeUnitNormalVctsToSegments(contactSegments);

%% 3. Solve the system and get the displacement field
[dHat, lambdaHat, nodeIDs_active, numIter, ~, minElSize] = ...
    solveSignoriniFrictionlessContact2D...
    (analysis, strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC,bodyForces, ...
    parameters, contactSegments, computeStiffMtxLoadVct, solve_LinearSystem, ...
    propNLinearAnalysis, propContact , propGaussInt, caseName, pathToOutput, ...
    isUnitTest, '');

%% 5. Compute the numerical solution in terms of the contact length and maximum contact pressure
[contactLength, contactForce, maxContactPressure] = ...
    computePostprocResultantsSignoriniFrictionlessContact2D ...
    (strMsh, parameters, dHat, lambdaHat, nodeIDs_active);

%% 6. Define the expected solution

% Expected solution in terms of the contact length
expContactLength = 2.484417940196065;

% Expected total contact force
expContactForce = 6.391759601518371e+05;

% Expected solution in terms of the maximum contact pressure
expMaxContactPressure = Inf;

% Expected number of contact interations
expNumIter = 14;

% Expected minimum element area size
expMinElSize = 0.079594036208751;

%% 7. Verify the results
testCase.verifyEqual(contactLength, expContactLength, 'AbsTol', absTol2);
testCase.verifyEqual(contactForce, expContactForce, 'AbsTol', absTol7);
testCase.verifyEqual(maxContactPressure, expMaxContactPressure, 'AbsTol', absTol0);
testCase.verifyEqual(numIter, expNumIter, 'AbsTol', absTol);
testCase.verifyEqual(minElSize, expMinElSize, 'AbsTol', absTol);

end
