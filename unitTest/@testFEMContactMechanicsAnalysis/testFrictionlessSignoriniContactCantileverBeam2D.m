function testFrictionlessSignoriniContactCantileverBeam2D(testCase)
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
%        cantilever beam which is fixed at its one end encountering a 
%        circular contact support.
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

% Define tolerances for both cases
absTol = 1e-15;
absTol2 = 1e-15*1e2;
absTol3 = absTol2*1e1;

%% 1. Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/';
caseName = 'unitTest_cantilever_beam';

% Parse the data from the GiD input file
[strMsh, homDBC, inhomDBC, valuesInhomDBC, propNBC, analysis, parameters, ...
    propNLinearAnalysis, ~, propGaussInt, propContact] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'');

%% 2. UI

% On the body forces
bodyForces = @computeConstantVerticalStructureBodyForceVct;

% Amplitude of the externally applied boundary traction
propNBC.tractionLoadVct = [0; - 1e3; 0];

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
center = [2,-4.1];
radius = 4;
startAngle = 3*pi/4;
endAngle = pi/4;
numSegments = 30;
contactSegments = createCircleSegments ...
    (center, radius, startAngle, endAngle, numSegments);

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
expContactLength = 0.101340127927814;

% Expected total contact force
expContactForce = 2.040613795404637e+02;

% Expected solution in terms of the maximum contact pressure
expMaxContactPressure = 1.006814298112094e+03;

% Expected number of contact interations
expNumIter = 18;

% Expected minimum element area size
expMinElSize = 0.073932639612014;

%% 7. Verify the results
testCase.verifyEqual(contactLength, expContactLength, 'AbsTol', absTol);
testCase.verifyEqual(contactForce, expContactForce, 'AbsTol', absTol2);
testCase.verifyEqual(maxContactPressure, expMaxContactPressure, 'AbsTol', absTol3);
testCase.verifyEqual(numIter, expNumIter, 'AbsTol', absTol);
testCase.verifyEqual(minElSize, expMinElSize, 'AbsTol', absTol);

end