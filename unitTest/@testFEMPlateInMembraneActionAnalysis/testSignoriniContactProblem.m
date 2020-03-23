function testSignoriniContactProblem(testCase)
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
% Date : 03.03.2020
%
%% Function layout
%
% 0. Parse data from GiD input file
%
% 1. Define rigid wall- line
%
% 2. Compute the load vector
%
% 3. Solve the system and get the displacement field
%
% 4. Get the length of the contact area and the reaction force on the contact
%
% 5. Define the expected solution
%
% 6. Verify the results
%
%% 0. Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMContactLinearPlateInMembraneAction/';
caseName = 'unitTest_HertzContact';

% Parse the data from the GiD input file
[strMsh,homDBC,~,~,NBC,analysis,parameters,~,~,gaussInt,contactNodes] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'');

%% 1. Define rigid wall- line  [(x0,y0) ; (x1,y1)]

% define bottom contact line segment
wall_1 = [5, -1; 5, 5];
    
% add a wall to the segments of points
segments.points(:,:,1) = wall_1;    

%% 2. Compute the load vector
time = 0;
F = computeLoadVctFEMPlateInMembraneAction(strMsh,NBC,time,gaussInt,'');

%% 3. Solve the system and get the displacement field

maxIteration = 30;

[displacement,lagrange] = solveSignoriniLagrange_1(strMsh,homDBC,contactNodes,F,segments,parameters,analysis,maxIteration,''); 
%[displacement,lagrange] = solveSignoriniLagrange_2(strMsh,homDBC,contactNodes,F,segments,parameters,analysis,maxIteration);

%% 4. Get the length of the contact area and the reaction force on the contact

% Numerical calculation of contact resultants
[contactLength,~,maxContactPressure] = computeContactResultants(strMsh,displacement,lagrange,parameters);

% define variables
radius = 5;
force = sum(F);

% analytical calculation of contact resultants
hertzContactLength = sqrt(4*(2*force)*radius*((1-parameters.nue^2)/parameters.E)/(pi*parameters.t));
hertzPressure = 2*(2*force)/(parameters.t*pi*hertzContactLength);

%% 5. Define the expected solution

% Define the expected computed solutions from previous tests
expComputedPressure = 10726.32519035318;
expComputedContactLength = 0.876276913636778;

%% 6. Verify the results

% Define tolerances for both cases
absTol = 1e-9;
relTol = 1e-1;

% Compare expected values
testCase.verifyEqual(contactLength,expComputedContactLength,'AbsTol',absTol);
testCase.verifyEqual(maxContactPressure,expComputedPressure,'AbsTol',absTol);
testCase.verifyEqual(contactLength,hertzContactLength,'RelTol',relTol);
testCase.verifyEqual(maxContactPressure,hertzPressure,'RelTol',relTol);

end