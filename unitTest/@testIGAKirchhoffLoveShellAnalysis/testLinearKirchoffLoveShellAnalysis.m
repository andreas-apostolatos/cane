function testLinearKirchoffLoveShellAnalysis(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function definition
%
% Test the linear Kirchoff-Love shell isogeometric analysis over the single 
% patch Scordelis-Lo roof benchmark example
%
% Function layout :
%
% 1. Define NURBS parameters
%
% 2. Define material constants
%
% 3. UI
%
% 4. Define h- and p-refinement
%
% 5. Define Dirichlet and Neumann boundary conditions 
%
% 6. Create the B-Spline patch array
%
% 7. Solve the system applying linear analysis
%
% 8. Verify the result
%
%% Function main body

%% 0. Read input

% Define a tolerance for the verification of the results
absTol = 1e-15;
absTolRelaxed = absTol*1.0e+06;

% Define reference solution on the displacements
expSolDisp = [                 0
                               0
                               0
               0.000324616053150
               0.001843253177329
              -0.017879413643329
               0.005445653483464
               0.001192111934897
              -0.030514095388845
               0.010566690913778
               0.001843253177330
              -0.017879413643329
               0.010891306966928
                               0
                               0
               0.008236541510349
                               0
                               0
               0.008070176342455
               0.000058760971878
              -0.015349374046370
               0.005445653483464
              -0.000029018307808
              -0.028612204272582
               0.002821130624472
               0.000058760971878
              -0.015349374046371
               0.002654765456578
                               0
                               0
               0.008236541510349
                               0
                               0
               0.008070176342455
              -0.000058760971878
              -0.015349374046370
               0.005445653483464
               0.000029018307809
              -0.028612204272580
               0.002821130624473
              -0.000058760971878
              -0.015349374046369
               0.002654765456578
                               0
                               0
               0.000000000000000
                               0
                               0
               0.000324616053150
              -0.001843253177329
              -0.017879413643328
               0.005445653483464
              -0.001192111934896
              -0.030514095388842
               0.010566690913778
              -0.001843253177328
              -0.017879413643327
               0.010891306966927
                               0
                               0];
               
% Define the reference solution on the forces
expSolForces = 1.0e+04 * ...
              [0.000000000000052
               3.495824938997198
               2.074953864755903
               0.000000000000001
              -0.000000000000006
              -0.592564162207258
              -0.000000000000024
              -0.000000000000009
              -0.888846243310872
              -0.000000000000023
               0.000000000000007
              -0.592564162207252
              -0.000000000000002
               3.495824938997232
               2.074953864755928
               0.000000000000010
               2.132271765737141
               0.979197782141091
               0.000000000000000
               0.000000000000018
              -1.152665350305322
              -0.000000000000006
               0.000000000000005
              -1.728998025457970
               0.000000000000020
               0.000000000000001
              -1.152665350305312
               0.000000000000009
               2.132271765737080
               0.979197782141083
              -0.000000000000018
              -2.132271765737145
               0.979197782141099
               0.000000000000009
              -0.000000000000019
              -1.152665350305313
              -0.000000000000036
               0.000000000000024
              -1.728998025457986
              -0.000000000000003
               0.000000000000006
              -1.152665350305301
               0.000000000000001
              -2.132271765737174
               0.979197782141098
               0.000000000000004
              -3.495824938997196
               2.074953864755885
               0.000000000000014
               0.000000000000008
              -0.592564162207253
               0.000000000000015
              -0.000000000000015
              -0.888846243310895
               0.000000000000006
              -0.000000000000019
              -0.592564162207250
               0.000000000000012
              -3.495824938997167
               2.074953864755860];
           
% Define the reference solution on the minimum element area size
expSolMinElArea = 50;
   
%% 1. Define NURBS parameters

% Global variables
Length = 50;
Radius = 25;

% Polynomial degrees
p = 1;
q = 2;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP(:, :, 1) = [-Length/2 -Length/2 -Length/2
               Length/2  Length/2  Length/2];
         
% y-coordinates
CP(:, :, 2) = [-Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)
               -Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)];
         
% z-coordinates
CP(:, :, 3) = [Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)
               Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)];
       
% Weights
weight = cos(2*pi/9);
CP(:, :, 4) = [1 weight 1
               1 weight 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = false;
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));
for i = 1:nxi
    for j = 1:neta
        if CP(i, j, 4) ~= 1
            isNURBS = true;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% 2. Define material constants

% Young's modulus
parameters.E = 4.32e8;

% Poisson ratio
parameters.nue = .0;

% Thickness of the shell
parameters.t = .25;

% Density of the shell (used only for dynamics)
parameters.rho = 7850;

%% 3. UI

% On the analysis
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'user';
if strcmp(int.type,'user')
    int.xiNGP = ceil(p + 2);
    int.etaNGP = ceil(q + 1);
    int.xiNGPForLoad = ceil((p + 2)/2);
    int.etaNGPForLoad = ceil((q + 1)/2);
end

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

%% 4. Define h- and p-refinement

% Degree by which to elevate
tp = 1;
tq = 0;
[Xi, Eta, CP, p, q] = degreeElevateBSplineSurface ...
    (p, q, Xi, Eta, CP, tp, tq, '');

% Number of knots to exist in both directions
scaling = 1;
edgeRatio = ceil(Length/Radius/(sin(4*pi/9)));
refXi = edgeRatio*scaling;
refEta = ceil(4/3)*scaling;
[Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, refXi, refEta, '');

%% 5. Define Dirichlet and Neumann boundary conditions 

% supports (Dirichlet boundary conditions)

% back and front curved edges are a rigid diaphragm
homDOFs = [];
xiSup = [0 0];
etaSup = [0 1];
for dirSupp = [2 3]
    homDOFs = findDofs3D(homDOFs, xiSup, etaSup, dirSupp, CP);
end
xiSup = [1 1];
etaSup = [0 1];
for dirSupp = [2 3]
    homDOFs = findDofs3D(homDOFs, xiSup, etaSup, dirSupp, CP);
end

% Fix the back left corner of the shell to avoid rigid body motions
xiSup = [0 0];
etaSup = [0 0];
dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSup, etaSup, dirSupp, CP);

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables = [];

% load (Neuman boundary conditions)
FAmp = - 9e1;
NBC.noCnd = 1;
xib = [0 1];
etab = [0 1];
dirForce = 'z';
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC.isFollower(1, 1) = false;
NBC.isTimeDependent(1, 1) = false;

%% 6. Create the B-Spline patch array
BSplinePatch = fillUpPatch ...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, weakDBC, cables, NBC, [], [], [], [], [], int);

%% 7. Solve the system applying linear analysis
[dHatLinear, F, minElArea] = solve_IGAKirchhoffLoveShellLinear ...
    (BSplinePatch, solve_LinearSystem, '');

%% 8. Verify the result
testCase.verifyEqual(dHatLinear, expSolDisp, 'AbsTol', absTolRelaxed);
testCase.verifyEqual(F, expSolForces, 'AbsTol', absTolRelaxed);
testCase.verifyEqual(minElArea, expSolMinElArea, 'AbsTol', absTol);

end
