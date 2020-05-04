function testCantileverBeamSubjectToDistributedLoad(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the solution of a cantilever beam subject to constant distributed
% load for both the linear Bernoulli and Timoshenko beam analysis.
% Analytical solutions for this problem exist using both analysis types,
% namely,
%
%        Bernoulli :
%             w(x) = p*x^2*(6*L^2 - 4L*x + x^2)/24/E/I
%
%       Timoshenko :
%             w(x) = (p*L*x-p*x^2/2)/G/Aq - (-p*x^4/24+p*L*x^3/6-p*L^2*x^2/4)/E/I
%          beta(x) = - p*x^3/6/E/I
%               Aq = alpha*A 
%
% Function layout :
%
% 0. Read input
%
% 1. Define the geometry using NURBS
%
% 2. Define the material constants
%
% 3. GUI
%
% 4. Define the h- and p- refinement
%
% 5. Apply boundary conditions corresponding to the Benoulli beam analysis
%
% 6. Define the analytical solution corresponding to the Bernoulli beam theory
%
% 7. Define the analytical solution corresponding to the Timoshenko beam theory
%
% 8. Solve the Bernoulli beam problem
%
% 9. Compute the relative error with respect to the analytical solution for the Bernoulli beam analysis
%
% 10. Apply boundary conditions corresponding to the Benoulli beam analysis
%
% 11. Solve the Timoshenko beam problem
%
% 12. Compute the relative error with respect to the analytical solution for the Timoshenko beam analysis
%
% 13. Verify the results
%
%% Function main body

%% 0. Read input

% Define tolerances
absTol = 1e-12;
absTolRelaxed4 = absTol*1e4;

%% 1. Define the geometry using NURBS

% Polynomial degree
p = 1;

% Knot vector
Xi = [0 0 1 1];

% Beam's length
L = 10;

% Control Point coordinates and weights

% x-coordinates
CP(:,1) = [0 1]*L;

% y-coordinates
CP(:,2) = [0 0]*L;      

% z-coordinates
CP(:,3) = [0 0];

% weights
CP(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
for i=length(CP(:,1));
    if CP(i,4)~=1
        isNURBS = 1;
        break;
    end
end

%% 2. Define the material constants

% Young's modulus (connected to epsilon_11 = E)
parameters.EYoung = 4e6;

% Poisson ration
parameters.Nu = 0;

% shear modulus (connected to epsilon_12 = E/(1+nu))
parameters.GShear = parameters.EYoung/(2*(1+parameters.Nu));

% shear correction factor
parameters.alpha = 5/6;

% width of the beam
parameters.b = 1;

% height of the beam
parameters.h = 1;

% cross sectional area
parameters.A = parameters.b*parameters.h;

% Moment of inertia (I_z = I for a simple 2D case)
parameters.I = parameters.b*(parameters.h^3)/12;

% shear cross sectional area
parameters.Aq = parameters.alpha*parameters.A;

%% 3. GUI

% Integration
% .type = 'automatic' : Default choice of Gauss points
% .type = 'manual' : Manual choice of Gauss points
int.type = 'default';
intError.type = 'user';

% Choose linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Number of Gauss points
int.noGP = 6;
int.noGPLoad = 6;
intError.noGP = 12;

% General problem settings
pLoad = - 1e3;
problemSettings.Length = L;
problemSettings.pressure = pLoad;
problemSettings.EYoung = parameters.EYoung;
problemSettings.I = parameters.I;
problemSettings.GShear = parameters.GShear;
problemSettings.Aq = parameters.Aq;

%% 4. Define the h- and p- refinement

% Order elevation
tp = 1;
[Xi,CP,p] = degreeElevateBSplineCurve(p,Xi,CP,tp,'');

% Knot insertion
n = 5;
[Xi,CP] = knotRefineUniformlyBSplineCurve(n,p,Xi,CP,'');

%% 5. Apply boundary conditions corresponding to the Benoulli beam analysis

% Dirichlet boundary conditions

% Clamp the left edge of the beam
homDOFs = [];
xib = [Xi(1) Xi(p+1)]; dir = 1;
homDOFs = findDofsForBernoulliBeams2D(homDOFs,xib,dir,CP);
xib = [Xi(1) Xi(p+2)]; dir = 2;
homDOFs = findDofsForBernoulliBeams2D(homDOFs,xib,dir,CP);

% Neumann boundary conditions
NBC.noCnd = 1;
xib = [0 1];
NBC.xiLoadExtension = {xib};
NBC.loadAmplitude(1,1) = pLoad;
loadDir = 2;
NBC.loadDirection(1,1) = loadDir;
NBC.computeLoadVct = {'computeLoadVctLinePressureVectorForIGABernoulliBeam2D'};

%% 6. Define the analytical solution corresponding to the Bernoulli beam theory

% Define the expected solution in terms of the displacement field
expSolDispBernoulli = [                    0
                                           0
                                           0
                                           0
                                           0
                          -0.487999999999998
                                           0
                          -1.271999999999993
                                           0
                          -2.207999999999986
                                           0
                          -3.199999999999977
                                           0
                          -3.699999999999973];
                      
% Define the expected solution in terms of the complete force vector
expSolFBernoulli = 1e4*[                   0
                          -4.066666666666649
                                           0
                           4.866666666666648
                                           0
                          -0.200000000000000
                                           0
                          -0.200000000000012
                                           0
                          -0.200000000000000
                                           0
                          -0.133333333333349
                                           0
                          -0.066666666666663];
                      
% Define the expected solution in terms of the minimum element edge size
expSolMinElEdgeSizeBernoulli = 1.999999999999999;

% Define the expected solution in terms of the relative error in the 
% solution wrt to the analytical one
expSolDispRelErrL2Bernouli = 0.019076812686924;

%% 7. Define the analytical solution corresponding to the Timoshenko beam theory

% Define the expected solution in terms of the displacement field
expSolDispTimoshenko = [                   0
                                           0
                                           0
                                           0
                          -0.014126449541159
                          -0.144251907859262
                                           0
                          -0.520170804869023
                          -0.336640767688687
                                           0
                          -1.319850756888020
                          -0.443453021657367
                                           0
                          -2.265851909784059
                          -0.489847041933921
                                           0
                          -3.263045396500013
                          -0.501214633445362
                                           0
                          -3.763837050966113
                          -0.499999999999997];
                      
% Define the expected solution in terms of the complete force vector
expSolFTimoshenko = 1e4*[                  0
                           0.933333333333340
                           4.999999999999961
                                           0
                          -0.133333333333339
                          -0.000000000000009
                                           0
                          -0.200000000000001
                          -0.000000000000006
                                           0
                          -0.200000000000043
                          -0.000000000000010
                                           0
                          -0.200000000000044
                          -0.000000000000046
                                           0
                          -0.133333333333355
                          -0.000000000000035
                                           0
                          -0.066666666666558
                          -0.000000000000047];
                      
% Define the expected solution in terms of the minimum element edge size
expSolMinElEdgeSizeTimoshenko = 2.000000000000000;

% Define the expected solution in terms of the relative error in the 
% displacement wrt to the analytical one
expSolRelDisplErrL2Timoshenko = 0.006356040665311;

% Define the expected solution in terms of the relative error in the 
% rotation wrt to the analytical one
expSolRelRotlErrL2Timoshenko = 1.419396956304162;

%% 8. Solve the Bernoulli beam problem
analysis.type = 'Bernoulli';
[dHatBernoulli,FBernoulli,minElEdgeSizeBernoulli] = solve_IGABeamLinear2D...
    (analysis,p,Xi,CP,homDOFs,NBC,parameters,isNURBS,solve_LinearSystem,int,'');

%% 9. Compute the relative error with respect to the analytical solution for the Bernoulli beam analysis
[dispRelErrL2Bernouli,~] = computeErrIGABernoulliBeam2D(p,Xi,CP,isNURBS,...
    dHatBernoulli,@computeExactDispl4BernoulliCantileverBeamInUniformPressure,...
    problemSettings,intError,'');

%% 10. Apply boundary conditions corresponding to the Benoulli beam analysis

% Dirichlet boundary conditions

% Clamp the left edge of the beam (3 DoFs two translations and 1 rotation)
homDOFs = [];
xib = [Xi(1) Xi(p+1)]; dir = 1;
homDOFs = findDofsForTimoshenkoBeams2D(homDOFs,xib,dir,CP);
xib = [Xi(1) Xi(p+1)]; dir = 2;
homDOFs = findDofsForTimoshenkoBeams2D(homDOFs,xib,dir,CP);
xib = [Xi(1) Xi(p+1)]; dir = 3;
homDOFs = findDofsForTimoshenkoBeams2D(homDOFs,xib,dir,CP);

% Neumann boundary conditions
NBC.noCnd = 1;
xib = [0 1];
NBC.xiLoadExtension = {xib};
NBC.loadAmplitude(1,1) = pLoad;
loadDir = 2;
NBC.loadDirection(1,1) = loadDir;
NBC.computeLoadVct = {'computeLoadVctLinePressureVectorForIGATimoshenkoBeam2D'};

%% 11. Solve the Timoshenko beam problem
analysis.type = 'Timoshenko';
[dHatTimoshenko,FTimoshenko,minElEdgeSizeTimoshenko] = solve_IGABeamLinear2D...
    (analysis,p,Xi,CP,homDOFs,NBC,parameters,isNURBS,solve_LinearSystem,int,'');

%% 12. Compute the relative error with respect to the analytical solution for the Timoshenko beam analysis
[relDisplErrL2Timoshenko,relRotErrL2Timoshenko] = computeErrIGATimoshenkoBeam2D...
    (p,Xi,CP,isNURBS,dHatTimoshenko,...
    @computeExactDispl4TimoshenkoCantileverBeamInUniformPressure,...
    problemSettings,intError,'');

%% 13. Verify the results
testCase.verifyEqual(dHatBernoulli,expSolDispBernoulli,'AbsTol',absTol);
testCase.verifyEqual(FBernoulli,expSolFBernoulli,'AbsTol',absTolRelaxed4);
testCase.verifyEqual(minElEdgeSizeBernoulli,expSolMinElEdgeSizeBernoulli,'AbsTol',absTol);
testCase.verifyEqual(minElEdgeSizeBernoulli,expSolMinElEdgeSizeBernoulli,'AbsTol',absTol);
testCase.verifyEqual(dispRelErrL2Bernouli,expSolDispRelErrL2Bernouli,'AbsTol',absTol);
testCase.verifyEqual(dHatTimoshenko,expSolDispTimoshenko,'AbsTol',absTol);
testCase.verifyEqual(FTimoshenko,expSolFTimoshenko,'AbsTol',absTolRelaxed4);
testCase.verifyEqual(minElEdgeSizeTimoshenko,expSolMinElEdgeSizeTimoshenko,'AbsTol',absTol);
testCase.verifyEqual(relDisplErrL2Timoshenko,expSolRelDisplErrL2Timoshenko,'AbsTol',absTol);
testCase.verifyEqual(relRotErrL2Timoshenko,expSolRelRotlErrL2Timoshenko,'AbsTol',absTol);

end
