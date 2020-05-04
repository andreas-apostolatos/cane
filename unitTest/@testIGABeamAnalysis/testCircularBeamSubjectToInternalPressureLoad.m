function testCircularBeamSubjectToInternalPressureLoad(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the solution to a circular isogeometric beam problem subject to
% internal pressure load. The beam is reduced by applying symmetry boundary
% conditions at its one quarter.
%
% Function layout :
%
% 0. Read input
%
% 1. Define the geometry in terms of NURBS
%
% 2. Define the material constants
%
% 3. GUI
%
% 4. Define the h- and the p- refinement
%
% 5. Define the boundary conditions
%
% 6. Define the expected solution
%
% 7. Solve the Timoshenko beam problem
%
% 8. Verify the results
%
%% Function main body

%% 0. Read input

% Define the absolute tolerances
absTol = 1e-11;
absTolRelaxed6 = absTol*1e6;

%% 1. Define the geometry in terms of NURBS

% Polynomial degree
p = 2;

% Knot vector
Xi = [0 0 0 1 1 1];

% Beam's length
L = 10;

% Control Point coordinates and weights

% x-coordinates
CP(:,1) = [0 1 1]*L;

% y-coordinates
CP(:,2) = [1 1 0]*L;

% z-coordinates
CP(:,3) = [0 0 0];

% weights
CP(:,4) = [1 1/sqrt(2) 1];

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
% .type = 'default' : Default choice of Gauss points
% .type = 'user' : Manual choice of Gauss points
int.type = 'default';

% Number of Gauss points
int.noGP = 6;

% int.ngaussLoad = ceil(0.5*(p+1));
int.noGPLoad = 6;

% Choose linear equation solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

%% 4. Define the h- and the p- refinement

% Order elevation
tp = 1;
[Xi,CP,p] = degreeElevateBSplineCurve(p,Xi,CP,tp,'');

% Knot insertion
n = 11;
[Xi,CP] = knotRefineUniformlyBSplineCurve(n,p,Xi,CP,'');

%% 5. Define the boundary conditions

% Dirichlet boundary conditions
homDOFs = [];

% Fix displacents and rotations accordingly at each end of the beam
    
% Clamp the left edge
xib = [Xi(1) Xi(p+1)]; dir = 1;
homDOFs = findDofsForTimoshenkoBeams2D(homDOFs,xib,dir,CP);
xib = [Xi(1) Xi(p+1)]; dir = 3;
homDOFs = findDofsForTimoshenkoBeams2D(homDOFs,xib,dir,CP);

% Clamp the right edge
xib = [Xi(length(Xi)-p) Xi(length(Xi))]; dir = 2;
homDOFs = findDofsForTimoshenkoBeams2D(homDOFs,xib,dir,CP);
xib = [Xi(length(Xi)-p) Xi(length(Xi))]; dir = 3;
homDOFs = findDofsForTimoshenkoBeams2D(homDOFs,xib,dir,CP);

% Neumann boundary conditions

% On the application of a pressure load on the beam
NBC.noCnd = 1;
xib = [0 1];
NBC.xiLoadExtension = {xib};
pLoad = 1e5;
NBC.loadAmplitude(1,1) = pLoad;
loadDir = 2;
NBC.loadDirection(1,1) = loadDir;
NBC.computeLoadVct = {'computeLoadVctLinePressureVectorForIGATimoshenkoBeam2D'};

%% 6. Define the expected solution

% Define the expected solution in terms of the discrete solution field
expSolDispRotTimoshenko = [0
                           2.500431802405122
                                           0
                           0.108941378790603
                           2.500503458563283
                           0.000156346073462
                           0.332120275119736
                           2.486077061566117
                          -0.000894884726075
                           0.670105437525216
                           2.415624844276520
                          -0.001108852304008
                           1.007357494682567
                           2.294846948844185
                          -0.000633519461504
                           1.332386823512453
                           2.122017539790580
                          -0.000516297579015
                           1.633427640089347
                           1.899590899635212
                          -0.000146319267477
                           1.899590899636696
                           1.633427640088420
                           0.000146319267990
                           2.122017539792662
                           1.332386823511964
                           0.000516297579476
                           2.294846948846820
                           1.007357494682365
                           0.000633519461875
                           2.415624844279582
                           0.670105437525163
                           0.001108852304264
                           2.486077061569448
                           0.332120275119733
                           0.000894884726211
                           2.500503458566675
                           0.108941378790602
                          -0.000156346073420
                           2.500431802408513
                                           0
                                           0];
                               
% Define the expected solution in terms of the complete force vector
expSolFTimoshenko = 1e5*[  -9.990239548610477
                           0.330307253823363
                           0.083768682670374
                           0.052355048826295
                           0.666923578011784
                           0.000000000000002
                           0.165841996610161
                           1.010443248935539
                           0.000000000000001
                           0.378746149803587
                           1.347814300033753
                          -0.000000000000005
                           0.583638019857851
                           1.319527948248408
                          -0.000000000000008
                           0.786062721547095
                           1.246987874704027
                          -0.000000000000000
                           0.972278358104104
                           1.129313050104494
                          -0.000000000000001
                           1.129313050104493
                           0.972278358104104
                           0.000000000000001
                           1.246987874704029
                           0.786062721547094
                          -0.000000000000003
                           1.319527948248415
                           0.583638019857851
                           0.000000000000002
                           1.347814300033759
                           0.378746149803591
                          -0.000000000000001
                           1.010443248935540
                           0.165841996610166
                          -0.000000000000003
                           0.666923578011785
                           0.052355048826294
                          -0.000000000000002
                           0.330307253823373
                          -9.990239548610434
                          -0.083768682670719];
              
% Define the expected solution in terms of the minimum element edge size
expSolMinElEdgeSizeTimoshenko = 1.335822696345342;

%% 7. Solve the Timoshenko beam problem
analysis.type = 'Timoshenko';
[dHatTimoshenko,FTimoshenko,minElEdgeSizeTimoshenko] = solve_IGABeamLinear2D...
    (analysis,p,Xi,CP,homDOFs,NBC,parameters,isNURBS,solve_LinearSystem,int,'');

%% 8. Verify the results
testCase.verifyEqual(dHatTimoshenko,expSolDispRotTimoshenko,'AbsTol',absTol);
testCase.verifyEqual(FTimoshenko,expSolFTimoshenko,'AbsTol',absTolRelaxed6);
testCase.verifyEqual(minElEdgeSizeTimoshenko,expSolMinElEdgeSizeTimoshenko,'AbsTol',absTolRelaxed6);

end
