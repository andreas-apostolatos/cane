function testNonlinearKirchoffLoveShellMultipatchAnalysis(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function definition
%
% Test the geometrically nonlinear Kirchoff-Love shell isogeometric 
% analysis over the 4-patch semispherical shell
%
% Function layout :
%
% 1. Define the multipatch geometry
%
% 2. Define the material constants
%
% 3. GUI
%
% 4. Refinement
%
% 5. Define boundary conditions
%
% 6. Create the patches and the Lagrange Multiplier fields
%
% 7. Define the expected solutions for the penalty method
%
% 8. Define the expected solutions for the Lagrange Multipliers method
%
% 9. Define parameters related to the nonlinear analysis
%
% 10. Solve the nonlinear coupled system using the penalty method
%
% 11. Solve the nonlinear coupled system using the augmented Lagrange Multipliers method
%
% 12. Verify the solution
%
%% Function main body

%% 0. Read input

% Define tolerances for the check
relaxationFactor = 1e2;
absTol = relaxationFactor*1e-15;
absTolRelaxed4 = relaxationFactor*absTol*1e4;
absTolRelaxed7 = relaxationFactor*absTol*1e7;

%% 1. Define the multipatch geometry

% Global variables
Radius = 0.075000000000000;

% Patch 1 :
% _________

% Polynomial degrees
p1 = 2;
q1 = 2;

% Knot vectors
Xi1 = [0 0 0 1 1 1];
Eta1 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP1(:,:,1) = [0      0      0
              Radius Radius Radius
              Radius Radius Radius];
         
% y-coordinates
CP1(:,:,2) = [-Radius -Radius 0
              -Radius -Radius 0
              0       0       0];
         
% z-coordinates
CP1(:,:,3) = [0 Radius Radius
              0 Radius Radius
              0 0      0];
       
% Weights
weight = sqrt(2)/2;
CP1(:,:,4) = [1      weight  1
              weight weight^2 weight 
              1      weight  1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = 0;
nxi1 = length(CP1(:,1,1));
neta1 = length(CP1(1,:,1));
for i = 1:nxi1
    for j = 1:neta1
        if CP1(i,j,4)~=1
            isNURBS1 = 1;
            break;
        end
    end
    if isNURBS1
        break;
    end
end

% Patch 2 :
% _________

% Polynomial degrees
p2 = 2;
q2 = 2;

% Knot vectors
Xi2 = [0 0 0 1 1 1];
Eta2 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP2(:,:,1) = [0      0      0
              Radius Radius Radius
              Radius Radius Radius];
         
% y-coordinates
CP2(:,:,2) = [Radius Radius 0
              Radius Radius 0
              0      0      0];
         
% z-coordinates
CP2(:,:,3) = [0 Radius Radius
              0 Radius Radius
              0 0      0];
       
% Weights
weight = sqrt(2)/2;
CP2(:,:,4) = [1      weight  1
              weight weight^2 weight 
              1      weight  1];
          
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS2 = 0;
nxi2 = length(CP2(:,1,1));
neta2 = length(CP2(1,:,1));
for i = 1:nxi2
    for j = 1:neta2
        if CP2(i,j,4)~=1
            isNURBS2 = 1;
            break;
        end
    end
    if isNURBS2
        break;
    end
end

% Patch 3 :
% _________

% Polynomial degrees
p3 = 2;
q3 = 2;

% Knot vectors
Xi3 = [0 0 0 1 1 1];
Eta3 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP3(:,:,1) = [0       0       0
              -Radius -Radius -Radius
              -Radius -Radius -Radius];
         
% y-coordinates
CP3(:,:,2) = [Radius Radius 0
              Radius Radius 0
              0      0      0];
         
% z-coordinates
CP3(:,:,3) = [0 Radius Radius
              0 Radius Radius
              0 0      0];
       
% Weights
weight = sqrt(2)/2;
CP3(:,:,4) = [1      weight  1
              weight weight^2 weight 
              1      weight  1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS3 = 0;
nxi3 = length(CP3(:,1,1));
neta3 = length(CP3(1,:,1));
for i = 1:nxi3
    for j = 1:neta3
        if CP3(i,j,4)~=1
            isNURBS3 = 1;
            break;
        end
    end
    if isNURBS3
        break;
    end
end

% Patch 4 :
% _________

% Polynomial degrees
p4 = 2;
q4 = 2;

% Knot vectors
Xi4 = [0 0 0 1 1 1];
Eta4 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP4(:,:,1) = [0       0       0
              -Radius -Radius -Radius
              -Radius -Radius -Radius];
         
% y-coordinates
CP4(:,:,2) = [-Radius -Radius 0
              -Radius -Radius 0
              0       0       0];
         
% z-coordinates
CP4(:,:,3) = [0 Radius Radius
              0 Radius Radius
              0 0      0];
       
% Weights
weight = sqrt(2)/2;
CP4(:,:,4) = [1      weight  1
              weight weight^2 weight 
              1      weight  1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS4 = 0;
nxi4 = length(CP4(:,1,1));
neta4 = length(CP4(1,:,1));
for i = 1:nxi4
    for j = 1:neta4
        if CP4(i,j,4)~=1
            isNURBS4 = 1;
            break;
        end
    end
    if isNURBS4
        break;
    end
end

%%% Coupling between patch 1 and patch 2 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

% Polynomial degree
pLambda12 = 1;

% Knot vector
XiLambda12 = [0 0 1 1];

% Control points weights
CPLambda12(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSLambda12 = 0;
nxiLambda12 = length(CPLambda12(:,1,1));
for i = 1:nxiLambda12
    if CPLambda12(i,4)~=1
        isNURBSLambda12 = 1;
        break;
    end
end

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

% Polynomial degree
pMu12 = 1;

% Knot vector
XiMu12 = [0 0 1 1];

% Control point weights
CPMu12(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu12 = 0;
nxiMu12 = length(CPMu12(:,1,1));
for i = 1:nxiMu12
    if CPMu12(i,4)~=1
        isNURBSMu12 = 1;
        break;
    end
end

%%% Coupling between patch 2 and patch 3 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

% Polynomial degree
pLambda23 = 1;

% Knot vector
XiLambda23 = [0 0 1 1];

% Control points weights
CPLambda23(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSLambda23 = 0;
nxiLambda23 = length(CPLambda23(:,1,1));
for i = 1:nxiLambda23
    if CPLambda23(i,4)~=1
        isNURBSLambda23 = 1;
        break;
    end
end

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

% Polynomial degree
pMu23 = 1;

% Knot vector
XiMu23 = [0 0 1 1];

% Control point weights
CPMu23(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu23 = 0;
nxiMu23 = length(CPMu23(:,1,1));
for i = 1:nxiMu23
    if CPMu23(i,4)~=1
        isNURBSMu23 = 1;
        break;
    end
end

%%% Coupling between patch 3 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

% Polynomial degree
pLambda34 = 1;

% Knot vector
XiLambda34 = [0 0 1 1];

% Control points weights
CPLambda34(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSLambda34 = 0;
nxiLambda34 = length(CPLambda34(:,1,1));
for i = 1:nxiLambda34
    if CPLambda34(i,4)~=1
        isNURBSLambda34 = 1;
        break;
    end
end

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

% Polynomial degree
pMu34 = 1;

% Knot vector
XiMu34 = [0 0 1 1];

% Control point weights
CPMu34(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu34 = 0;
nxiMu34 = length(CPMu34(:,1,1));
for i = 1:nxiMu34
    if CPMu34(i,4)~=1
        isNURBSMu34 = 1;
        break;
    end
end

%%% Coupling between patch 1 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

% Polynomial degree
pLambda14 = 1;

% Knot vector
XiLambda14 = [0 0 1 1];

% Control points weights
CPLambda14(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSLambda14 = 0;
nxiLambda14 = length(CPLambda14(:,1,1));
for i = 1:nxiLambda14
    if CPLambda14(i,4)~=1
        isNURBSLambda14 = 1;
        break;
    end
end

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

% Polynomial degree
pMu14 = 1;

% Knot vector
XiMu14 = [0 0 1 1];

% Control point weights
CPMu14(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu14 = 0;
nxiMu14 = length(CPMu14(:,1,1));
for i = 1:nxiMu14
    if CPMu14(i,4)~=1
        isNURBSMu14 = 1;
        break;
    end
end

%% 2. Define the material constants

% general parameters
EYoung = 2.1e6;
nue = 0.0;
thickness = .03;

% Patch 1 :
% _________

% Young's modulus
parameters1.E = EYoung;

% Poisson ratio
parameters1.nue = nue;

% Thickness of the plate
parameters1.t = thickness;

% Patch 2 :
% _________

% Young's modulus
parameters2.E = EYoung;

% Poisson ratio
parameters2.nue = nue;

% Thickness of the plate
parameters2.t = thickness;

% Patch 3 :
% _________

% Young's modulus
parameters3.E = EYoung;

% Poisson ratio
parameters3.nue = nue;

% Thickness of the plate
parameters3.t = .03;

% Patch 4 :
% _________

% Young's modulus
parameters4.E = EYoung;

% Poisson ratio
parameters4.nue = nue;

% Thickness of the plate
parameters4.t = thickness;

%% 3. GUI

% Analysis type
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

% Define linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points

% Patch 1 :
% _________

int1.type = 'user';
if strcmp(int1.type,'user')
    int1.xiNGP = p1 + 1;
    int1.etaNGP = q1 + 1;
    int1.xiNGPForLoad = ceil((p1 + 1)/2);
    int1.etaNGPForLoad = ceil((q1 + 1)/2);
end

% Patch 2 :
% _________

int2.type = 'user';
if strcmp(int2.type,'user')
    int2.xiNGP = p2 + 1;
    int2.etaNGP = q2 + 1;
    int2.xiNGPForLoad = ceil((p2 + 1)/2);
    int2.etaNGPForLoad = ceil((q2 + 1)/2);
end

% Patch 3 :
% _________

int3.type = 'user';
if strcmp(int3.type,'user')
    int3.xiNGP = p3 + 1;
    int3.etaNGP = q3 + 1;
    int3.xiNGPForLoad = ceil((p3 + 1)/2);
    int3.etaNGPForLoad = ceil((q3 + 1)/2);
end

% Patch 2 :
% _________

int4.type = 'user';
if strcmp(int4.type,'user')
    int4.xiNGP = p4 + 1;
    int4.etaNGP = q4 + 1;
    int4.xiNGPForLoad = ceil((p4 + 1)/2);
    int4.etaNGPForLoad = ceil((q4 + 1)/2);
end

% Interface integration :
% _______________________

intC.type = 'default';
intC.method = 'Nitsche';
if strcmp(intC.type,'user')
    if strcmp(intC.method,'lagrangeMultipliers')
        intC.noGP1 = 12;
        intC.nGP2 = 12;
    else
        intC.noGPs = 12;
    end
    intC.nGPError = 12;
end

% On the coupling

% Compute the material matrices for the membrane and the bending part
Dm = EYoung*thickness/(1-nue^2)*...
      [1   nue  0
       nue 1    0
       0   0    (1-nue)/2];
Db = thickness^3/12*Dm;

%% 4. Refinement

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

tp1 = 0;
tq1 = 0;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% Patch 2 :
% _________

tp2 = 0;
tq2 = 0;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');

% Patch 3 :
% _________

tp3 = 0;
tq3 = 0;
[Xi3,Eta3,CP3,p3,q3] = degreeElevateBSplineSurface(p3,q3,Xi3,Eta3,CP3,tp3,tq3,'');

% Patch 2 :
% _________

tp4 = 0;
tq4 = 0;
[Xi4,Eta4,CP4,p4,q4] = degreeElevateBSplineSurface(p4,q4,Xi4,Eta4,CP4,tp4,tq4,'');

%%% Coupling between patch 1 and patch 2 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

tpLambda = 0;
[XiLambda12,CPLambda12,pLambda12] = degreeElevateBSplineCurve(pLambda12,XiLambda12,CPLambda12,tpLambda,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

tpMu = 0;
[XiMu12,CPMu12,pMu12] = degreeElevateBSplineCurve(pMu12,XiMu12,CPMu12,tpMu,'');

%%% Coupling between patch 2 and patch 3 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

tpLambda = 0;
[XiLambda23,CPLambda23,pLambda23] = degreeElevateBSplineCurve(pLambda23,XiLambda23,CPLambda23,tpLambda,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

tpMu = 0;
[XiMu23,CPMu23,pMu23] = degreeElevateBSplineCurve(pMu23,XiMu23,CPMu23,tpMu,'');

%%% Coupling between patch 3 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

tpLambda = 0;
[XiLambda34,CPLambda34,pLambda34] = degreeElevateBSplineCurve(pLambda34,XiLambda34,CPLambda34,tpLambda,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

tpMu = 0;
[XiMu34,CPMu34,pMu34] = degreeElevateBSplineCurve(pMu34,XiMu34,CPMu34,tpMu,'');

%%% Coupling between patch 1 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

tpLambda = 0;
[XiLambda14,CPLambda14,pLambda14] = degreeElevateBSplineCurve(pLambda14,XiLambda14,CPLambda14,tpLambda,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

tpMu = 0;
[XiMu14,CPMu14,pMu14] = degreeElevateBSplineCurve(pMu14,XiMu14,CPMu14,tpMu,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

refXi1 = 2;
refEta1 = 2;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface(p1,Xi1,q1,Eta1,CP1,refXi1,refEta1,'');

% Patch 2 :
% _________

refXi2 = 5;
refEta2 = 5;
[Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface(p2,Xi2,q2,Eta2,CP2,refXi2,refEta2,'');

% Patch 3 :
% _________

refXi3 = 4;
refEta3 = 4;
[Xi3,Eta3,CP3] = knotRefineUniformlyBSplineSurface(p3,Xi3,q3,Eta3,CP3,refXi3,refEta3,'');

% Patch 4 :
% _________

refXi4 = 3;
refEta4 = 3;
[Xi4,Eta4,CP4] = knotRefineUniformlyBSplineSurface(p4,Xi4,q4,Eta4,CP4,refXi4,refEta4,'');

%%% Coupling between patch 1 and patch 2 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

nLambda12 = min(refEta1,refEta2);
[XiLambda12,CPLambda12] = knotRefineUniformlyBSplineCurve(nLambda12,pLambda12,XiLambda12,CPLambda12,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

nMu12 = ceil(nLambda12/2);
[XiMu12,CPMu12] = knotRefineUniformlyBSplineCurve(nMu12,pMu12,XiMu12,CPMu12,'');

%%% Coupling between patch 2 and patch 3 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

nLambda23 = min(refEta2,refEta3);
[XiLambda23,CPLambda23] = knotRefineUniformlyBSplineCurve(nLambda23,pLambda23,XiLambda23,CPLambda23,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

nMu23 = ceil(nLambda23/2);
[XiMu23,CPMu23] = knotRefineUniformlyBSplineCurve(nMu23,pMu23,XiMu23,CPMu23,'');

%%% Coupling between patch 3 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

nLambda34 = min(refEta3,refEta4);
[XiLambda34,CPLambda34] = knotRefineUniformlyBSplineCurve(nLambda34,pLambda34,XiLambda34,CPLambda34,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

nMu34 = ceil(nLambda34/2);
[XiMu34,CPMu34] = knotRefineUniformlyBSplineCurve(nMu34,pMu34,XiMu34,CPMu34,'');

%%% Coupling between patch 1 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

nLambda14 = min(refEta1,refEta4);
[XiLambda14,CPLambda14] = knotRefineUniformlyBSplineCurve(nLambda14,pLambda14,XiLambda14,CPLambda14,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

nMu14 = ceil(nLambda14/2);
[XiMu14,CPMu14] = knotRefineUniformlyBSplineCurve(nMu14,pMu14,XiMu14,CPMu14,'');

%% 5. Define boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs1 = [];
xisup1 = [0 1];   etasup1 = [0 0];
for dir = 1:3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak Dirichlet boundary conditions
weakDBC1 = [];

% Embedded cables
cables1 = [];

% Patch 2 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs2 = [];
xisup2 = [0 1];   etasup2 = [0 0];
for dir = 1:3
    homDOFs2 = findDofs3D...
        (homDOFs2,xisup2,etasup2,dir,CP2);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs2 = [];
valuesInhomDOFs2 = [];

% Weak Dirichlet boundary conditions
weakDBC2 = [];

% Embedded cables
cables2 = [];

% Patch 3 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs3 = [];
xisup3 = [0 1];   etasup3 = [0 0];
for dir = 1:3
    homDOFs3 = findDofs3D...
        (homDOFs3,xisup3,etasup3,dir,CP3);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs3 = [];
valuesInhomDOFs3 = [];

% Weak Dirichlet boundary conditions
weakDBC3 = [];

% Embedded cables
cables3 = [];

% Patch 4 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs4 = [];
xisup4 = [0 1];   etasup4 = [0 0];
for dir = 1:3
    homDOFs4 = findDofs3D...
        (homDOFs4,xisup4,etasup4,dir,CP4);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs4 = [];
valuesInhomDOFs4 = [];

% Weak Dirichlet boundary conditions
weakDBC4 = [];

% Embedded cables
cables4 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parameter
loadAmplitude = 1e5;

% Patch 1 :
% _________

FAmp1 = loadAmplitude;
NBC1.noCnd = 1;
xib1 = [0 1];   etab1 = [0 1];   dirForce1 = 'normal';
NBC1.xiLoadExtension = {xib1};
NBC1.etaLoadExtension = {etab1};
NBC1.loadAmplitude = {FAmp1};
NBC1.loadDirection = {dirForce1};
NBC1.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC1.isFollower(1,1) = true;
NBC1.isTimeDependent(1,1) = false;

% Patch 2 :
% _________

FAmp2 = - loadAmplitude;
NBC2.noCnd = 1;
xib2 = [0 1];   etab2 = [0 1];   dirForce2 = 'normal';
NBC2.xiLoadExtension = {xib2};
NBC2.etaLoadExtension = {etab2};
NBC2.loadAmplitude = {FAmp2};
NBC2.loadDirection = {dirForce2};
NBC2.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC2.isFollower(1,1) = true;
NBC2.isTimeDependent(1,1) = false;

% Patch 3 :
% _________

FAmp3 = loadAmplitude;
NBC3.noCnd = 1;
xib3 = [0 1];   etab3 = [0 1];   dirForce3 = 'normal';
NBC3.xiLoadExtension = {xib3};
NBC3.etaLoadExtension = {etab3};
NBC3.loadAmplitude = {FAmp3};
NBC3.loadDirection = {dirForce3};
NBC3.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC3.isFollower(1,1) = true;
NBC3.isTimeDependent(1,1) = false;

% Patch 4 :
% _________

FAmp4 = - loadAmplitude;
NBC4.noCnd = 1;
xib4 = [0 1];   etab4 = [0 1];   dirForce4 = 'normal';
NBC4.xiLoadExtension = {xib4};
NBC4.etaLoadExtension = {etab4};
NBC4.loadAmplitude = {FAmp4};
NBC4.loadDirection = {dirForce4};
NBC4.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC4.isFollower(1,1) = true;
NBC4.isTimeDependent(1,1) = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface parametrizations     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Connection with patch 2:
xicoup12 = [0 1];   etacoup12 = [1 1];

% Connection with patch 4:
xicoup14 = [0 0];   etacoup14 = [0 1];

% Patch 2 :
% _________

% Connection with patch 1:
xicoup21 = [0 1];   etacoup21 = [1 1];

% Connection with patch 3:
xicoup23 = [0 0];   etacoup23 = [0 1];

% Patch 3 :
% _________

% Connection with patch 2:
xicoup32 = [0 0];   etacoup32 = [0 1];

% Connection with patch 4:
xicoup34 = [0 1];   etacoup34 = [1 1];

% Patch 4 :
% _________

% Connection with patch 3:
xicoup43 = [0 1];   etacoup43 = [1 1];

% Connection with patch 1:
xicoup41 = [0 0];   etacoup41 = [0 1];

% Define connections :
% ____________________

% Number of connections
noConnections = 4;

% Define connections by patch numbers
connections.No = noConnections;
connections.xiEtaCoup = zeros(noConnections,10);
connections.xiEtaCoup(1,:) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21];
connections.xiEtaCoup(2,:) = [2 3 xicoup23 etacoup23 xicoup32 etacoup32];
connections.xiEtaCoup(3,:) = [3 4 xicoup34 etacoup34 xicoup43 etacoup43];
connections.xiEtaCoup(4,:) = [1 4 xicoup14 etacoup14 xicoup41 etacoup41];

%% 6. Create the patches and the Lagrange Multiplier fields

% Patch 1 :
% _________

patch1 = fillUpPatch(analysis,p1,Xi1,q1,Eta1,CP1,isNURBS1,parameters1,homDOFs1,....
    inhomDOFs1,valuesInhomDOFs1,weakDBC1,cables1,NBC1,[],[],[],[],[],int1);

% Patch 2 :
% _________

patch2 = fillUpPatch(analysis,p2,Xi2,q2,Eta2,CP2,isNURBS2,parameters2,homDOFs2,...
    inhomDOFs2,valuesInhomDOFs2,weakDBC2,cables2,NBC2,[],[],[],[],[],int2);

% Patch 3 :
% _________

patch3 = fillUpPatch(analysis,p3,Xi3,q3,Eta3,CP3,isNURBS3,parameters3,homDOFs3,...
    inhomDOFs3,valuesInhomDOFs3,weakDBC3,cables3,NBC3,[],[],[],[],[],int3);

% Patch 4 :
% _________

patch4 = fillUpPatch(analysis,p4,Xi4,q4,Eta4,CP4,isNURBS4,parameters4,homDOFs4,...
    inhomDOFs4,valuesInhomDOFs4,weakDBC4,cables4,NBC4,[],[],[],[],[],int4);

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch1 patch2 patch3 patch4};

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

lambda12 = fillUpLagrangeMultipliers(pLambda12,XiLambda12,CPLambda12,isNURBSLambda12);
lambda23 = fillUpLagrangeMultipliers(pLambda23,XiLambda23,CPLambda23,isNURBSLambda23);
lambda34 = fillUpLagrangeMultipliers(pLambda34,XiLambda34,CPLambda34,isNURBSLambda34);
lambda14 = fillUpLagrangeMultipliers(pLambda14,XiLambda14,CPLambda14,isNURBSLambda14);
connections.lambda = {lambda12 lambda23 lambda34 lambda14};

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

mu12 = fillUpLagrangeMultipliers(pMu12,XiMu12,CPMu12,isNURBSMu12);
mu23 = fillUpLagrangeMultipliers(pMu23,XiMu23,CPMu23,isNURBSMu23);
mu34 = fillUpLagrangeMultipliers(pMu34,XiMu34,CPMu34,isNURBSMu34);
mu14 = fillUpLagrangeMultipliers(pMu14,XiMu14,CPMu14,isNURBSMu14);
connections.mu = {mu12 mu23 mu34 mu14};

%% 7. Define the expected solutions for the penalty method

% Defined the expected solution in terms of the Control Point Displacements
expSolDispPenalty = [  0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
      -0.000007638609990
      -0.003042503089102
       0.002048039431318
       0.001578214812278
      -0.003148849501743
       0.001724802971797
       0.001473221080284
      -0.000262260775442
       0.000663407515352
      -0.000014300388260
       0.000090323748985
       0.000002262184116
       0.000009927280605
      -0.002190227732446
       0.005979988533284
       0.002061232875911
      -0.002015753542842
       0.005851767340769
       0.003329776788432
      -0.000480084417598
       0.001962805419541
       0.000004473011057
       0.000118202027670
      -0.000039391158678
       0.000020934945755
       0.000023944711226
       0.006071423213423
       0.002087281845421
       0.000023816752789
       0.005979456342863
       0.003274255174133
      -0.000020614505244
       0.002077337535140
       0.000015409300458
       0.000011355076435
      -0.000015826141894
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
       0.000005975636372
       0.001099650195613
       0.000663821115353
       0.000168494430948
       0.001112737143957
       0.000575826605946
       0.000496479217367
       0.000980219248284
       0.000509310457897
       0.000634125985517
       0.000613785695410
       0.000399430430457
       0.000483920857986
       0.000216643154960
       0.000255443776007
       0.000170438018473
      -0.000004583627155
       0.000093903795205
       0.000001429554129
      -0.000000530911628
       0.000003234849559
       0.000024655551579
       0.002593813223589
       0.002354306546575
       0.000488574507185
       0.002642904474926
       0.002268751970155
       0.001323826346388
       0.002376456648657
       0.002014289312001
       0.001749328173809
       0.001558477848248
       0.001503040259743
       0.001399902423135
       0.000585033362127
       0.000860959026062
       0.000512632290150
      -0.000006923093350
       0.000279030870584
       0.000003872039221
      -0.000000937793767
       0.000010011878276
       0.000027087992279
       0.002856992507254
       0.004098983231617
       0.000651016031353
       0.002847101412370
       0.004058408949370
       0.001794819570535
       0.002575572134884
       0.003651314008788
       0.002465170082587
       0.001778355981145
       0.002730188477462
       0.002112858614614
       0.000725886536793
       0.001513302485195
       0.000827170267633
       0.000000826252004
       0.000452616589963
       0.000005447513738
      -0.000001328109932
       0.000016062057221
       0.000027844201460
       0.002042339905632
       0.005388709362656
       0.000702975550375
       0.002019169705250
       0.005313132352809
       0.001959908363356
       0.001819846364756
       0.004793629352848
       0.002767676421458
       0.001300769920568
       0.003638435354368
       0.002523044191083
       0.000575000228406
       0.002045134667910
       0.001059428703782
       0.000010877966192
       0.000599287134433
       0.000006159978015
      -0.000001075090753
       0.000015500042005
       0.000032613718486
       0.000748942452865
       0.005995973972903
       0.000710503462972
       0.000764881228202
       0.005909143829947
       0.001977103793715
       0.000690201110724
       0.005313443447573
       0.002805901785155
       0.000499565316963
       0.004057712309117
       0.002671313895147
       0.000238161198116
       0.002347125211425
       0.001176051482370
       0.000013981193513
       0.000706685302322
       0.000007343847846
       0.000002637712134
       0.000014548826111
       0.000019651395014
       0.000081683599145
       0.006076518204528
       0.000776055832970
       0.000091977489522
       0.006025330850915
       0.001971200366977
       0.000082086471205
       0.005389405013482
       0.002807790595752
       0.000060346307253
       0.004077711371213
       0.002658986619796
       0.000035528176739
       0.002394549723747
       0.001180525328510
       0.000015630577986
       0.000737784377111
       0.000006906757787
       0.000004438217177
      -0.000002397010160
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
      -0.000008772219104
       0.001414597568851
       0.000846544213938
      -0.000284418248028
       0.001448181779739
       0.000722103314558
      -0.000759955418836
       0.001125102347955
       0.000605133733422
      -0.000742444436815
       0.000442873490468
       0.000398944581782
      -0.000273316506032
      -0.000005594422457
       0.000144943786379
      -0.000002378882765
      -0.000003414602505
       0.000006317810273
      -0.000024674127640
       0.002946117534143
       0.003049393165292
      -0.000730604441601
       0.002991302775926
       0.002979904142657
      -0.001904030122778
       0.002437502122870
       0.002467829203075
      -0.002028370121421
       0.001071008753032
       0.001435193760932
      -0.000812756315148
       0.000003234653716
       0.000416499382537
      -0.000005941642167
      -0.000008267025005
       0.000022210885816
      -0.000039221191421
       0.002489972956132
       0.005027673801717
      -0.000882641062622
       0.002453422408846
       0.004967990472791
      -0.002374917458132
       0.002046442594735
       0.004184060369836
      -0.002755685009451
       0.001013224493366
       0.002457609630324
      -0.001256398630647
       0.000027315946027
       0.000650702771184
      -0.000007933139200
      -0.000010122778198
       0.000036762453771
      -0.000032324377812
       0.000938459930676
       0.005987224985405
      -0.000911364511434
       0.000924964811321
       0.005862551598484
      -0.002439967021808
       0.000778209505543
       0.004970086387195
      -0.002945114848887
       0.000423543133835
       0.003012817773997
      -0.001488336035584
       0.000023069761911
       0.000812492420048
      -0.000010674050883
      -0.000002908331160
       0.000051771304466
      -0.000049911486441
       0.000072125469628
       0.006070144841010
      -0.000905048557043
       0.000049646361865
       0.005997417900833
      -0.002496796631402
       0.000060454224240
       0.005012844341552
      -0.002900604696204
       0.000048270798919
       0.003069709223999
      -0.001496892514295
       0.000012592266609
       0.000838556789515
      -0.000011133659967
       0.000001421577735
       0.000048081408311
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
      -0.000031087848774
      -0.001959485964026
       0.001263849977379
      -0.000581272897285
      -0.002043756187787
       0.000997717569295
      -0.001197056622342
      -0.001115626952136
       0.000720735097759
      -0.000528612256965
       0.000022956342928
       0.000264191842090
      -0.000001925996158
       0.000004220615809
       0.000009097133198
      -0.000062792393481
      -0.002975804904194
       0.004149202695872
      -0.001204573853604
      -0.003059434215118
       0.004176533282383
      -0.002780414352553
      -0.001954731726037
       0.002841293563189
      -0.001507811395346
      -0.000047458925732
       0.000727749629683
      -0.000003162498096
       0.000024328124268
       0.000029778337673
      -0.000053379462805
      -0.001306367403613
       0.006032320528022
      -0.001233524477654
      -0.001182585209176
       0.005853622143460
      -0.003082975924514
      -0.000819378283682
       0.004195452635505
      -0.002053736788907
      -0.000067074716847
       0.001103296957258
      -0.000008925944951
       0.000020187769168
       0.000041136074713
      -0.000042146919590
      -0.000003126200209
       0.006064704910330
      -0.001260634637054
      -0.000015578987124
       0.005987530588905
      -0.003079292851636
      -0.000000003354142
       0.004215055354100
      -0.002045280667765
      -0.000011261196043
       0.001175620491581
      -0.000013185036205
       0.000005175922995
       0.000034018776760];

% Define the expected solution in terms of the deformation
expSolCPHistPenalty = struct([]);
for counterPatches = 1:4
    for counterComp = 1:4
        expSolCPHistPenalty{counterPatches}(:,:,counterComp,1) = ...
            BSplinePatches{counterPatches}.CP(:,:,counterComp);
    end
end

% Patch 1 :
% _________

expSolCPHistPenalty{1}(:,:,1,2) = [0  -0.000007638609990   0.000009927280605   0.000020934945755
                                   0.031066017177982   0.032644231990261   0.033127250053893   0.033153299023403
                                   0.075000000000000   0.076473221080284   0.078329776788432   0.078274255174133
                                   0.075000000000000   0.074985699611740   0.075004473011057   0.075015409300458];
expSolCPHistPenalty{1}(:,:,2,2) = [   -0.075000000000000  -0.078042503089102  -0.033256244910428   0.000023944711226
                                      -0.075000000000000  -0.078148849501743  -0.033081770720824   0.000023816752789
                                      -0.031066017177982  -0.031328277953424  -0.013348050061634  -0.000020614505244
                                                       0   0.000090323748985   0.000118202027670   0.000011355076435];
expSolCPHistPenalty{1}(:,:,3,2) = [0   0.033114056609301   0.080979988533284   0.081071423213423
                                   0   0.032790820149779   0.080851767340769   0.080979456342863
                                   0   0.013531373159388   0.033028822597523   0.033143354713122
                                   0   0.000002262184116  -0.000039391158678  -0.000015826141894];
expSolCPHistPenalty{1}(:,:,4,2) = [1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000
                                   0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                                   0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                                   1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000];

% Patch 2 :
% _________

expSolCPHistPenalty{2}(:,:,1,2) = [                0   0.000005975636372   0.000024655551579   0.000027087992279   0.000027844201460   0.000032613718486   0.000019651395014
                                   0.011266582861675   0.011435077292624   0.011755157368860   0.011917598893028   0.011969558412050   0.011977086324647   0.012042638694646
                                   0.033674222389275   0.034170701606642   0.034998048735664   0.035469041959811   0.035634130752631   0.035651326182990   0.035645422756253
                                   0.053765960839350   0.054400086824867   0.055515289013159   0.056231130921937   0.056533637260807   0.056571862624504   0.056573751435101
                                   0.068112370398213   0.068596291256198   0.069512272821347   0.070225229012826   0.070635414589295   0.070783684293359   0.070771357018008
                                   0.075000000000000   0.075170438018473   0.075512632290150   0.075827170267633   0.076059428703782   0.076176051482370   0.076180525328510
                                   0.075000000000000   0.075001429554129   0.075003872039221   0.075005447513738   0.075006159978015   0.075007343847846   0.075006906757787];
expSolCPHistPenalty{2}(:,:,2,2) = [0.075000000000000   0.076099650195613   0.070706183621802   0.056622953346604   0.035716562294908   0.012015525314540   0.000081683599145
                                   0.075000000000000   0.076112737143957   0.070755274873138   0.056613062251719   0.035693392094526   0.012031464089877   0.000091977489522
                                   0.068112370398213   0.069092589646497   0.064233723332168   0.051403932661625   0.032401594474757   0.010922116644011   0.000082086471205
                                   0.053765960839350   0.054379746534760   0.050386838374989   0.040322069914191   0.025441128884266   0.008576347356085   0.000060346307253
                                   0.033674222389275   0.033890865544236   0.031166781472129   0.024866245500491   0.015694376942038   0.005296740089466   0.000035528176739
                                   0.011266582861675   0.011261999234520   0.010224992439938   0.008077608291126   0.005069456857542   0.001706459718566   0.000015630577986
                                                   0  -0.000000530911628  -0.000000937793767  -0.000001328109932  -0.000001075090753   0.000002637712134   0.000004438217177];
expSolCPHistPenalty{2}(:,:,3,2) = [0   0.011930403977028   0.036028528935850   0.057864944070967   0.073501079760869   0.080995973972903   0.081076518204528
                                   0   0.011842409467621   0.035942974359430   0.057824369788719   0.073425502751022   0.080909143829947   0.081025330850915
                                   0   0.010741225991185   0.032596037422002   0.052479674535529   0.066650896036358   0.073425813845786   0.073501775411695
                                   0   0.008476212469578   0.025643399223441   0.041273902410508   0.052466795881109   0.057823673148467   0.057843672210563
                                   0   0.005314022667357   0.015980335739693   0.025653661448893   0.032626882777911   0.036021347600701   0.036068772113023
                                   0   0.001786382320259   0.005337609761934   0.008529398629085   0.010831202667721   0.011973268163997   0.012004367238786
                                   0   0.000003234849559   0.000010011878276   0.000016062057221   0.000015500042005   0.000014548826111  -0.000002397010160];
expSolCPHistPenalty{2}(:,:,4,2) = [1.000000000000000   0.941421356237309   0.871126983722081   0.847695526217005   0.871126983722081   0.941421356237309   1.000000000000000
                                   0.941421356237309   0.886274169979695   0.820097546470558   0.798038671967512   0.820097546470558   0.886274169979695   0.941421356237309
                                   0.871126983722081   0.820097546470558   0.758862221768731   0.738450446868122   0.758862221768731   0.820097546470558   0.871126983722081
                                   0.847695526217005   0.798038671967513   0.738450446868122   0.718587705168325   0.738450446868122   0.798038671967513   0.847695526217005
                                   0.871126983722081   0.820097546470558   0.758862221768731   0.738450446868122   0.758862221768731   0.820097546470558   0.871126983722081
                                   0.941421356237309   0.886274169979695   0.820097546470558   0.798038671967512   0.820097546470558   0.886274169979695   0.941421356237309
                                   1.000000000000000   0.941421356237309   0.871126983722081   0.847695526217005   0.871126983722081   0.941421356237309   1.000000000000000];

% Patch 3 :
% _________

expSolCPHistPenalty{3}(:,:,1,2) = [                    0  -0.000008772219104  -0.000024674127640  -0.000039221191421  -0.000032324377812  -0.000049911486441
                                      -0.014305767737291  -0.014590185985319  -0.015036372178892  -0.015188408799913  -0.015217132248724  -0.015210816294334
                                      -0.042049512883487  -0.042809468302323  -0.043953543006265  -0.044424430341618  -0.044489479905295  -0.044546309514889
                                      -0.064016504294496  -0.064758948731310  -0.066044874415917  -0.066772189303947  -0.066961619143382  -0.066917108990700
                                      -0.075000000000000  -0.075273316506032  -0.075812756315148  -0.076256398630647  -0.076488336035584  -0.076496892514295
                                      -0.075000000000000  -0.075002378882765  -0.075005941642167  -0.075007933139200  -0.075010674050883  -0.075011133659967];
expSolCPHistPenalty{3}(:,:,2,2) = [    0.075000000000000   0.076414597568851   0.066962621828639   0.044539485839619   0.015244227667967   0.000072125469628
                                       0.075000000000000   0.076448181779739   0.067007807070422   0.044502935292333   0.015230732548612   0.000049646361865
                                       0.064016504294496   0.065141606642451   0.057079006417366   0.037937946889231   0.012988946062747   0.000060454224240
                                       0.042049512883487   0.042492386373954   0.036962513047527   0.024588711609880   0.008444217330866   0.000048270798919
                                       0.014305767737291   0.014300173314834   0.012213971210921   0.008047990143059   0.002751802969289   0.000012592266609
                                                       0  -0.000003414602505  -0.000008267025005  -0.000010122778198  -0.000002908331160   0.000001421577735];
expSolCPHistPenalty{3}(:,:,3,2) = [0   0.015152311951229   0.045098906048778   0.069044178096213   0.080987224985405   0.081070144841010
                                   0   0.015027871051849   0.045029417026143   0.068984494767286   0.080862551598484   0.080997417900833
                                   0   0.012815870290627   0.038359333497570   0.058825564664332   0.068986590681691   0.069029348636048
                                   0   0.008419618778814   0.025010680877445   0.038349113924819   0.045062330657483   0.045119222107486
                                   0   0.002873676993757   0.008437173579569   0.012861439328389   0.015118260157339   0.015144324526806
                                   0   0.000006317810273   0.000022210885816   0.000036762453771   0.000051771304466   0.000048081408311];
expSolCPHistPenalty{3}(:,:,4,2) = [1.000000000000000   0.926776695296637   0.853553390593274   0.853553390593274   0.926776695296637   1.000000000000000
                                   0.926776695296637   0.858915042944955   0.791053390593274   0.791053390593274   0.858915042944955   0.926776695296637
                                   0.853553390593274   0.791053390593274   0.728553390593274   0.728553390593274   0.791053390593274   0.853553390593274
                                   0.853553390593274   0.791053390593274   0.728553390593274   0.728553390593274   0.791053390593274   0.853553390593274
                                   0.926776695296637   0.858915042944955   0.791053390593274   0.791053390593274   0.858915042944955   0.926776695296637
                                   1.000000000000000   0.926776695296637   0.853553390593274   0.853553390593274   0.926776695296637   1.000000000000000];

% Patch 4 :
% _________                           
    
expSolCPHistPenalty{4}(:,:,1,2) = [                    0  -0.000031087848774  -0.000062792393481  -0.000053379462805  -0.000042146919590
                                      -0.019590290622281  -0.020171563519566  -0.020794864475885  -0.020823815099934  -0.020850925259335
                                      -0.055094310254260  -0.056291366876602  -0.057874724606813  -0.058177286178774  -0.058173603105896
                                      -0.075000000000000  -0.075528612256965  -0.076507811395346  -0.077053736788907  -0.077045280667765
                                      -0.075000000000000  -0.075001925996158  -0.075003162498096  -0.075008925944951  -0.075013185036205];
expSolCPHistPenalty{4}(:,:,2,2) = [   -0.075000000000000  -0.076959485964026  -0.058070115158454  -0.020896658025894  -0.000003126200209
                                      -0.075000000000000  -0.077043756187787  -0.058153744469378  -0.020772875831456  -0.000015578987124
                                      -0.055094310254260  -0.056209937206396  -0.042426505357940  -0.015210225610549  -0.000000003354142
                                      -0.019590290622281  -0.019567334279353  -0.014438306252600  -0.005184134539053  -0.000011261196043
                                                       0   0.000004220615809   0.000024328124268   0.000020187769168   0.000005175922995];
expSolCPHistPenalty{4}(:,:,3,2) = [0   0.020854140599660   0.059243512950132   0.081032320528022   0.081064704910330
                                   0   0.020588008191575   0.059270843536643   0.080853622143460   0.080987530588905
                                   0   0.015111582424627   0.043313067195092   0.059289762889765   0.059309365608360
                                   0   0.005381251664295   0.015118596956551   0.020693587579539   0.020765911113862
                                   0   0.000009097133198   0.000029778337673   0.000041136074713   0.000034018776760];
expSolCPHistPenalty{4}(:,:,4,2) = [1.000000000000000   0.902368927062183   0.837281545103638   0.902368927062183   1.000000000000000
                                   0.902368927062183   0.814269680527354   0.755536849504136   0.814269680527354   0.902368927062183
                                   0.837281545103638   0.755536849504136   0.701040385771135   0.755536849504136   0.837281545103638
                                   0.902368927062183   0.814269680527354   0.755536849504136   0.814269680527354   0.902368927062183
                                   1.000000000000000   0.902368927062183   0.837281545103638   0.902368927062183   1.000000000000000];
                          
% Define the expected solution in terms of the residual history
expSolResHistoryPenaltyNLinear = 1e2*[ 3.948973554405663
                                       1.241784204262515
                                       0.515881258734050
                                       0.006559936144054
                                       0.000000909238555
                                       0.000000000105359
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0
                                                       0];
                               
% Define the expected solution in terms of the minimum element area size
epxSolMinElASize = [0.001293936786475
                    0.000074036937538
                    0.000148876380671
                    0.000367163557641];
                       
%% 8. Define the expected solutions for the Lagrange Multipliers method

% Defined the expected solution in terms of the Control Point Displacements
expSolDispLM = 1.0e+05 * [     0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
              -0.000000000321466
              -0.000000038103840
               0.000000020597901
               0.000000015979612
              -0.000000030958758
               0.000000018042036
               0.000000015374517
              -0.000000002860006
               0.000000006371317
               0.000000000015581
               0.000000000939444
              -0.000000000098036
              -0.000000000140109
              -0.000000022350228
               0.000000064037675
               0.000000018743696
              -0.000000019957523
               0.000000057110155
               0.000000037998114
              -0.000000005301393
               0.000000021736129
               0.000000000077546
               0.000000001801045
              -0.000000001690704
              -0.000000000116315
              -0.000000000316848
               0.000000069920614
               0.000000022901937
               0.000000001523025
               0.000000063261039
               0.000000039540449
              -0.000000001062107
               0.000000022974750
              -0.000000000067797
               0.000000000988893
              -0.000000002024875
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
              -0.000000000018657
               0.000000014879325
               0.000000006054174
               0.000000001774867
               0.000000012901381
               0.000000005899535
               0.000000004921978
               0.000000009888655
               0.000000005174359
               0.000000006202519
               0.000000006016772
               0.000000004092102
               0.000000004950236
               0.000000002179427
               0.000000002577404
               0.000000002030017
              -0.000000000023655
               0.000000000908510
               0.000000000004491
              -0.000000000031228
               0.000000000010146
               0.000000000112552
               0.000000030216716
               0.000000024870861
               0.000000004921343
               0.000000029186439
               0.000000023525107
               0.000000013048174
               0.000000023902829
               0.000000020462840
               0.000000017208561
               0.000000015229988
               0.000000015213848
               0.000000014454384
               0.000000005912480
               0.000000008755202
               0.000000006116123
              -0.000000000021066
               0.000000002686013
               0.000000000007955
              -0.000000000097509
               0.000000000045786
              -0.000000000101419
               0.000000032506920
               0.000000044095456
               0.000000006529876
               0.000000030486683
               0.000000042457123
               0.000000017690387
               0.000000025666231
               0.000000037040273
               0.000000024642530
               0.000000017394429
               0.000000027690184
               0.000000022305180
               0.000000007416695
               0.000000015488424
               0.000000009904720
               0.000000000068282
               0.000000004311900
              -0.000000000003915
              -0.000000000180101
               0.000000000077969
               0.000000000191790
               0.000000023271663
               0.000000058813681
               0.000000007119076
               0.000000021493648
               0.000000056362161
               0.000000019849150
               0.000000018024562
               0.000000049243115
               0.000000028753667
               0.000000012740478
               0.000000037508157
               0.000000027659274
               0.000000006037380
               0.000000021124232
               0.000000012752775
               0.000000000146533
               0.000000005597806
              -0.000000000028708
              -0.000000000251508
               0.000000000097205
              -0.000000000284502
               0.000000008525305
               0.000000068046730
               0.000000007459249
               0.000000007949263
               0.000000064697815
               0.000000021250506
               0.000000006636492
               0.000000056256082
               0.000000031298013
               0.000000004559656
               0.000000043199039
               0.000000030791818
               0.000000002771180
               0.000000024345045
               0.000000014249220
               0.000000000114022
               0.000000006447749
              -0.000000000059131
              -0.000000000293509
               0.000000000057107
               0.000000000103613
               0.000000000146610
               0.000000069910756
               0.000000008028247
               0.000000000699455
               0.000000067213328
               0.000000022300690
               0.000000000245619
               0.000000058018721
               0.000000032657255
              -0.000000000222448
               0.000000044409497
               0.000000031604911
               0.000000000882957
               0.000000024568135
               0.000000014318877
               0.000000000011765
               0.000000006702833
              -0.000000000064453
              -0.000000000284854
               0.000000000061798
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
              -0.000000000004687
               0.000000018780368
               0.000000007807620
              -0.000000002940164
               0.000000016107904
               0.000000007501893
              -0.000000007412148
               0.000000011024009
               0.000000006138291
              -0.000000007372084
               0.000000004346088
               0.000000004084399
              -0.000000003264183
              -0.000000000033095
               0.000000001432175
              -0.000000000007446
               0.000000000017935
              -0.000000000014166
               0.000000000097964
               0.000000033402489
               0.000000032756834
              -0.000000007295927
               0.000000032044058
               0.000000030753080
              -0.000000018553497
               0.000000023818419
               0.000000024852857
              -0.000000020427259
               0.000000010538931
               0.000000014600697
              -0.000000009756465
               0.000000000060095
               0.000000004210510
              -0.000000000004418
               0.000000000003728
              -0.000000000029293
              -0.000000000086452
               0.000000028620208
               0.000000054195316
              -0.000000008792045
               0.000000025577848
               0.000000051758914
              -0.000000023533397
               0.000000020054519
               0.000000042322424
              -0.000000028875797
               0.000000010109225
               0.000000025354408
              -0.000000015238279
               0.000000000222946
               0.000000006526798
               0.000000000037509
              -0.000000000047390
              -0.000000000089517
               0.000000000131329
               0.000000010407375
               0.000000067533317
              -0.000000009466218
               0.000000009758829
               0.000000063107968
              -0.000000026046578
               0.000000007626922
               0.000000052378514
              -0.000000033246287
               0.000000004236585
               0.000000031861537
              -0.000000018289111
               0.000000000099177
               0.000000007929054
               0.000000000097669
              -0.000000000034917
              -0.000000000231311
              -0.000000000209659
               0.000000000415305
               0.000000069911014
              -0.000000009891916
               0.000000000222078
               0.000000066555072
              -0.000000028365549
               0.000000000221987
               0.000000054818788
              -0.000000034481884
               0.000000000255881
               0.000000032660425
              -0.000000018462552
               0.000000000025296
               0.000000008070641
               0.000000000123643
              -0.000000000009526
              -0.000000000312101
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
              -0.000000000176724
              -0.000000024870839
               0.000000012891077
              -0.000000006174096
              -0.000000021178286
               0.000000010009405
              -0.000000011552525
              -0.000000010678955
               0.000000007366198
              -0.000000006112059
               0.000000000067553
               0.000000002611843
              -0.000000000013629
               0.000000000002485
               0.000000000012663
              -0.000000000347722
              -0.000000033803876
               0.000000043283293
              -0.000000011820426
              -0.000000031135660
               0.000000042302907
              -0.000000027245587
              -0.000000018847101
               0.000000028453222
              -0.000000017724538
              -0.000000000747699
               0.000000007609495
               0.000000000018428
               0.000000000204000
               0.000000000003015
               0.000000000075510
              -0.000000014367086
               0.000000066328201
              -0.000000012269559
              -0.000000012109063
               0.000000060763591
              -0.000000032844790
              -0.000000008129751
               0.000000043880570
              -0.000000024947814
              -0.000000000606329
               0.000000011213651
               0.000000000102797
               0.000000000174664
              -0.000000000152564
              -0.000000000307782
               0.000000000473444
               0.000000069913969
              -0.000000013546534
               0.000000000076119
               0.000000065572120
              -0.000000035859015
               0.000000000362172
               0.000000046117689
              -0.000000025196604
               0.000000000015886
               0.000000011030201
               0.000000000121325
               0.000000000017173
              -0.000000000059114
               0.051011922222595
              -0.003535736228393
               0.006370673117514
               0.003775470871255
              -0.041147589216426
               0.049573950374596
              -0.060440597470262
              -0.047641313498911
               0.050411039441831
              -0.000006215994394
               0.000133822174801
               0.000022597765540
               0.004529415953464
               0.050997045826744
              -0.010688959884007
               6.938578688394112
               0.044057732296866
               4.451861905651795
               4.838393313832296
               0.040711435064606
               7.572994587887181
               0.862877745897475
               0.024503432703226
               8.462959361529247
              -4.340808845958854
              -0.002046643865702
               6.922914546249020
              -9.626071385023932
               0.000008601619883
               0.004003475153182
              -0.000001641297187
              -0.350995071092061
               0.000002877932485
              -0.711965978743103
              -0.029156185818164
              -0.000722301847589
               0.001190756451748
              -0.022760592741357
               0.032244461980144
               0.027053763509412
               0.007317342629040
               0.043709395452941
               0.060848986747436
               0.062990988910191
               0.048723785355937
               0.053833414140360
               0.000000168142454
              -0.000069230171316
               0.000000898027678
               0.001751724100311
              -0.000018910450366
               0.004590218478746
               0.047999468282224
              -0.052756679840518
              -0.216682354668623
               0.040105509080368
               0.180956770593666
              -0.074532773698849
               0.004839263331163
               0.222212302413601
               0.218528134364214
               0.000024576819026
               0.004444776885759
              -0.000028297618435
              -0.015605972950425];
               
% Define the expected solution in terms of the deformation
expSolCPHistLM = struc([]);
for counterPatches = 1:4
    for counterComp = 1:4
        expSolCPHistLM{counterPatches}(:,:,counterComp,1) = ...
            BSplinePatches{counterPatches}.CP(:,:,counterComp);
    end
end

% Patch 1 :
% _________

expSolCPHistLM{1}(:,:,1,2) = [                 0  -0.000032146605525  -0.000014010938199  -0.000011631507319
                               0.031066017177982   0.032663978353931   0.032940386820309   0.033356210855395
                               0.075000000000000   0.076537451666027   0.078799811358669   0.078954044942631
                               0.075000000000000   0.075001558132746   0.075007754570372   0.074993220345232];
expSolCPHistLM{1}(:,:,2,2) = [-0.075000000000000  -0.078810383971464  -0.033301039997237  -0.000031684806260
                              -0.075000000000000  -0.078095875803864  -0.033061769460850   0.000152302452437
                              -0.031066017177982  -0.031352017801011  -0.013398104970748  -0.000106210725497
                                               0   0.000093944409887   0.000180104462573   0.000098889315314];
expSolCPHistLM{1}(:,:,3,2) = [ 0   0.033125807326082   0.081403767524589   0.081992061400157
                               0   0.032870220788683   0.080711015526880   0.081326103943167
                               0   0.013505097374888   0.033239630126625   0.033363492160224
                               0  -0.000009803573115  -0.000169070424814  -0.000202487509165];
expSolCPHistLM{1}(:,:,4,2) = [ 1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000
                               0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                               0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                               1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000];

% Patch 2 :
% _________

expSolCPHistLM{2}(:,:,1,2) = [                 0  -0.000001865710677   0.000011255228883  -0.000010141857720   0.000019179035653  -0.000028450194990   0.000010361280089
                               0.011266582861675   0.011444069591987   0.011758717121990   0.011919570481533   0.011978490477185   0.012012507732112   0.012069407545433
                               0.033674222389275   0.034166420221619   0.034979039746175   0.035443261130960   0.035659137421297   0.035799273032626   0.035904291364617
                               0.053765960839350   0.054386212777819   0.055486816898210   0.056230213866669   0.056641327545363   0.056895762111844   0.057031686356730
                               0.068112370398213   0.068607393981771   0.069557808817865   0.070342888350707   0.070878297797435   0.071191552157278   0.071272861490485
                               0.075000000000000   0.075203001710896   0.075611612338653   0.075990471953344   0.076275277539184   0.076424921978346   0.076431887729676
                               0.075000000000000   0.075000449105566   0.075000795524435   0.074999608535419   0.074997129156304   0.074994086935521   0.074993554660348];
expSolCPHistLM{2}(:,:,2,2) = [ 0.075000000000000   0.076487932527495   0.071134042034297   0.057016652872098   0.036001388668354   0.012119113392688   0.000014660972843
                               0.075000000000000   0.076290138148014   0.071031014269659   0.056814629121337   0.035823587206288   0.012061509121229   0.000069945535415
                               0.068112370398213   0.069101235857193   0.064247549591881   0.051394983618221   0.032384204284844   0.010895564751503   0.000024561873300
                               0.053765960839350   0.054367638013372   0.050351359289972   0.040283156878243   0.025414406747353   0.008532747603514  -0.000022244750841
                               0.033674222389275   0.033892165098419   0.031172996087091   0.024882028463211   0.015723114704633   0.005335696934305   0.000088295730290
                               0.011266582861675   0.011264217349930   0.010229808941800   0.008083610212222   0.005073232180442   0.001703880750745   0.000001176515842
                                               0  -0.000003122811791  -0.000009750906185  -0.000018010054567  -0.000025150847901  -0.000029350933826  -0.000028485353941];
expSolCPHistLM{2}(:,:,3,2) = [ 0   0.011872000309909   0.036161308491972   0.058175506432955   0.073993738533050   0.081804672960174   0.081991075639692
                               0   0.011856536394204   0.036026733054183   0.058011673149638   0.073748586454951   0.081469781542597   0.081721332780198
                               0   0.010749351438943   0.032628032061107   0.052532387858364   0.066781578215512   0.073737978576046   0.073914242515717
                               0   0.008485992252467   0.025661743759739   0.041312732356521   0.052579176263193   0.058085864743035   0.058206910552468
                               0   0.005316319315323   0.015994896941515   0.025689201411219   0.032694171279654   0.036108726844178   0.036131035908115
                               0   0.001783329502245   0.005327180219116   0.008507972079525   0.010791696177757   0.011911357745796   0.011936866170277
                               0   0.000001014618711   0.000004578624563   0.000007796935192   0.000009720454134   0.000005710683075   0.000006179816411];
expSolCPHistLM{2}(:,:,4,2) = [ 1.000000000000000   0.941421356237309   0.871126983722081   0.847695526217005   0.871126983722081   0.941421356237309   1.000000000000000
                               0.941421356237309   0.886274169979695   0.820097546470558   0.798038671967512   0.820097546470558   0.886274169979695   0.941421356237309
                               0.871126983722081   0.820097546470558   0.758862221768731   0.738450446868122   0.758862221768731   0.820097546470558   0.871126983722081
                               0.847695526217005   0.798038671967513   0.738450446868122   0.718587705168325   0.738450446868122   0.798038671967513   0.847695526217005
                               0.871126983722081   0.820097546470558   0.758862221768731   0.738450446868122   0.758862221768731   0.820097546470558   0.871126983722081
                               0.941421356237309   0.886274169979695   0.820097546470558   0.798038671967512   0.820097546470558   0.886274169979695   0.941421356237309
                               1.000000000000000   0.941421356237309   0.871126983722081   0.847695526217005   0.871126983722081   0.941421356237309   1.000000000000000];

% Patch 3 :
% _________

expSolCPHistLM{3}(:,:,1,2) = [                 0  -0.000000468702408   0.000009796374874  -0.000008645152161   0.000013132919013  -0.000020965900109
                              -0.014305767737291  -0.014599784177639  -0.015035360482529  -0.015184972253730  -0.015252389509049  -0.015294959344280
                              -0.042049512883487  -0.042790727730793  -0.043904862631903  -0.044402852621003  -0.044654170714259  -0.044886067737164
                              -0.064016504294496  -0.064753712681949  -0.066059230159016  -0.066904084002481  -0.067341133025780  -0.067464692667074
                              -0.075000000000000  -0.075326418250513  -0.075975646488869  -0.076523827926335  -0.076828911092694  -0.076846255166899
                              -0.075000000000000  -0.075000744600169  -0.075000441793360  -0.074996249112920  -0.074990233148135  -0.074987635709384];
expSolCPHistLM{3}(:,:,2,2) = [ 0.075000000000000   0.076878036842928   0.067356753147098   0.044911533697105   0.015346505260477   0.000041530501855
                               0.075000000000000   0.076610790449901   0.067220910125149   0.044607297716803   0.015281650606112   0.000022207767513
                               0.064016504294496   0.065118905177908   0.057023346240341   0.037896956158542   0.012973428758969   0.000022198695328
                               0.042049512883487   0.042484121704821   0.036945397382318   0.024586409581583   0.008444332649788   0.000025588059148
                               0.014305767737291   0.014302458197361   0.012216746092519   0.008042968775498   0.002738650861859   0.000002529604324
                                               0   0.000001793501582   0.000000372768458  -0.000004739023268  -0.000003491729636  -0.000000952643808];
expSolCPHistLM{3}(:,:,3,2) = [ 0   0.015086529751723   0.045325196276806   0.069436035860372   0.081753331664855   0.081991101434340
                               0   0.015055957028123   0.045124820853078   0.069192395711932   0.081310796829235   0.081655507228992
                               0   0.012824565658535   0.038376789995765   0.058873746711299   0.069254355732045   0.069498383047265
                               0   0.008429114058011   0.025035556776806   0.038426945067768   0.045235666592284   0.045315555336678
                               0   0.002871950744098   0.008441725241926   0.012863416325920   0.015098673157208   0.015112831844108
                               0  -0.000001416550498  -0.000002929266869  -0.000008951742206  -0.000023131098372  -0.000031210142254]; 
expSolCPHistLM{3}(:,:,4,2) = [ 1.000000000000000   0.926776695296637   0.853553390593274   0.853553390593274   0.926776695296637   1.000000000000000
                               0.926776695296637   0.858915042944955   0.791053390593274   0.791053390593274   0.858915042944955   0.926776695296637
                               0.853553390593274   0.791053390593274   0.728553390593274   0.728553390593274   0.791053390593274   0.853553390593274
                               0.853553390593274   0.791053390593274   0.728553390593274   0.728553390593274   0.791053390593274   0.853553390593274
                               0.926776695296637   0.858915042944955   0.791053390593274   0.791053390593274   0.858915042944955   0.926776695296637
                               1.000000000000000   0.926776695296637   0.853553390593274   0.853553390593274   0.926776695296637   1.000000000000000];

% Patch 4 :
% _________

expSolCPHistLM{4}(:,:,1,2) = [                 0  -0.000017672357994  -0.000034772204790   0.000007551007805  -0.000030778234052
                              -0.019590290622281  -0.020207700194503  -0.020772333184763  -0.020817246541445  -0.020944944068010
                              -0.055094310254260  -0.056249562766028  -0.057818868915565  -0.058378789208030  -0.058680211713852
                              -0.075000000000000  -0.075611205923038  -0.076772453756885  -0.077494781384902  -0.077519660406084
                              -0.075000000000000  -0.075001362938089  -0.074998157173618  -0.074989720323655  -0.074987867495569];
expSolCPHistLM{4}(:,:,2,2) = [-0.075000000000000  -0.077487083942449  -0.058474697834411  -0.021026999237985   0.000047344447225
                              -0.075000000000000  -0.077117828645025  -0.058207876253905  -0.020801196931501   0.000007611925565
                              -0.055094310254260  -0.056162205723327  -0.042356483696506  -0.015203822421602   0.000036217186428
                              -0.019590290622281  -0.019583535316476  -0.014465617242770  -0.005177692743359   0.000001588639542
                                               0   0.000000248507115   0.000020400015274   0.000017466376175   0.000001717306728];
expSolCPHistLM{4}(:,:,3,2) = [ 0   0.020879398310980   0.059422639508121   0.081632820090889   0.081991396907504
                               0   0.020591231170938   0.059324600934040   0.081076359064305   0.081557211951366
                               0   0.015127467077037   0.043317095829461   0.059482367261407   0.059706079124466
                               0   0.005378244167939   0.015151796825473   0.020711655737563   0.020693310761616
                               0   0.000001266308846   0.000000301531773  -0.000015256403378  -0.000005911359011];
expSolCPHistLM{4}(:,:,4,2) = [ 1.000000000000000   0.902368927062183   0.837281545103638   0.902368927062183   1.000000000000000
                               0.902368927062183   0.814269680527354   0.755536849504136   0.814269680527354   0.902368927062183
                               0.837281545103638   0.755536849504136   0.701040385771135   0.755536849504136   0.837281545103638
                               0.902368927062183   0.814269680527354   0.755536849504136   0.814269680527354   0.902368927062183
                               1.000000000000000   0.902368927062183   0.837281545103638   0.902368927062183   1.000000000000000];
    
% Define the expected solution in terms of the residual history
expSolResHistoryLMNLinear = 1e2*[  3.948973554405663
                                   1.154522711617360
                                   0.437364599395353
                                   0.002838880702911
                                   0.000000580619310
                                   0.000000000080675
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0
                                                   0];
                           
%% 9. Define parameters related to the nonlinear analysis
propNLinearAnalysis.method = 'newtonRapshon';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 1e-6;
propNLinearAnalysis.maxIter = 50;

%% 10. Solve the nonlinear coupled system using the penalty method
propCouplingPenalty.alphaD = zeros(connections.No,1);
propCouplingPenalty.alphaR = zeros(connections.No,1);
for iConnections = 1:connections.No
    propCouplingPenalty.alphaD(iConnections,1) = norm(Dm)*1e3;
    propCouplingPenalty.alphaR(iConnections,1) = norm(Db)*1e3;
end
propCouplingPenalty.intC = intC;
plot_IGANLinear = '';
[dHatPenaltyNLinear,CPHistoryPenaltyNLinear,resHistoryPenaltyNLinear,~] = ...
    solve_DDMPenaltyIGAKirchhoffLoveShellMultipatchesNLinear...
    (BSplinePatches,connections,propCouplingPenalty,propNLinearAnalysis,...
    solve_LinearSystem,plot_IGANLinear,[],'');

%% 11. Solve the nonlinear coupled system using the augmented Lagrange Multipliers method
propCouplingLM.alphaD = zeros(connections.No,1);
propCouplingLM.alphaR = zeros(connections.No,1);
for iConnections = 1:connections.No
    propCouplingLM.alphaD(iConnections,1) = norm(Dm)*1e1;
    propCouplingLM.alphaR(iConnections,1) = norm(Db)*1e1;
end
propCouplingLM.intC = intC;
plot_IGANLinear = '';
[dHatLMNLinear,CPHistoryLMNLinear,resHistoryLMNLinear,~] = ...
    solve_DDMLagrangeMultipliersIGAKLShellMultipatchesNLinear...
    (BSplinePatches,connections,propCouplingLM,propNLinearAnalysis,...
    solve_LinearSystem,plot_IGANLinear,[],'');

%% 12. Verify the solution
testCase.verifyEqual(dHatPenaltyNLinear,expSolDispPenalty,'AbsTol',absTol);
testCase.verifyEqual(CPHistoryPenaltyNLinear,expSolCPHistPenalty,'AbsTol',absTol);
testCase.verifyEqual(resHistoryPenaltyNLinear,expSolResHistoryPenaltyNLinear,'AbsTol',absTolRelaxed4);
testCase.verifyEqual([BSplinePatches{1}.minElArea;BSplinePatches{2}.minElArea;BSplinePatches{3}.minElArea;BSplinePatches{4}.minElArea],epxSolMinElASize,'AbsTol',absTol);
testCase.verifyEqual(dHatLMNLinear,expSolDispLM,'AbsTol',absTolRelaxed7);
testCase.verifyEqual(CPHistoryLMNLinear,expSolCPHistLM,'AbsTol',absTol);
testCase.verifyEqual(resHistoryLMNLinear,expSolResHistoryLMNLinear,'AbsTol',absTolRelaxed4);

end
