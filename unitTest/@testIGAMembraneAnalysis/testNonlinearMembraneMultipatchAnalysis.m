function testNonlinearMembraneMultipatchAnalysis(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Test the isogeometric membrane analysis over the semisphere modelled with
% 4 patches and solved using the penalty and the Lagrange Multipliers
% methods for the geometrically nonlinear case.
% 
% Function :
%
% 0. Read input
%
% 1. Define the multipatch geometry in terms of NURBS
%
% 2. Define the material parameters
%
% 3. GUI
%
% 4. Define the refinement
%
% 5. Define the boundary conditions
%
% 6. Create the patches and the Lagrange Multiplier fields
%
% 7. Define the expected solutions for the penalty method
%
% 8. Define the expected solutions for the Lagrange Multipliers method
%
% 9. Define the nonlinear analysis parameters
%
% 10. Solve the multipatch problem with the Penalty method
%
% 11. Solve the multipatch problem with the Lagrange Multipliers method
%
% 12. Verify the solution
%
%% Function main body

%% 0. Read input

% Define absolute tolerance
relaxed = 1e5;
absTol = 1e-15*relaxed;

%% 1. Define the multipatch geometry in terms of NURBS

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
pLambda12 = 0;

% Knot vector
XiLambda12 = [0 1];

% Control points weights
CPLambda12(:,4) = (1);

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSLambda12 = 0;
nxiLambda12 = length(CPLambda12(:,1,1));
for i=1:nxiLambda12
    if CPLambda12(i,4)~=1
        isNURBSLambda12 = 1;
        break;
    end
end

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

% Polynomial degree
pMu12 = 0;

% Knot vector
XiMu12 = [0 1];

% Control point weights
CPMu12(:,4) = (1);

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu12 = 0;
nxiMu12 = length(CPMu12(:,1,1));
for i=1:nxiMu12
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
pLambda23 = 0;

% Knot vector
XiLambda23 = [0 1];

% Control points weights
CPLambda23(:,4) = (1);

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSLambda23 = 0;
nxiLambda23 = length(CPLambda23(:,1,1));
for i=1:nxiLambda23
    if CPLambda23(i,4)~=1
        isNURBSLambda23 = 1;
        break;
    end
end

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

% Polynomial degree
pMu23 = 0;

% Knot vector
XiMu23 = [0 1];

% Control point weights
CPMu23(:,4) = (1);

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu23 = 0;
nxiMu23 = length(CPMu23(:,1,1));
for i=1:nxiMu23
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
pLambda34 = 0;

% Knot vector
XiLambda34 = [0 1];

% Control points weights
CPLambda34(:,4) = (1);

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSLambda34 = 0;
nxiLambda34 = length(CPLambda34(:,1,1));
for i=1:nxiLambda34
    if CPLambda34(i,4)~=1
        isNURBSLambda34 = 1;
        break;
    end
end

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

% Polynomial degree
pMu34 = 0;

% Knot vector
XiMu34 = [0 1];

% Control point weights
CPMu34(:,4) = (1);

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu34 = 0;
nxiMu34 = length(CPMu34(:,1,1));
for i=1:nxiMu34
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
pLambda14 = 0;

% Knot vector
XiLambda14 = [0 1];

% Control points weights
CPLambda14(:,4) = (1);

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSLambda14 = 0;
nxiLambda14 = length(CPLambda14(:,1,1));
for i=1:nxiLambda14
    if CPLambda14(i,4)~=1
        isNURBSLambda14 = 1;
        break;
    end
end

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

% Polynomial degree
pMu14 = 0;

% Knot vector
XiMu14 = [0 1];

% Control point weights
CPMu14(:,4) = (1);

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu14 = 0;
nxiMu14 = length(CPMu14(:,1,1));
for i=1:nxiMu14
    if CPMu14(i,4)~=1
        isNURBSMu14 = 1;
        break;
    end
end

%% 2. Define the material parameters

% general parameters
EYoung = 21e6;
nue = 0.0;
thickness = .03;
prestress.voigtVector = [10/thickness
                         10/thickness
                         0];
density = 7860;

% Patch 1 :
% _________

% Young's modulus
parameters1.E = EYoung;

% Poisson ratio
parameters1.nue = nue;

% Thickness of the membrane
parameters1.t = thickness;

% Density of the membrane
parameters1.rho = density;

% Prestress of the membrane
parameters1.prestress = prestress;

% Patch 2 :
% _________

% Young's modulus
parameters2.E = EYoung;

% Poisson ratio
parameters2.nue = nue;

% Thickness of the plate
parameters2.t = thickness;

% Density of the membrane
parameters2.rho = density;

% Prestress of the membrane
parameters2.prestress = prestress;

% Patch 3 :
% _________

% Young's modulus
parameters3.E = EYoung;

% Poisson ratio
parameters3.nue = nue;

% Thickness of the plate
parameters3.t = .03;

% Density of the membrane
parameters3.rho = density;

% Prestress of the membrane
parameters3.prestress = prestress;

% Patch 4 :
% _________

% Young's modulus
parameters4.E = EYoung;

% Poisson ratio
parameters4.nue = nue;

% Thickness of the plate
parameters4.t = thickness;

% Density of the membrane
parameters4.rho = density;

% Prestress of the membrane
parameters4.prestress = prestress;

%% 3. GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points

% Choose linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Patch 1 :
% _________

int1.type = 'user';
if strcmp(int1.type,'user')
    int1.xiNGP = p1 + 2;
    int1.etaNGP = q1 + 2;
    int1.xiNGPForLoad = ceil((p1 + 2)/2);
    int1.etaNGPForLoad = ceil((q1 + 2)/2);
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
    int3.xiNGP = p3 + 2;
    int3.etaNGP = q3 + 2;
    int3.xiNGPForLoad = ceil((p3 + 2)/2);
    int3.etaNGPForLoad = ceil((q3 + 2)/2);
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
        intC.nGP1 = 12;
        intC.nGP2 = 12;
    else
        intC.nGP = 12;
    end
    intC.nGPError = 12;
end

% On the coupling
% _______________

% Compute the material matrices for the membrane and the bending part
Dm = EYoung*thickness/(1-nue^2)*...
      [1   nue  0
       nue 1    0
       0   0    (1-nue)/2];
Db = thickness^3/12*Dm;

%% 4. Define the refinement

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

tp1 = 1;
tq1 = 1;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% Patch 2 :
% _________

tp2 = 0;
tq2 = 0;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');

% Patch 3 :
% _________

tp3 = 1;
tq3 = 1;
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

tpLambda12 = 0;
[XiLambda12,CPLambda12,pLambda12] = degreeElevateBSplineCurve...
    (pLambda12,XiLambda12,CPLambda12,tpLambda12,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

tpMu12 = 0;
[XiMu12,CPMu12,pMu12] = degreeElevateBSplineCurve...
    (pMu12,XiMu12,CPMu12,tpMu12,'');

%%% Coupling between patch 2 and patch 3 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

tpLambda23 = 0;
[XiLambda23,CPLambda23,pLambda23] = degreeElevateBSplineCurve...
    (pLambda23,XiLambda23,CPLambda23,tpLambda23,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

tpMu23 = 0;
[XiMu23,CPMu23,pMu23] = degreeElevateBSplineCurve...
    (pMu23,XiMu23,CPMu23,tpMu23,'');

%%% Coupling between patch 3 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

tpLambda34 = 0;
[XiLambda34,CPLambda34,pLambda34] = degreeElevateBSplineCurve...
    (pLambda34,XiLambda34,CPLambda34,tpLambda34,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

tpMu34 = 0;
[XiMu34,CPMu34,pMu34] = degreeElevateBSplineCurve...
    (pMu34,XiMu34,CPMu34,tpMu34,'');

%%% Coupling between patch 1 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

tpLambda14 = 0;
[XiLambda14,CPLambda14,pLambda14] = degreeElevateBSplineCurve...
    (pLambda14,XiLambda14,CPLambda14,tpLambda14,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

tpMu14 = 0;
[XiMu14,CPMu14,pMu14] = degreeElevateBSplineCurve...
    (pMu14,XiMu14,CPMu14,tpMu14,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

noKnotsXi1 = 1;
noKnotsEta1 = 1;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface(p1,Xi1,q1,Eta1,CP1,noKnotsXi1,noKnotsEta1,'');

% Patch 2 :
% _________

noKnotsXi2 = 2;
noKnotsEta2 = 2;
[Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface(p2,Xi2,q2,Eta2,CP2,noKnotsXi2,noKnotsEta2,'');

% Patch 3 :
% _________

noKnotsXi3 = 1;
noKnotsEta3 = 1;
[Xi3,Eta3,CP3] = knotRefineUniformlyBSplineSurface(p3,Xi3,q3,Eta3,CP3,noKnotsXi3,noKnotsEta3,'');

% Patch 4 :
% _________

noKnotsXi4 = 2;
noKnotsEta4 = 2;
[Xi4,Eta4,CP4] = knotRefineUniformlyBSplineSurface(p4,Xi4,q4,Eta4,CP4,noKnotsXi4,noKnotsEta4,'');

% Scale down the discretization of the Lagrange Multipliers field
scale = 1;

%%% Coupling between patch 1 and patch 2 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

% nLambda12 = min(nKnotsEta1,nKnotsEta2);
noLambda12 = ceil(min(noKnotsEta1,noKnotsEta2)/scale);
[XiLambda12,CPLambda12] = knotRefineUniformlyBSplineCurve...
    (noLambda12,pLambda12,XiLambda12,CPLambda12,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

noMu12 = ceil(noLambda12/2);
[XiMu12,CPMu12] = knotRefineUniformlyBSplineCurve...
    (noMu12,pMu12,XiMu12,CPMu12,'');

%%% Coupling between patch 2 and patch 3 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

% nLambda23 = min(nKnotsEta2,nKnotsEta3);
noLambda23 = ceil(min(noKnotsEta2,noKnotsEta3)/scale);
[XiLambda23,CPLambda23] = knotRefineUniformlyBSplineCurve...
    (noLambda23,pLambda23,XiLambda23,CPLambda23,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

noMu23 = ceil(noLambda23/2);
[XiMu23,CPMu23] = knotRefineUniformlyBSplineCurve...
    (noMu23,pMu23,XiMu23,CPMu23,'');

%%% Coupling between patch 3 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

% nLambda34 = min(nKnotsEta3,nKnotsEta4);
noLambda34 = ceil(min(noKnotsEta3,noKnotsEta4)/scale);
[XiLambda34,CPLambda34] = knotRefineUniformlyBSplineCurve...
    (noLambda34,pLambda34,XiLambda34,CPLambda34,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

noMu34 = ceil(noLambda34/2);
[XiMu34,CPMu34] = knotRefineUniformlyBSplineCurve...
    (noMu34,pMu34,XiMu34,CPMu34,'');

%%% Coupling between patch 1 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

% nLambda14 = min(nKnotsEta1,nKnotsEta4);
noLambda14 = ceil(min(noKnotsEta1,noKnotsEta4)/scale);
[XiLambda14,CPLambda14] = knotRefineUniformlyBSplineCurve...
    (noLambda14,pLambda14,XiLambda14,CPLambda14,'');

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

noMu14 = ceil(noLambda14/2);
[XiMu14,CPMu14] = knotRefineUniformlyBSplineCurve...
    (noMu14,pMu14,XiMu14,CPMu14,'');

%% 5. Define the boundary conditions

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
xisup1 = [1 1];   etasup1 = [0 1];
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
cables1.No = 0;

% Patch 2 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs2 = [];
xisup2 = [0 1];   etasup2 = [0 0];
for dir = 1:3
    homDOFs2 = findDofs3D...
        (homDOFs2,xisup2,etasup2,dir,CP2);
end
xisup2 = [1 1];   etasup2 = [0 1];
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
cables2.No = 0;

% Patch 3 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs3 = [];
xisup3 = [0 1];   etasup3 = [0 0];
for dir = 1:3
    homDOFs3 = findDofs3D...
        (homDOFs3,xisup3,etasup3,dir,CP3);
end
xisup3 = [1 1];   etasup3 = [0 1];
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
cables3.No = 0;

% Patch 4 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs4 = [];
xisup4 = [0 1];   etasup4 = [0 0];
for dir = 1:3
    homDOFs4 = findDofs3D...
        (homDOFs4,xisup4,etasup4,dir,CP4);
end
xisup4 = [1 1];   etasup4 = [0 1];
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
cables4.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parameter
loadAmplitude = 1e1;

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

%%% Coupling between patch 1 and patch 2 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

lambda12 = fillUpLagrangeMultipliers...
    (pLambda12,XiLambda12,CPLambda12,isNURBSLambda12);

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

mu12 = fillUpLagrangeMultipliers...
    (pMu12,XiMu12,CPMu12,isNURBSMu12);

%%% Coupling between patch 2 and patch 3 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

lambda23 = fillUpLagrangeMultipliers...
    (pLambda23,XiLambda23,CPLambda23,isNURBSLambda23);

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

mu23 = fillUpLagrangeMultipliers...
    (pMu23,XiMu23,CPMu23,isNURBSMu23);

%%% Coupling between patch 3 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

lambda34 = fillUpLagrangeMultipliers...
    (pLambda34,XiLambda34,CPLambda34,isNURBSLambda34);

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

mu34 = fillUpLagrangeMultipliers...
    (pMu34,XiMu34,CPMu34,isNURBSMu34);

%%% Coupling between patch 1 and patch 4 %%%
%%% ------------------------------------ %%%

% Lagrange multipliers field for the traction forces :
% ____________________________________________________

lambda14 = fillUpLagrangeMultipliers...
    (pLambda14,XiLambda14,CPLambda14,isNURBSLambda14);

% Lagrange multipliers field for the traction moments :
% _____________________________________________________

mu14 = fillUpLagrangeMultipliers...
    (pMu14,XiMu14,CPMu14,isNURBSMu14);

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch1 patch2 patch3 patch4};

% Collect all the Lagrange Multipliers into the connections :
% __________________________________________________________

connections.lambda = {lambda12 lambda23 lambda34 lambda14};
connections.mu = {mu12 mu23 mu34 mu14};

%% 7. Define the expected solutions for the penalty method

% Define the expected solution for the Penalty method in terms of the
% Control Point displacements
expSolDispPenalty = 1e-5*[             0
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
                      -0.000662826391263
                       0.104426857965901
                      -0.080638560672446
                      -0.067040880612507
                       0.247715622750825
                      -0.075206507361432
                      -0.389044215642538
                       0.084171045512374
                      -0.070025089655701
                                       0
                                       0
                                       0
                      -0.001708414895743
                       0.084603085324231
                      -0.141274127631648
                      -0.073435541460009
                       0.039275690869899
                      -0.186616512669172
                      -0.078434222395972
                       0.014766253697020
                      -0.069848558933217
                                       0
                                       0
                                       0
                      -0.000095201909076
                       0.000148706929923
                      -0.140218564557358
                      -0.075111393104354
                       0.001480026458863
                      -0.132443564645838
                      -0.090652576933179
                       0.000978527665412
                      -0.069027317852856
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
                       0.000653005888039
                      -0.071966389240653
                      -0.057267603225024
                      -0.046292782355552
                      -0.124018607378110
                      -0.061513674927784
                      -0.159375743571629
                      -0.045611131102642
                      -0.033530188490378
                                       0
                                       0
                                       0
                       0.001065368123393
                      -0.061829458236549
                      -0.140700105369012
                      -0.051506160047071
                      -0.044018151767629
                      -0.143203842973418
                      -0.072015274201766
                      -0.010033498035540
                      -0.051809897989086
                                       0
                                       0
                                       0
                       0.000174624448081
                       0.000194832004259
                      -0.140286396195981
                      -0.054679043161273
                       0.001005450358074
                      -0.134210672523780
                      -0.062629819577832
                       0.000803058612658
                      -0.049303629549547
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
                       0.000662826391257
                      -0.104426857965875
                      -0.080638560672434
                       0.067040880612463
                      -0.247715622750863
                      -0.075206507361514
                       0.389044215642645
                      -0.084171045512383
                      -0.070025089655678
                                       0
                                       0
                                       0
                       0.001708414895718
                      -0.084603085324230
                      -0.141274127631615
                       0.073435541460033
                      -0.039275690869791
                      -0.186616512669154
                       0.078434222395781
                      -0.014766253697014
                      -0.069848558933159
                                       0
                                       0
                                       0
                       0.000095201909016
                      -0.000148706929902
                      -0.140218564557359
                       0.075111393104302
                      -0.001480026458797
                      -0.132443564645761
                       0.090652576933161
                      -0.000978527665414
                      -0.069027317852803
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
                      -0.000653005888042
                       0.071966389240669
                      -0.057267603225035
                       0.046292782355568
                       0.124018607378163
                      -0.061513674927828
                       0.159375743571583
                       0.045611131102626
                      -0.033530188490363
                                       0
                                       0
                                       0
                      -0.001065368123429
                       0.061829458236558
                      -0.140700105369034
                       0.051506160046966
                       0.044018151767783
                      -0.143203842973316
                       0.072015274201827
                       0.010033498035489
                      -0.051809897989103
                                       0
                                       0
                                       0
                      -0.000174624448140
                      -0.000194832004239
                      -0.140286396195981
                       0.054679043161220
                      -0.001005450358027
                      -0.134210672523728
                       0.062629819577818
                      -0.000803058612654
                      -0.049303629549507
                                       0
                                       0
                                       0];

% Define the expected solution in terms of the deformation of the Control
% Points
expSolCPHistPenalty = struc([]);
for counterPatches = 1:4
    for counterComp = 1:4
        expSolCPHistPenalty{counterPatches}(:,:,counterComp,1) = ...
            BSplinePatches{counterPatches}.CP(:,:,counterComp);
    end
end

% Patch 1 :
% _________

expSolCPHistPenalty{1}(:,:,1,2) = [                0  -0.000000006628264  -0.000000017084149  -0.000000000952019
                                   0.043933982822018   0.043933312413212   0.043933248466603   0.043933231708087
                                   0.075000000000000   0.074996109557844   0.074999215657776   0.074999093474231
                                   0.075000000000000   0.075000000000000   0.075000000000000   0.075000000000000];
expSolCPHistPenalty{1}(:,:,2,2) = [-0.075000000000000  -0.074998955731420  -0.043933136791165   0.000000001487069
                                   -0.075000000000000  -0.074997522843772  -0.043933590065109   0.000000014800265
                                   -0.043933982822018  -0.043933141111563  -0.025735783625535   0.000000009785277
                                                       0                   0                   0                   0];
expSolCPHistPenalty{1}(:,:,3,2) = [0   0.043933176436411   0.074998587258724   0.074998597814354
                                   0   0.043933230756944   0.074998133834873   0.074998675564354
                                   0   0.025735231037175   0.043933284336429   0.043933292548839
                                   0                   0                   0                   0];
expSolCPHistPenalty{1}(:,:,4,2) = [1.000000000000000   0.804737854124365   0.804737854124365   1.000000000000000
                                   0.804737854124365   0.647603013860688   0.647603013860688   0.804737854124365
                                   0.804737854124365   0.647603013860688   0.647603013860688   0.804737854124365
                                   1.000000000000000   0.804737854124365   0.804737854124365   1.000000000000000];

% Patch 2 :
% _________

expSolCPHistPenalty{2}(:,:,1,2) = [0   0.000000006530059   0.000000010653681   0.000000001746244
                                   0.031066017177982   0.031065554250159   0.031065502116382   0.031065470387551
                                   0.075000000000000   0.074998406242564   0.074999279847258   0.074999373701804
                                   0.075000000000000   0.075000000000000   0.075000000000000   0.075000000000000];
expSolCPHistPenalty{2}(:,:,2,2) = [0.075000000000000   0.074999280336108   0.031065398883400   0.000000001948320
                                   0.075000000000000   0.074998759813926   0.031065576996464   0.000000010054504
                                   0.031066017177982   0.031065561066671   0.012867865309055   0.000000008030586
                                                   0                   0                   0                   0];
expSolCPHistPenalty{2}(:,:,3,2) = [0   0.031065444501950   0.074998592998946   0.074998597136038
                                   0   0.031065402041233   0.074998567961570   0.074998657893275
                                   0   0.012867630342151   0.031065499079002   0.031065524141687
                                   0                   0                   0                   0];
expSolCPHistPenalty{2}(:,:,4,2) = [1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000
                                   0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                                   0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                                   1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000];

% Patch 3 :
% _________

expSolCPHistPenalty{3}(:,:,1,2) = [                    0   0.000000006628264   0.000000017084149   0.000000000952019
                                      -0.043933982822018  -0.043933312413212  -0.043933248466603  -0.043933231708087
                                      -0.075000000000000  -0.074996109557844  -0.074999215657776  -0.074999093474231
                                      -0.075000000000000  -0.075000000000000  -0.075000000000000  -0.075000000000000];
expSolCPHistPenalty{3}(:,:,2,2) = [0.075000000000000   0.074998955731420   0.043933136791165  -0.000000001487069
                                   0.075000000000000   0.074997522843772   0.043933590065109  -0.000000014800265
                                   0.043933982822018   0.043933141111563   0.025735783625535  -0.000000009785277
                                                   0                   0                   0                   0];
expSolCPHistPenalty{3}(:,:,3,2) = [0   0.043933176436411   0.074998587258724   0.074998597814354
                                   0   0.043933230756944   0.074998133834873   0.074998675564354
                                   0   0.025735231037175   0.043933284336429   0.043933292548839
                                   0                   0                   0                   0];
expSolCPHistPenalty{3}(:,:,4,2) = [1.000000000000000   0.804737854124365   0.804737854124365   1.000000000000000
                                   0.804737854124365   0.647603013860688   0.647603013860688   0.804737854124365
                                   0.804737854124365   0.647603013860688   0.647603013860688   0.804737854124365
                                   1.000000000000000   0.804737854124365   0.804737854124365   1.000000000000000];
                                   
% Patch 4 :
% _________

expSolCPHistPenalty{4}(:,:,1,2) = [0                   -0.000000006530059  -0.000000010653681  -0.000000001746244
                                   -0.031066017177982  -0.031065554250159  -0.031065502116382  -0.031065470387551
                                   -0.075000000000000  -0.074998406242564  -0.074999279847258  -0.074999373701804
                                   -0.075000000000000  -0.075000000000000  -0.075000000000000  -0.075000000000000];
expSolCPHistPenalty{4}(:,:,2,2) = [   -0.075000000000000  -0.074999280336108  -0.031065398883400  -0.000000001948320
                                      -0.075000000000000  -0.074998759813926  -0.031065576996464  -0.000000010054504
                                      -0.031066017177982  -0.031065561066671  -0.012867865309055  -0.000000008030586
                                                       0                   0                   0                   0];
expSolCPHistPenalty{4}(:,:,3,2) = [0   0.031065444501950   0.074998592998946   0.074998597136038
                                   0   0.031065402041233   0.074998567961570   0.074998657893275
                                   0   0.012867630342151   0.031065499079002   0.031065524141687
                                   0                   0                   0                   0];
expSolCPHistPenalty{4}(:,:,4,2) = [1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000
                                   0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                                   0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                                   1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000];
                                                   
% Define the expected solution in terms of the residual history
expSolResHistoryPenalty = [    1.715150786432061
                               0.000089565287372
                               0.000000000011649
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
expSolMinElAreaSize = [ 0.008835393905557
                        0.001293936786475
                        0.008835393905557
                        0.001293936786475];

%% 8. Define the expected solutions for the Lagrange Multipliers method

% Define the expected solution in terms of Control Point displacements
expSolDispLM = [                   0
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
                   0.000000547209939
                   0.000001312790214
                  -0.000000561898614
                  -0.000001127890778
                   0.000002909720421
                  -0.000001095229781
                  -0.000002777349543
                   0.000000179623011
                  -0.000000448517542
                                   0
                                   0
                                   0
                  -0.000000672888419
                   0.000001105209563
                  -0.000002097080435
                  -0.000000743230127
                   0.000000612715573
                  -0.000001806248887
                  -0.000001373417341
                  -0.000000201994066
                  -0.000000900964308
                                   0
                                   0
                                   0
                   0.000000129771187
                   0.000000430060748
                  -0.000001926859725
                  -0.000000590868840
                   0.000000049682014
                  -0.000001239391211
                  -0.000000509592594
                  -0.000000571390066
                  -0.000000411047103
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
                   0.000000291105286
                  -0.000000801164770
                  -0.000000423769265
                  -0.000000553221263
                  -0.000001546031064
                  -0.000000732295892
                  -0.000001242736079
                  -0.000000200163173
                  -0.000000264925624
                                   0
                                   0
                                   0
                  -0.000000361547408
                  -0.000001022439295
                  -0.000002032554788
                  -0.000000716012214
                  -0.000000689460209
                  -0.000001698671092
                  -0.000000769277943
                   0.000000172486147
                  -0.000000473686811
                                   0
                                   0
                                   0
                   0.000000105933546
                  -0.000000464103592
                  -0.000001908885590
                  -0.000000403407956
                  -0.000000165763803
                  -0.000001520400804
                  -0.000000327675814
                   0.000000383018060
                  -0.000000220058079
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
                  -0.000000547209939
                  -0.000001312790214
                  -0.000000561898614
                   0.000001127890778
                  -0.000002909720421
                  -0.000001095229781
                   0.000002777349543
                  -0.000000179623011
                  -0.000000448517542
                                   0
                                   0
                                   0
                   0.000000672888419
                  -0.000001105209563
                  -0.000002097080435
                   0.000000743230127
                  -0.000000612715573
                  -0.000001806248887
                   0.000001373417341
                   0.000000201994066
                  -0.000000900964308
                                   0
                                   0
                                   0
                  -0.000000129771187
                  -0.000000430060748
                  -0.000001926859725
                   0.000000590868840
                  -0.000000049682014
                  -0.000001239391211
                   0.000000509592594
                   0.000000571390066
                  -0.000000411047103
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
                  -0.000000291105286
                   0.000000801164770
                  -0.000000423769265
                   0.000000553221263
                   0.000001546031064
                  -0.000000732295892
                   0.000001242736079
                   0.000000200163173
                  -0.000000264925624
                                   0
                                   0
                                   0
                   0.000000361547408
                   0.000001022439295
                  -0.000002032554788
                   0.000000716012214
                   0.000000689460209
                  -0.000001698671092
                   0.000000769277943
                  -0.000000172486147
                  -0.000000473686811
                                   0
                                   0
                                   0
                  -0.000000105933546
                   0.000000464103592
                  -0.000001908885590
                   0.000000403407956
                   0.000000165763803
                  -0.000001520400804
                   0.000000327675814
                  -0.000000383018060
                  -0.000000220058079
                                   0
                                   0
                                   0
                  -0.093998527118543
                  -0.171116839265407
                   0.205702249286556
                  -0.000273222941707
                   0.006810108902067
                   0.295065272906847
                  -0.018406737191299
                   0.156201548284016
                  -0.004642652030288
                   0.008603600522321
                   0.093998527118543
                   0.171116839265402
                   0.205702249286554
                  -0.000273222941707
                   0.006810108902067
                   0.295065272906846
                  -0.018406737191296
                  -0.156201548284020
                  -0.004642652030288
                   0.008603600522321];
              
% Define the expected solution in terms of the Control Point deformation
expSolCPHistLM = struc([]);
for counterPatches = 1:4
    for counterComp = 1:4
        expSolCPHistLM{counterPatches}(:,:,counterComp,1) = ...
            BSplinePatches{counterPatches}.CP(:,:,counterComp);
    end
end

% Patch 1 :
% _________

expSolCPHistLM{1}(:,:,1,2) = [                 0   0.000000547209933  -0.000000672888423   0.000000129771196
                               0.043933982822018   0.043932854931235   0.043933239591887   0.043933391953188
                               0.075000000000000   0.074997222650462   0.074998626582670   0.074999490407421
                               0.075000000000000   0.075000000000000   0.075000000000000   0.075000000000000];
expSolCPHistLM{1}(:,:,2,2) = [-0.075000000000000  -0.074998687209772  -0.043932877612445   0.000000430060757
                              -0.075000000000000  -0.074997090279579  -0.043933370106449   0.000000049682008
                              -0.043933982822018  -0.043933803199010  -0.025736133282142  -0.000000571390071
                                               0                   0                   0                   0];
expSolCPHistLM{1}(:,:,3,2) = [ 0   0.043933420923401   0.074997902919556   0.074998073140275
                               0   0.043932887592240   0.074998193751111   0.074998760608799
                               0   0.025735482770528   0.043933081857712   0.043933571774919
                               0                   0                   0                   0];
expSolCPHistLM{1}(:,:,4,2) = [ 1.000000000000000   0.804737854124365   0.804737854124365   1.000000000000000
                               0.804737854124365   0.647603013860688   0.647603013860688   0.804737854124365
                               0.804737854124365   0.647603013860688   0.647603013860688   0.804737854124365
                               1.000000000000000   0.804737854124365   0.804737854124365   1.000000000000000];

% Patch 2 :
% _________

expSolCPHistLM{2}(:,:,1,2) = [                 0   0.000000291105281  -0.000000361547409   0.000000105933555
                               0.031066017177982   0.031065463956716   0.031065301165769   0.031065613770036
                               0.075000000000000   0.074998757263923   0.074999230722065   0.074999672324196
                               0.075000000000000   0.075000000000000   0.075000000000000   0.075000000000000];
expSolCPHistLM{2}(:,:,2,2) = [ 0.075000000000000   0.074999198835221   0.031064994738676  -0.000000464103601
                               0.075000000000000   0.074998453968934   0.031065327717772  -0.000000165763802
                               0.031066017177982   0.031065817014810   0.012868138130186   0.000000383018063
                                               0                   0                   0                   0];
expSolCPHistLM{2}(:,:,3,2) = [ 0   0.031065593408715   0.074997967445205   0.074998091114409
                               0   0.031065284882090   0.074998301328908   0.074998479599203
                               0   0.012867700718411   0.031065543491173   0.031065797119906
                               0                   0                   0                   0];
expSolCPHistLM{2}(:,:,4,2) = [ 1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000
                               0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                               0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                               1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000];

% Patch 3 :
% _________

expSolCPHistLM{3}(:,:,1,2) = [                 0  -0.000000547209933   0.000000672888423  -0.000000129771196
                              -0.043933982822018  -0.043932854931235  -0.043933239591887  -0.043933391953188
                              -0.075000000000000  -0.074997222650462  -0.074998626582670  -0.074999490407421
                              -0.075000000000000  -0.075000000000000  -0.075000000000000  -0.075000000000000];
expSolCPHistLM{3}(:,:,2,2) = [ 0.075000000000000   0.074998687209772   0.043932877612445  -0.000000430060757
                               0.075000000000000   0.074997090279579   0.043933370106449  -0.000000049682008
                               0.043933982822018   0.043933803199010   0.025736133282142   0.000000571390071
                                               0                   0                   0                   0];
expSolCPHistLM{3}(:,:,3,2) = [ 0   0.043933420923401   0.074997902919556   0.074998073140275
                               0   0.043932887592240   0.074998193751111   0.074998760608799
                               0   0.025735482770528   0.043933081857712   0.043933571774919
                               0                   0                   0                   0];
expSolCPHistLM{3}(:,:,4,2) = [ 1.000000000000000   0.804737854124365   0.804737854124365   1.000000000000000
                               0.804737854124365   0.647603013860688   0.647603013860688   0.804737854124365
                               0.804737854124365   0.647603013860688   0.647603013860688   0.804737854124365
                               1.000000000000000   0.804737854124365   0.804737854124365   1.000000000000000];
                                   
% Patch 4 :
% _________

expSolCPHistLM{4}(:,:,1,2) = [                 0  -0.000000291105281   0.000000361547409  -0.000000105933555
                              -0.031066017177982  -0.031065463956716  -0.031065301165769  -0.031065613770036
                              -0.075000000000000  -0.074998757263923  -0.074999230722065  -0.074999672324196
                              -0.075000000000000  -0.075000000000000  -0.075000000000000  -0.075000000000000];
expSolCPHistLM{4}(:,:,2,2) = [-0.075000000000000  -0.074999198835221  -0.031064994738676   0.000000464103601
                              -0.075000000000000  -0.074998453968934  -0.031065327717772   0.000000165763802
                              -0.031066017177982  -0.031065817014810  -0.012868138130186  -0.000000383018063
                                               0                   0                   0                   0];
expSolCPHistLM{4}(:,:,3,2) = [ 0   0.031065593408715   0.074997967445205   0.074998091114409
                               0   0.031065284882090   0.074998301328908   0.074998479599203
                               0   0.012867700718411   0.031065543491173   0.031065797119906
                               0                   0                   0                   0];
expSolCPHistLM{4}(:,:,4,2) = [ 1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000
                               0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                               0.853553390593274   0.728553390593274   0.728553390593274   0.853553390593274
                               1.000000000000000   0.853553390593274   0.853553390593274   1.000000000000000];
                                                   

% Define the expected solution in terms of the residual history
expSolResHistoryLM = [ 1.715150786432061
                       0.000091487966429
                       0.000000000069326
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

%% 9. Define the nonlinear analysis parameters
propNLinearAnalysis.scheme = 'newtonRapshon';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 10e-6;
propNLinearAnalysis.maxIter = 120;

%% 10. Solve the multipatch problem with the Penalty method
propCouplingPenalty.method = 'penalty';
propCouplingPenalty.alphaD = zeros(connections.No,1);
propCouplingPenalty.alphaR = zeros(connections.No,1);
for iConnections = 1:connections.No
    propCouplingPenalty.alphaD(iConnections,1) = norm(Dm)/0.740369375376674/1.0e-04;
    propCouplingPenalty.alphaR(iConnections,1) = norm(Db)/0.740369375376674/1.0e-04*1e-1;
end
propCouplingPenalty.intC = intC;
plot_IGANLinear = '';
[dHatNonlinearPenalty,CPHistoryPenalty,resHistoryPenalty,minElASizePenalty] = ...
    solve_DDMIGAMembraneMultipatchesNLinear...
    (BSplinePatches,connections,propCouplingPenalty,propNLinearAnalysis,...
    solve_LinearSystem,plot_IGANLinear,[],'');

%% 11. Solve the multipatch problem with the Lagrange Multipliers method
propCouplingLM.method = 'lagrangeMultipliers';
propCouplingLM.alphaD = zeros(connections.No,1);
propCouplingLM.alphaR = zeros(connections.No,1);
for iConnections = connections.No
    propCouplingLM.alphaD(iConnections,1) = 0;
    propCouplingLM.alphaR(iConnections,1) = 0;
    propCouplingLM.intC = intC;
end
plot_IGANLinear = '';
[dHatNonlinearLM,CPHistoryLM,resHistoryLM,minElASizeLM] = ...
    solve_DDMIGAMembraneMultipatchesNLinear...
    (BSplinePatches,connections,propCouplingLM,propNLinearAnalysis,...
    solve_LinearSystem,plot_IGANLinear,[],'');

%% 12. Verify the solution
testCase.verifyEqual(dHatNonlinearPenalty,expSolDispPenalty,'AbsTol',absTol);
testCase.verifyEqual(CPHistoryPenalty,expSolCPHistPenalty,'AbsTol',absTol);
testCase.verifyEqual(resHistoryPenalty,expSolResHistoryPenalty,'AbsTol',absTol);
testCase.verifyEqual([BSplinePatches{1}.minElArea;BSplinePatches{2}.minElArea;BSplinePatches{3}.minElArea;BSplinePatches{4}.minElArea],expSolMinElAreaSize,'AbsTol',absTol);
testCase.verifyEqual(dHatNonlinearLM,expSolDispLM,'AbsTol',absTol);
testCase.verifyEqual(CPHistoryLM,expSolCPHistLM,'AbsTol',absTol);
testCase.verifyEqual(resHistoryLM,expSolResHistoryLM,'AbsTol',absTol);

end
