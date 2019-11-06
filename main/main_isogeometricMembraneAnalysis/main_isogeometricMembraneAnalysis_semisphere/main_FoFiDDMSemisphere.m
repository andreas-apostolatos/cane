%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Form-finding analysis over a 4-patch circular plate into a
%        semisphere and subsequently saving the results for further
%        steady-state or transient analysis.
%
% Date : 09.11.2016
%
%% Preamble
clear;
clc;

%% Includes 

% Add general math functions
addpath('../../../generalMath/');

% Add general auxiliary functions
addpath('../../../auxiliary/');

% Add system solvers
addpath('../../../equationSystemSolvers/');

% Add efficient computation functions
addpath('../../../efficientComputation/');

% Add transient analysis solvers
addpath('../../../transientAnalysis/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../../CAGDKernel/CAGDKernel_graphics/',...
        '../../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the isogeometric Kirchhoff-Love shell formulation
addpath('../../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../../isogeometricThinStructureAnalysis/graphicsMultipatches/',...
        '../../../isogeometricThinStructureAnalysis/loads/',...
        '../../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../../isogeometricThinStructureAnalysis/solvers/',...
        '../../../isogeometricThinStructureAnalysis/metrics/',...
        '../../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../../isogeometricThinStructureAnalysis/penaltyDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/penaltyDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../../isogeometricThinStructureAnalysis/precomputedData/',...
        '../../../isogeometricThinStructureAnalysis/output/',...
        '../../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/');
    
% Add functions related to the visualization of the multipatch geometry
% related to the isogeometric mortar-based mapping
addpath('../../../isogeometricMortarBasedMappingAnalysis/graphics/');
    
%% CAD modelling via NURBS

% Global variables
Radius = 0.075;

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
% CP1(:,:,1) = [0      0      0
%               Radius Radius Radius
%               Radius Radius Radius];
CP1(:,:,1) = [0      0        0
              Radius Radius/2 Radius
              Radius Radius   Radius];
          
% y-coordinates
% CP1(:,:,2) = [-Radius -Radius 0
%               -Radius -Radius 0
%               0       0       0];
CP1(:,:,2) = [-Radius -Radius/2 0
              -Radius -Radius/2 0
              0       0         0];
          
% z-coordinates for an exact semisphere
% CP1(:,:,3) = [0 Radius Radius
%               0 Radius Radius
%               0 0      0];
CP1(:,:,3) = [0 0 0
              0 0 0
              0 0 0];
          
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
% CP2(:,:,1) = [0      0      0
%               Radius Radius Radius
%               Radius Radius Radius];
CP2(:,:,1) = [0      0        0
              Radius Radius/2 Radius
              Radius Radius   Radius];
         
% y-coordinates
% CP2(:,:,2) = [Radius Radius 0
%               Radius Radius 0
%               0      0      0];
CP2(:,:,2) = [Radius Radius/2 0
              Radius Radius/2 0
              0      0        0];
         
% z-coordinates for an exact semisphere
% CP2(:,:,3) = [0 Radius Radius
%               0 Radius Radius
%               0 0      0];
CP2(:,:,3) = [0 0 0
              0 0 0
              0 0 0];
          
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
% CP3(:,:,1) = [0       0       0
%               -Radius -Radius -Radius
%               -Radius -Radius -Radius];
CP3(:,:,1) = [0       0         0
              -Radius -Radius/2 -Radius
              -Radius -Radius   -Radius];
         
% y-coordinates
% CP3(:,:,2) = [Radius Radius 0
%               Radius Radius 0
%               0      0      0];
CP3(:,:,2) = [Radius Radius/2 0
              Radius Radius/2 0
              0      0        0];
         
% z-coordinates for an exact semisphere
% CP3(:,:,3) = [0 Radius Radius
%               0 Radius Radius
%               0 0      0];
CP3(:,:,3) = [0 0 0
              0 0 0
              0 0 0];
       
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
% CP4(:,:,1) = [0       0       0
%               -Radius -Radius -Radius
%               -Radius -Radius -Radius];
CP4(:,:,1) = [0       0         0
              -Radius -Radius/2 -Radius
              -Radius -Radius   -Radius];
         
% y-coordinates
% CP4(:,:,2) = [-Radius -Radius 0
%               -Radius -Radius 0
%               0       0       0];
CP4(:,:,2) = [-Radius -Radius/2 0
              -Radius -Radius/2 0
              0       0         0];
         
% z-coordinates for an exact semisphere
% CP4(:,:,3) = [0 Radius Radius
%               0 Radius Radius
%               0 0      0];
CP4(:,:,3) = [0 0 0
              0 0 0
              0 0 0];
          
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
CPLambda12(:,4) = [1];

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
CPMu12(:,4) = [1];

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
CPLambda23(:,4) = [1];

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
CPMu23(:,4) = [1];

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
CPLambda34(:,4) = [1];

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
CPMu34(:,4) = [1];

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
CPLambda14(:,4) = [1];

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
CPMu14(:,4) = [1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBSMu14 = 0;
nxiMu14 = length(CPMu14(:,1,1));
for i=1:nxiMu14
    if CPMu14(i,4)~=1
        isNURBSMu14 = 1;
        break;
    end
end

%% Material constants

% Force amplitude
FAmp = +19;

% Parameters
EYoung = 7.0e5;
poissonRatio = .49;
thickness = 1.65e-4;
density = 1050;
prestress.computeParametricCoordinates = @(X) [X(1,1)^2 + X(2,1)^2
                                               atan(X(2,1)/X(1,1))];
prestress.computeBaseVectors = @(theta1,theta2) [cos(theta2) -sin(theta2)
                                                 sin(theta2) cos(theta2)
                                                 0           0];
prestress.voigtVector = [abs(FAmp)*Radius/2/thickness
                         abs(FAmp)*Radius/2/thickness
                         0];

% Patch 1 :
% _________

% Young's modulus
parameters1.E = EYoung;

% Poisson ratio
parameters1.nue = poissonRatio;

% Thickness of the plate
parameters1.t = thickness;

% Density of the membrane (used only for dynamics)
parameters1.rho = density;

% Prestress for the membrane
parameters1.prestress = prestress;

% Patch 2 :
% _________

% Young's modulus
parameters2.E = EYoung;

% Poisson ratio
parameters2.nue = poissonRatio;

% Thickness of the plate
parameters2.t = thickness;

% Density of the membrane (used only for dynamics)
parameters2.rho = density;

% Prestress for the membrane
parameters2.prestress = prestress;

% Patch 3 :
% _________

% Young's modulus
parameters3.E = EYoung;

% Poisson ratio
parameters3.nue = poissonRatio;

% Thickness of the plate
parameters3.t = thickness;

% Density of the membrane (used only for dynamics)
parameters3.rho = density;

% Prestress for the membrane
parameters3.prestress = prestress;

% Patch 4 :
% _________

% Young's modulus
parameters4.E = EYoung;

% Poisson ratio
parameters4.nue = poissonRatio;

% Thickness of the plate
parameters4.t = thickness;

% Density of the membrane (used only for dynamics)
parameters4.rho = density;

% Prestress for the membrane
parameters4.prestress = prestress;

%% GUI

% Case name
caseName = 'DDM4PatchesSemisphereistributedLoad';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points

% Patch 1 :
% _________

int1.type = 'default';
if strcmp(int1.type,'user')
    int1.xiNGP = 6;
    int1.etaNGP = 6;
    int1.xiNGPForLoad = 6;
    int1.etaNGPForLoad = 6;
    int1.nGPForLoad = 6;
    int1.nGPError = 12;
end

% Patch 2 :
% _________

int2.type = 'default';
if strcmp(int2.type,'user')
    int2.xiNGP = 6;
    int2.etaNGP = 6;
    int2.xiNGPForLoad = 6;
    int2.etaNGPForLoad = 6;
    int2.nGPForLoad = 6;
    int2.nGPError = 12;
end

% Patch 3 :
% _________

int3.type = 'default';
if strcmp(int3.type,'user')
    int3.xiNGP = 6;
    int3.etaNGP = 6;
    int3.xiNGPForLoad = 6;
    int3.etaNGPForLoad = 6;
    int3.nGPForLoad = 6;
    int3.nGPError = 12;
end

% Patch 4 :
% _________

int4.type = 'default';
if strcmp(int4.type,'user')
    int4.xiNGP = 6;
    int4.etaNGP = 6;
    int4.xiNGPForLoad = 6;
    int4.etaNGPForLoad = 6;
    int4.nGPForLoad = 6;
    int4.nGPError = 12;
end

% Interface integration :
% _______________________

intC.type = 'user';
intC.method = 'Nitsche';
if strcmp(intC.type,'user')
    if strcmp(intC.method,'lagrangeMultipliers')
        intC.nGP1 = 16;
        intC.nGP2 = 16;
    else
        intC.noGPs = 16;
    end
    intC.nGPError = 16;
end

% Patch 1 :
% _________

Dm1 = parameters1.E*parameters1.t/(1-parameters1.nue^2)*...
      [1              parameters1.nue 0
       parameters1.nue 1              0
       0               0              (1-parameters1.nue)/2];
Db1 = parameters1.t^2/12*Dm1;
   
% Patch 2 :
% _________

Dm2 = parameters2.E*parameters2.t/(1-parameters2.nue^2)*...
      [1              parameters2.nue 0
       parameters2.nue 1              0
       0               0              (1-parameters2.nue)/2];
Db2 = parameters2.t^2/12*Dm2;

% Patch 3 :
% _________

Dm3 = parameters3.E*parameters3.t/(1-parameters3.nue^2)*...
      [1              parameters3.nue 0
       parameters3.nue 1              0
       0               0              (1-parameters3.nue)/2];
Db3 = parameters3.t^2/12*Dm3;

% Patch 4 :
% _________

Dm4 = parameters4.E*parameters4.t/(1-parameters4.nue^2)*...
      [1              parameters4.nue 0
       parameters4.nue 1              0
       0               0              (1-parameters4.nue)/2];
Db4 = parameters4.t^2/12*Dm4;

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','force'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x','y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

%% Refinement

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = 2; % 2
tp1 = a;
tq1 = a;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% Patch 2 :
% _________

b = 1; % 1
tp2 = b;
tq2 = b;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');

% Patch 3 :
% _________

c = 2; % 2
tp3 = c;
tq3 = c;
[Xi3,Eta3,CP3,p3,q3] = degreeElevateBSplineSurface(p3,q3,Xi3,Eta3,CP3,tp3,tq3,'');

% Patch 4 :
% _________

d = 1; % 1
tp4 = d;
tq4 = d;
[Xi4,Eta4,CP4,p4,q4] = degreeElevateBSplineSurface(p4,q4,Xi4,Eta4,CP4,tp4,tq4,'');

%%% Coupling between patch 1 and patch 2 %%%
%%% ------------------------------------ %%%

pLM = min(p1,p2);
if pLM > 0
    clear pLambda12 XiLambda12 CPLambda12;

    % Polynomial order
    pLambda12 = 1;

    % Knot vector
    XiLambda12 = [0 0 1 1];

    % Control Point coordinates and weights
    CPLambda12(:,4) = [1 1];

    % Find whether the basis is NURBS or a B-Spline
    isNURBSLambda12 = 0;
    nxiLambda12 = length(CPLambda12(:,1,1));
    for i = 1:nxiLambda12
        if CPLambda12(i,4) ~= 1
            isNURBSLambda12 = 1;
            break;
        end
    end

    % Perform accordingly a p-refinement
    tpLambda12 = pLM;
    [XiLambda12,CPLambda12,pLambda12] = degreeElevateBSplineCurve...
        (pLambda12,XiLambda12,CPLambda12,tpLambda12,'');
end

%%% Coupling between patch 2 and patch 3 %%%
%%% ------------------------------------ %%%

pLM = min(p2,p3);
if pLM > 0
    clear pLambda23 XiLambda23 CPLambda23;

    % Polynomial order
    pLambda23 = 1;

    % Knot vector
    XiLambda23 = [0 0 1 1];

    % Control Point coordinates and weights
    CPLambda23(:,4) = [1 1];

    % Find whether the basis is NURBS or a B-Spline
    isNURBSLambda23 = 0;
    nxiLambda23 = length(CPLambda23(:,1,1));
    for i = 1:nxiLambda23
        if CPLambda23(i,4) ~= 1
            isNURBSLambda23 = 1;
            break;
        end
    end

    % Perform accordingly a p-refinement
    tpLambda23 = pLM;
    [XiLambda23,CPLambda23,pLambda23] = degreeElevateBSplineCurve...
        (pLambda23,XiLambda23,CPLambda23,tpLambda23,'');
end

%%% Coupling between patch 3 and patch 4 %%%
%%% ------------------------------------ %%%

pLM = min(p3,p4);
if pLM > 0
    clear pLambda34 XiLambda34 CPLambda34;

    % Polynomial order
    pLambda34 = 1;

    % Knot vector
    XiLambda34 = [0 0 1 1];

    % Control Point coordinates and weights
    CPLambda34(:,4) = [1 1];

    % Find whether the basis is NURBS or a B-Spline
    isNURBSLambda34 = 0;
    nxiLambda34 = length(CPLambda34(:,1,1));
    for i = 1:nxiLambda34
        if CPLambda34(i,4) ~= 1
            isNURBSLambda34 = 1;
            break;
        end
    end

    % Perform accordingly a p-refinement
    tpLambda34 = pLM;
    [XiLambda34,CPLambda34,pLambda34] = degreeElevateBSplineCurve...
        (pLambda34,XiLambda34,CPLambda34,tpLambda34,'');
end

%%% Coupling between patch 1 and patch 4 %%%
%%% ------------------------------------ %%%

pLM = min(p1,p4);
if pLM > 0
    clear pLambda14 XiLambda14 CPLambda14;

    % Polynomial order
    pLambda14 = 1;

    % Knot vector
    XiLambda14 = [0 0 1 1];

    % Control Point coordinates and weights
    CPLambda14(:,4) = [1 1];

    % Find whether the basis is NURBS or a B-Spline
    isNURBSLambda14 = 0;
    nxiLambda14 = length(CPLambda14(:,1,1));
    for i = 1:nxiLambda14
        if CPLambda14(i,4) ~= 1
            isNURBSLambda14 = 1;
            break;
        end
    end

    % Perform accordingly a p-refinement
    tpLambda14 = pLM;
    [XiLambda14,CPLambda14,pLambda14] = degreeElevateBSplineCurve...
        (pLambda14,XiLambda14,CPLambda14,tpLambda14,'');
end

%%%%%%%%%%%%%%%%%%%
% Knot insertion  %
%%%%%%%%%%%%%%%%%%%

% Refinement parameters
noHRefStart = 8; % 8
iHRef = 8;

% Patch 1 :
% _________

noKnotsXi1 = ceil(5/9*(iHRef + noHRefStart - 1));
% noKnotsXi1 = 0;
noKnotsEta1 = noKnotsXi1;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface...
    (p1,Xi1,q1,Eta1,CP1,noKnotsXi1,noKnotsEta1,'');

% Patch 2 :
% _________

noKnotsXi2 = ceil(6/9*(iHRef + noHRefStart - 1));
% noKnotsXi2 = 0;
noKnotsEta2 = noKnotsXi2;
[Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface...
    (p2,Xi2,q2,Eta2,CP2,noKnotsXi2,noKnotsEta2,'');

% Patch 3 :
% _________

noKnotsXi3 = ceil(4/9*(iHRef + noHRefStart - 1));
% noKnotsXi3 = 0;
noKnotsEta3 = noKnotsXi3;
[Xi3,Eta3,CP3] = knotRefineUniformlyBSplineSurface...
    (p3,Xi3,q3,Eta3,CP3,noKnotsXi3,noKnotsEta3,'');

% Patch 4 :
% _________

noKnotsXi4 = ceil(7/9*(iHRef + noHRefStart - 1));
% noKnotsXi4 = 0;
noKnotsEta4 = noKnotsXi4;
[Xi4,Eta4,CP4] = knotRefineUniformlyBSplineSurface...
    (p4,Xi4,q4,Eta4,CP4,noKnotsXi4,noKnotsEta4,'');

% Scale the discretization of the Lagrange Multipliers field
% according to the discretization of the interface of the
% neighbouring patches
scaleLM = 1.0;
scaleALM = .6;

%%% Coupling between patch 1 and patch 2 %%%
%%% ------------------------------------ %%%

noLambda12 = ceil(min(noKnotsEta1,noKnotsEta2)*scaleLM);
[XiLambda12,CPLambda12] = knotRefineUniformlyBSplineCurve...
    (noLambda12,pLambda12,XiLambda12,CPLambda12,'');

%%% Coupling between patch 2 and patch 3 %%%
%%% ------------------------------------ %%%

noLambda23 = ceil(min(noKnotsEta2,noKnotsEta3)*scaleLM);
[XiLambda23,CPLambda23] = knotRefineUniformlyBSplineCurve...
    (noLambda23,pLambda23,XiLambda23,CPLambda23,'');

%%% Coupling between patch 3 and patch 4 %%%
%%% ------------------------------------ %%%

noLambda34 = ceil(min(noKnotsEta3,noKnotsEta4)*scaleLM);
[XiLambda34,CPLambda34] = knotRefineUniformlyBSplineCurve...
    (noLambda34,pLambda34,XiLambda34,CPLambda34,'');

%%% Coupling between patch 1 and patch 4 %%%
%%% ------------------------------------ %%%

noLambda14 = ceil(min(noKnotsEta1,noKnotsEta4)*scaleLM);
[XiLambda14,CPLambda14] = knotRefineUniformlyBSplineCurve...
    (noLambda14,pLambda14,XiLambda14,CPLambda14,'');

%% Fill up the arrays for the Lagrange Multipliers

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

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous DOFs
homDOFs1 = [];
% xiSup1 = [0 1];
% etaSup1 = [0 0];
xiSup1 = {[1 1]};
etaSup1 = {[0 1]};
noCnd1 = length(xiSup1);
for iCnd = 1:noCnd1
    for dirSupp = 1:3
        homDOFs1 = findDofs3D(homDOFs1,xiSup1{iCnd},etaSup1{iCnd},dirSupp,CP1);
    end
end

% Inhomogeneous DOFs
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak boundary conditions
weakDBC1.noCnd = 1;
weakDBC1.xiExtension = {[0 1]};
weakDBC1.etaExtension = {[0 0]};
weakDBC1.int.type = 'default';
weakDBC1.int.noGPs = 16;

% Embedded cables
cables1.No = 0;

% Patch 2 :
% _________

% Homogeneous DOFs
homDOFs2 = [];
% xiSup2 = [0 1];
% etaSup2 = [0 0];
xiSup2 = {[1 1]};
etaSup2 = {[0 1]};
noCnd2 = length(xiSup2);
for iCnd = 1:noCnd2
    for dirSupp = 1:3
        homDOFs2 = findDofs3D(homDOFs2,xiSup2{iCnd},etaSup2{iCnd},dirSupp,CP2);
    end
end

% Inhomogeneous DOFs
inhomDOFs2 = [];
valuesInhomDOFs2 = [];

% Weak boundary conditions
weakDBC2.noCnd = 1;
weakDBC2.xiExtension = {[0 1]};
weakDBC2.etaExtension = {[0 0]};
weakDBC2.int.type = 'default';
weakDBC2.int.noGPs = 16;

% Embedded cables
cables2.No = 0;

% Patch 3 :
% _________

% Homogeneous DOFs
homDOFs3 = [];
% xiSup3 = [0 1];
% etaSup3 = [0 0];
xiSup3 = {[1 1]};
etaSup3 = {[0 1]};
noCnd3 = length(xiSup3);
for iCnd = 1:noCnd3
    for dirSupp = 1:3
        homDOFs3 = findDofs3D(homDOFs3,xiSup3{iCnd},etaSup3{iCnd},dirSupp,CP3);
    end
end

% Inhomogeneous DOFs
inhomDOFs3 = [];
valuesInhomDOFs3 = [];

% Weak boundary conditions
weakDBC3.noCnd = 1;
weakDBC3.xiExtension = {[0 1]};
weakDBC3.etaExtension = {[0 0]};
weakDBC3.int.type = 'default';
weakDBC3.int.noGPs = 16;

% Embedded cables
cables3.No = 0;

% Patch 4 :
% _________

% Homogeneous DOFs
homDOFs4 = [];
% xiSup4 = [0 1];  
% etaSup4 = [0 0];
xiSup4 = {[1 1]};
etaSup4 = {[0 1]};
noCnd4 = length(xiSup4);
for iCnd = 1:noCnd4
    for dirSupp = 1:3
        homDOFs4 = findDofs3D(homDOFs4,xiSup4{iCnd},etaSup4{iCnd},dirSupp,CP4);
    end
end

% Inhomogeneous DOFs
inhomDOFs4 = [];
valuesInhomDOFs4 = [];

% Weak boundary conditions
weakDBC4.noCnd = 1;
weakDBC4.xiExtension = {[0 1]};
weakDBC4.etaExtension = {[0 0]};
weakDBC4.int.type = 'default';
weakDBC4.int.noGPs = 16;

% Embedded cables
cables4.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loadAmplitude = FAmp;
FoFiStep = -1;
loadAmplitude1 = @(x,y,z,t) loadAmplitude*(t<=FoFiStep)*sin(pi/2*t/FoFiStep) + loadAmplitude*(t>FoFiStep);
loadAmplitude2 = @(x,y,z,t) - loadAmplitude*(t<=FoFiStep)*sin(pi/2*t/FoFiStep) - loadAmplitude*(t>FoFiStep);
dirForce = 'normal';
isFollower = true;

% Patch 1 :
% _________

NBC1.noCnd = 1;
xib1 = [0 1];   etab1 = [0 1];
NBC1.xiLoadExtension = {xib1};
NBC1.etaLoadExtension = {etab1};
NBC1.loadAmplitude = {loadAmplitude1};
% NBC1.loadAmplitude = {@(x,y,z,t) (t<=5)*FAmp + (t>5)*FAmp*9e-1};
NBC1.loadDirection = {dirForce};
NBC1.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC1.isFollower(1,1) = isFollower;
NBC1.isTimeDependent(1,1) = false;
        
% Patch 2 :
% _________

NBC2.noCnd = 1;
xib2 = [0 1];   etab2 = [0 1];
NBC2.xiLoadExtension = {xib2};
NBC2.etaLoadExtension = {etab2};
NBC2.loadAmplitude = {loadAmplitude2};
% NBC2.loadAmplitude = {@(x,y,z,t) - (t<=5)*FAmp - (t>5)*FAmp*9e-1};
NBC2.loadDirection = {dirForce};
NBC2.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC2.isFollower(1,1) = isFollower;
NBC2.isTimeDependent(1,1) = false;

% Patch 3 :
% _________

NBC3.noCnd = 1;
xib3 = [0 1];   etab3 = [0 1];
NBC3.xiLoadExtension = {xib3};
NBC3.etaLoadExtension = {etab3};
NBC3.loadAmplitude = {loadAmplitude1};
% NBC3.loadAmplitude = {@(x,y,z,t) (t<=5)*FAmp + (t>5)*FAmp*9e-1};
NBC3.loadDirection = {dirForce};
NBC3.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC3.isFollower(1,1) = isFollower;
NBC3.isTimeDependent(1,1) = false;

% Patch 4 :
% _________

NBC4.noCnd = 1;
xib4 = [0 1];   etab4 = [0 1];
NBC4.xiLoadExtension = {xib4};
NBC4.etaLoadExtension = {etab4};
NBC4.loadAmplitude = {loadAmplitude2};
% NBC4.loadAmplitude = {@(x,y,z,t) - (t<=5)*FAmp - (t>5)*FAmp*9e-1};
NBC4.loadDirection = {dirForce};
NBC4.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC4.isFollower(1,1) = isFollower;
NBC4.isTimeDependent(1,1) = false;

% Collect all the Neumann boundary conditions into an arra<y
NBC = {NBC1 NBC2 NBC3 NBC4};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface parametrization     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Connection with patch 2:
xicoup12 = [0 1];
etacoup12 = [1 1];

% Connection with patch 4:
xicoup14 = [0 0];
etacoup14 = [0 1];

% Collect all interfaces into arrays:
xicoup1 = [xicoup12
           xicoup14];
etacoup1 = [etacoup12
            etacoup14];

% Patch 2 :
% _________

% Connection with patch 1:
xicoup21 = [0 1];
etacoup21 = [1 1];

% Connection with patch 3:
xicoup23 = [0 0];
etacoup23 = [0 1];

% Collect all interfaces into arrays:
xicoup2 = [xicoup21
           xicoup23];
etacoup2 = [etacoup21
            etacoup23];

% Patch 3 :
% _________

% connection with patch 2:
xicoup32 = [0 0];
etacoup32 = [0 1];

% Connection with patch 4:
xicoup34 = [0 1];
etacoup34 = [1 1];

% Collect all interfaces into arrays:
xicoup3 = [xicoup32
           xicoup34];
etacoup3 = [etacoup32
            etacoup34];

% Patch 4 :
% _________

% Connection with patch 3:
xicoup43 = [0 1];
etacoup43 = [1 1];

% Connection with patch 1:
xicoup41 = [0 0];
etacoup41 = [0 1];

% Collect all interfaces into arrays:
xicoup4 = [xicoup43
           xicoup41];
etacoup4 = [etacoup43
            etacoup41];

% Define connections by patch numbers
connections.xiEtaCoup(:,:) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21
                              2 3 xicoup23 etacoup23 xicoup32 etacoup32
                              3 4 xicoup34 etacoup34 xicoup43 etacoup43
                              4 1 xicoup41 etacoup41 xicoup14 etacoup14];
connections.No = length(connections.xiEtaCoup(:,1));

%% Create the patches and the Lagrange Multiplier fields

% 1st patch :
% ___________

patch1 = fillUpPatch...
    (analysis,p1,Xi1,q1,Eta1,CP1,isNURBS1,parameters1,homDOFs1,...
    inhomDOFs1,valuesInhomDOFs1,weakDBC1,cables1,NBC1,[],[],[],xicoup1,etacoup1,int1);

% 2nd patch :
% ___________

patch2 = fillUpPatch...
    (analysis,p2,Xi2,q2,Eta2,CP2,isNURBS2,parameters2,homDOFs2,...
    inhomDOFs2,valuesInhomDOFs2,weakDBC2,cables2,NBC2,[],[],[],xicoup2,etacoup2,int2);

% 3rd patch :
% ___________

patch3 = fillUpPatch...
    (analysis,p3,Xi3,q3,Eta3,CP3,isNURBS3,parameters3,homDOFs3,...
    inhomDOFs3,valuesInhomDOFs3,weakDBC3,cables3,NBC3,[],[],[],xicoup3,etacoup3,int3);

% 4th patch :
% ___________

patch4 = fillUpPatch...
    (analysis,p4,Xi4,q4,Eta4,CP4,isNURBS4,parameters4,homDOFs4,...
    inhomDOFs4,valuesInhomDOFs4,weakDBC4,cables4,NBC4,[],[],[],xicoup4,etacoup4,int4);

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch1 patch2 patch3 patch4};
noPatches = length(BSplinePatches);

% Collect all the Lagrange Multipliers into the connections :
% __________________________________________________________

connectionsLM = connections;
connectionsLM.lambda = {lambda12 lambda23 lambda34 lambda14};
% connectionsLM.mu = {mu12 mu23 mu34 mu14};

%% Plot the distribution of the determinant of the geometrical Jacobian for the multipatch geometry
% figure(graph.index)
% for iPatches = 1:4
%     plot_postprocBSplineSurfaceGeometricJacobian(BSplinePatches{iPatches}.p,...
%         BSplinePatches{iPatches}.q,BSplinePatches{iPatches}.Xi,...
%         BSplinePatches{iPatches}.Eta,BSplinePatches{iPatches}.CP,...
%         BSplinePatches{iPatches}.isNURBS);
%     hold on;
% end
% hold off;
% shading interp;
% colormap('jet');
% camlight left;
% lighting phong;
% colorbar;
% axis equal;
% graph.index = graph.index + 1;

%% Plot the multipatch geometry with the patch numbering
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,graph);

%% Compute the load vector for the visualization of the reference configuration
for counterPatches = 1:noPatches
    BSplinePatches{counterPatches}.FGamma = ...
        zeros(3*BSplinePatches{counterPatches}.noCPs,1);
%     for counterNBC = 1:NBC{counterPatches}.noCnd
%         funcHandle = str2func(NBC{counterPatches}.computeLoadVct{counterNBC});
%         BSplinePatches{counterPatches}.FGamma = funcHandle...
%             (BSplinePatches{counterPatches}.FGamma,...
%             BSplinePatches{counterPatches},...
%             NBC{counterPatches}.xiLoadExtension{counterNBC},...
%             NBC{counterPatches}.etaLoadExtension{counterNBC},...
%             NBC{counterPatches}.loadAmplitude{counterNBC},...
%             NBC{counterPatches}.loadDirection{counterNBC},...
%             NBC{counterPatches}.isFollower(counterNBC,1),0,...
%             BSplinePatches{counterPatches}.int,'outputEnabled');
%     end
end

%% Plot the reference configuration for the multipatch geometry before the form-finding analysis
color = [.85098 .8549 .85882];
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
% save data_FoFiDDMExactSemisphere;
% return;

%% Perform a form finding analysis to find the shape of the multipatch membrane with the Nitsche method

% Properties for the form-finding analysis
% ----------------------------------------

propFormFinding.tolerance = 1e-4;
propFormFinding.maxNoIter = 400;
propFormFinding.minNoIter = 1;
penaltyPrmScale = 1e6;

% Assign the parameters for the application of weak DBC
% -----------------------------------------------------

% Assing the B-Spline patches for the Penalty method
if exist('BSplinePatchesNitsche','var');
    clear BSplinePatchesNitsche;
end
BSplinePatchesNitsche = BSplinePatches;

% Properties for the weak Dirichlet boundary conditions
for iPatches = 1:noPatches
    BSplinePatchesNitsche{iPatches}.weakDBC.method = 'nitsche';
    BSplinePatchesNitsche{iPatches}.weakDBC.estimationStabilPrm = true;
end

% Assign the parameters for multipatch coupling
% ---------------------------------------------

propCouplingNitsche.method = 'nitsche';
propCouplingNitsche.estimationStabilPrm = true;
propCouplingNitsche.gammaTilde = .5;
propCouplingNitsche.intC = intC;

% Solve the problem using the Penalty method
% ------------------------------------------

[BSplinePatches,CPHistoryMultipatch,propCouplingNitsche,resHistory,...
    hasConverged,noIter] = ...
    solve_DDMFormFindingIGAMembrane(BSplinePatchesNitsche,connections,...
    propCouplingNitsche,propFormFinding,solve_LinearSystem,'outputEnabled');
save data_FoFiDDMSemisphere;

%% Perform a form finding analysis to find the shape of the multipatch membrane with the Penalty method

% % Properties for the form-finding analysis
% % ----------------------------------------
% 
% propFormFinding.tolerance = 1e-6;
% propFormFinding.maxNoIter = 400;
% propFormFinding.minNoIter = 1;
% penaltyPrmScale = 1e0;
% 
% % Assign the parameters for the application of weak DBC
% % -----------------------------------------------------
% 
% % Assing the B-Spline patches for the Penalty method
% if exist('BSplinePatchesPenalty','var');
%     clear BSplinePatchesPenalty;
% end
% BSplinePatchesPenalty = BSplinePatches;
% 
% % Properties for the weak Dirichlet boundary conditions
% for iPatches = 1:noPatches
%     % Assign the method name
%     BSplinePatchesPenalty{iPatches}.weakDBC.method = 'penalty';
% 
%     % Get the polynomial order along the Dirichlet boundary
%     isOnXi = false;
%     if BSplinePatchesPenalty{iPatches}.weakDBC.etaExtension{1}(1,1) == ...
%             BSplinePatchesPenalty{iPatches}.weakDBC.etaExtension{1}(1,2)
%         isOnXi = true;
%     end
%     if isOnXi
%         polOrder = BSplinePatchesPenalty{iPatches}.p;
%     else
%         polOrder = BSplinePatchesPenalty{iPatches}.q;
%     end
% 
%     % Assign the penalty parameter
%     BSplinePatchesPenalty{iPatches}.weakDBC.alpha = ...
%         norm(eval(['Dm' num2str(iPatches)]))*polOrder/...
%         BSplinePatches{iPatches}.minElArea*...
%         penaltyPrmScale;
% end
% 
% % Assign the parameters for multipatch coupling
% % ---------------------------------------------
% 
% propCouplingPenalty.method = 'penalty';
% propCouplingPenalty.intC = intC;
% propCouplingPenalty.alphaD = zeros(connections.No,1);
% propCouplingPenalty.alphaR = zeros(connections.No,1);
% for iConnections = 1:connections.No
%     % Get the id's of the patches
%     IDPatchI = connections.xiEtaCoup(iConnections,1);
%     IDPatchJ = connections.xiEtaCoup(iConnections,2);
%     
%     % Get the connection information
%     xicoupIJ = connections.xiEtaCoup(iConnections,3:4);
%     etacoupIJ = connections.xiEtaCoup(iConnections,5:6);
%     xicoupJI = connections.xiEtaCoup(iConnections,7:8);
%     etacoupJI = connections.xiEtaCoup(iConnections,9:10);
% 
%     % Get the mean polynomial order between the patches
%     isOnXiI = false;
%     if etacoupIJ(1,1) == etacoupIJ(1,2)
%         isOnXiI = true;
%     end
%     if isOnXiI
%         polOrderI = BSplinePatchesPenalty{IDPatchI}.p;
%     else
%         polOrderI = BSplinePatchesPenalty{IDPatchI}.q;
%     end
%     isOnXiJ = false;
%     if etacoupJI(1,1) == etacoupJI(1,2)
%         isOnXiJ = true;
%     end
%     if isOnXiJ
%         polOrderJ = BSplinePatchesPenalty{IDPatchJ}.p;
%     else
%         polOrderJ = BSplinePatchesPenalty{IDPatchJ}.q;
%     end
%     polOrderMean = mean([polOrderI polOrderJ]);
% 
%     % Assign the penalty parameters
%     propCouplingPenalty.alphaD(iConnections,1) = ...
%         max([norm(eval(['Dm' num2str(IDPatchI)])) ...
%         norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
%         min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea])*...
%         penaltyPrmScale;
%     propCouplingPenalty.alphaR(iConnections,1) = 0;
% end
% 
% % Solve the problem using the Penalty method
% % ------------------------------------------
% 
% [BSplinePatches,CPHistoryMultipatch,propCoupling,resHistory,...
%     hasConverged,noIter] = ...
%     solve_DDMFormFindingIGAMembrane(BSplinePatchesPenalty,connections,...
%     propCouplingPenalty,propFormFinding,solve_LinearSystem,'outputEnabled');
% save data_FoFiDDMSemisphere;

%% Perform a form finding analysis to find the shape of the multipatch membrane with the Lagrange Multipliers method

% % Scale the penalty parameters
% penaltyPrmScale = 0;%1e-2;
% 
% % Properties for the form-finding analysis
% % ----------------------------------------
% 
% propFormFinding.tolerance = 1e-5;
% propFormFinding.maxNoIter = 1e4;
% propFormFinding.minNoIter = 1;
% 
% % Assign the parameters for the application of weak DBC
% % -----------------------------------------------------
% 
% % Assing the B-Spline patches for the Penalty method
% if exist('BSplinePatchesLM','var');
%     clear BSplinePatchesLM;
% end
% BSplinePatchesLM = BSplinePatches;
% 
% % Properties for the weak Dirichlet boundary conditions
% for iPatches = 1:noPatches
%     BSplinePatchesLM{iPatches}.weakDBC.method = 'lagrangeMultipliers';
%     
%     % Get the polynomial order along the Dirichlet boundary
%     isOnXi = false;
%     if BSplinePatchesLM{iPatches}.weakDBC.etaExtension{1}(1,1) == ...
%             BSplinePatchesLM{iPatches}.weakDBC.etaExtension{1}(1,2)
%         isOnXi = true;
%     end
%     if isOnXi
%         polOrder = BSplinePatchesLM{iPatches}.p;
%     else
%         polOrder = BSplinePatchesLM{iPatches}.q;
%     end
% 
%     % Assign the penalty parameter
%     BSplinePatchesLM{iPatches}.weakDBC.alpha = ...
%         norm(eval(['Dm' num2str(iPatches)]))*polOrder/...
%         BSplinePatches{iPatches}.minElArea*...
%         penaltyPrmScale;
% 
%     % Find along which parametric line the weak Dirichlet 
%     % condition is to be imposed
%     isOnXi = false;
%     if BSplinePatchesLM{iPatches}.weakDBC.xiExtension{1}(1,1) == ...
%          BSplinePatchesLM{iPatches}.weakDBC.xiExtension{1}(1,2)
%         isOnXi = true;
%     end
% 
%     % Make a Lagrange Multipliers discretization
%     clear pLambda XiLambda CPLambda isNURBSLambda; 
% 
%     % Find out up to which polynomial degree the Lagrange
%     % Multipliers discretization needs to be increased
%     if isOnXi
%         polOrderPatch =  BSplinePatchesLM{iPatches}.p;
%     else
%         polOrderPatch =  BSplinePatchesLM{iPatches}.q;
%     end
%     pLM = polOrderPatch - 1;
% 
%     if pLM <= 0
%         pLambda = 0;
%         XiLambda = [0 1];
%         CPLambda = zeros(1,4);
%     else
%         pLambda = 1;
%         XiLambda = [0 0 1 1];
%         CPLambda = zeros(2,4);
% 
%         % Peform a p-refinement
%         tpLambda = polOrderPatch - pLM;
%         [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
%             (pLambda,XiLambda,CPLambda,tpLambda,'');
%     end
%     isNURBSLambda = 0;
%     nxiLambda = length(CPLambda(:,1));
%     for i = 1:nxiLambda
%         if CPLambda(i,4) ~= 1
%             isNURBSLambda = 1;
%             break;
%         end
%     end
% 
%     % Perform an h-refinement
%     percentage = 1.0;
%     if isOnXi
%         Rxi = unique(BSplinePatchesLM{iPatches}.Xi);
%     else
%         Rxi = unique(BSplinePatchesLM{iPatches}.Eta);
%     end
%     noXi = ceil(percentage*(length(Rxi) - 2));
%     [XiLambda,CPLambda] = knotRefineUniformlyBSplineCurve...
%         (noXi,pLambda,XiLambda,CPLambda,'');
% 
%     % Create the corresponding Lagrange Multipliers structure
%     BSplinePatchesLM{iPatches}.weakDBC.lambda{1} = fillUpLagrangeMultipliers...
%      (pLambda,XiLambda,CPLambda,isNURBSLambda);
% end
% 
% % Assign the parameters for multipatch coupling
% % ---------------------------------------------
% 
% propCouplingLM.method = 'lagrangeMultipliers';
% propCouplingLM.alphaD = zeros(connectionsLM.No,1);
% propCouplingLM.alphaR = zeros(connectionsLM.No,1);
% propCouplingLM.intC = intC;
% for iConnections = 1:connectionsLM.No
%     % Get the id's of the patches
%     IDPatchI = connectionsLM.xiEtaCoup(iConnections,1);
%     IDPatchJ = connectionsLM.xiEtaCoup(iConnections,2);
% 
%     % Get the connection information
%     xicoupIJ = connectionsLM.xiEtaCoup(iConnections,3:4);
%     etacoupIJ = connectionsLM.xiEtaCoup(iConnections,5:6);
%     xicoupJI = connectionsLM.xiEtaCoup(iConnections,7:8);
%     etacoupJI = connectionsLM.xiEtaCoup(iConnections,9:10);
%     
%     % Get the mean polynomial order between the patches
%     isOnXiI = false;
%     if etacoupIJ(1,1) == etacoupIJ(1,2)
%         isOnXiI = true;
%     end
%     if isOnXiI
%         polOrderI = BSplinePatchesLM{IDPatchI}.p;
%     else
%         polOrderI = BSplinePatchesLM{IDPatchI}.q;
%     end
%     isOnXiJ = false;
%     if etacoupJI(1,1) == etacoupJI(1,2)
%         isOnXiJ = true;
%     end
%     if isOnXiJ
%         polOrderJ = BSplinePatchesLM{IDPatchJ}.p;
%     else
%         polOrderJ = BSplinePatchesLM{IDPatchJ}.q;
%     end
%     polOrderMean = mean([polOrderI polOrderJ]);
% 
%     % Assign the penalty parameters
%     propCouplingLM.alphaD(iConnections,1) = ...
%         max([norm(eval(['Dm' num2str(IDPatchI)])) ...
%         norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
%         min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea])*...
%         penaltyPrmScale;
%     propCouplingLM.alphaR(iConnections,1) = 0;
% end
% fprintf('\n');
% 
% % Solve the problem using the Lagrange Multipliers method
% % -------------------------------------------------------
% 
% [BSplinePatches,CPHistoryMultipatch,propCouplingLM,resHistory,...
%     hasConverged,noIter] = solve_DDMFormFindingIGAMembrane...
%     (BSplinePatchesLM,connectionsLM,propCouplingLM,propFormFinding,...
%     solve_LinearSystem,'outputEnabled');
% save data_FoFiDDMSemisphere;

%% Compute the load vector for the visualization of the reference configuration
for counterPatches = 1:noPatches
    BSplinePatches{counterPatches}.FGamma = ...
        zeros(3*BSplinePatches{counterPatches}.noCPs,1);
%     for counterNBC = 1:NBC{counterPatches}.noCnd
%         funcHandle = str2func(NBC{counterPatches}.computeLoadVct{counterNBC});
%         BSplinePatches{counterPatches}.FGamma = funcHandle...
%             (BSplinePatches{counterPatches}.FGamma,...
%             BSplinePatches{counterPatches},...
%             NBC{counterPatches}.xiLoadExtension{counterNBC},...
%             NBC{counterPatches}.etaLoadExtension{counterNBC},...
%             NBC{counterPatches}.loadAmplitude{counterNBC},...
%             NBC{counterPatches}.loadDirection{counterNBC},...
%             NBC{counterPatches}.isFollower(counterNBC,1),0,...
%             BSplinePatches{counterPatches}.int,'outputEnabled');
%     end
end

%% Plot the deviation from an exact semisphere
propShpere.center = [0;0;0];
propShpere.radius = Radius;
int.type = 'default';
dHat = 'undefined';
[relErrorInL2Domain,relErrorInL2WeakDBC,errorInL2Interface,graph.index] = ...
    plot_deviationFromSphereBSPlineSurface...
    (BSplinePatches,dHat,connections,propShpere,int,graph,'outputEnabled');

%% Plot the reference configuration for the multipatch geometry after the form-finding analysis
% color = [217 218 219]/255;
color = 'none';
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = -310; % -40
el = 40;
axis off;
view(az,el);
camlight(0,0);
title('');
% print('semisphereSetting','-dpng','-r300');

%% END