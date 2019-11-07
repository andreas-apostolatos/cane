%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Modal analysis for a form-found semisphere, whose geometry is
%        directly read from a file.
%
% Date : 13.11.2016
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

%% Read the geometry from the results of a form finding analysis

% Read the data
% FoFiGeo = importdata('./data_FoFiPool/data_FoFiDDMSemisphere.mat');
FoFiGeo = importdata('./data_FoFiPool/data_FoFiDDMExactSemisphere.mat');

% Patch 1 :
% _________

% Polynomial orders
p1 = FoFiGeo.BSplinePatches{1}.p;
q1 = FoFiGeo.BSplinePatches{1}.q;

% Knot vectors
Xi1 = FoFiGeo.BSplinePatches{1}.Xi;
Eta1 = FoFiGeo.BSplinePatches{1}.Eta;

% Control Point coordinates and weights
CP1 = FoFiGeo.BSplinePatches{1}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS1 = FoFiGeo.BSplinePatches{1}.isNURBS;

% Patch 2 :
% _________

% Polynomial orders
p2 = FoFiGeo.BSplinePatches{2}.p;
q2 = FoFiGeo.BSplinePatches{2}.q;

% Knot vectors
Xi2 = FoFiGeo.BSplinePatches{2}.Xi;
Eta2 = FoFiGeo.BSplinePatches{2}.Eta;

% Control Point coordinates and weights
CP2 = FoFiGeo.BSplinePatches{2}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS2 = FoFiGeo.BSplinePatches{2}.isNURBS;

% Patch 3 :
% _________

% Polynomial orders
p3 = FoFiGeo.BSplinePatches{3}.p;
q3 = FoFiGeo.BSplinePatches{3}.q;

% Knot vectors
Xi3 = FoFiGeo.BSplinePatches{3}.Xi;
Eta3 = FoFiGeo.BSplinePatches{3}.Eta;

% Control Point coordinates and weights
CP3 = FoFiGeo.BSplinePatches{3}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS3 = FoFiGeo.BSplinePatches{3}.isNURBS;

% Patch 4 :
% _________

% Polynomial orders
p4 = FoFiGeo.BSplinePatches{4}.p;
q4 = FoFiGeo.BSplinePatches{4}.q;

% Knot vectors
Xi4 = FoFiGeo.BSplinePatches{4}.Xi;
Eta4 = FoFiGeo.BSplinePatches{4}.Eta;

% Control Point coordinates and weights
CP4 = FoFiGeo.BSplinePatches{4}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS4 = FoFiGeo.BSplinePatches{4}.isNURBS;

%% Material constants

% Parameters
EYoung = 7.0e5;
poissonRatio = .49;
thickness = 1.65e-4;
density = 1050;
% prestress.computeParametricCoordinates = @(X) [X(1,1)^2 + X(2,1)^2
%                                                atan(X(2,1)/X(1,1))];
% prestress.computeBaseVectors = @(theta1,theta2) [cos(theta2) -sin(theta2)
%                                                  sin(theta2) cos(theta2)
%                                                  0           0];
prestress.voigtVector = [7794.5
                         7794.5
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
graph.postprocConfig = 'current';

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

a = 1; % || Coarse: 1 || || Fine: 2 ||
tp1 = a;
tq1 = a;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% Patch 2 :
% _________

b = 0; % || Coarse: 0 || || Fine: 1 ||
tp2 = b;
tq2 = b;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');

% Patch 3 :
% _________

c = 1; % || Coarse: 1 || || Fine: 2 ||
tp3 = c;
tq3 = c;
[Xi3,Eta3,CP3,p3,q3] = degreeElevateBSplineSurface(p3,q3,Xi3,Eta3,CP3,tp3,tq3,'');

% Patch 4 :
% _________

d = 0; % || Coarse: 0 || || Fine: 1 ||
tp4 = d;
tq4 = d;
[Xi4,Eta4,CP4,p4,q4] = degreeElevateBSplineSurface(p4,q4,Xi4,Eta4,CP4,tp4,tq4,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Refinement parameters
noHRefStart = 8; % 8
iHRef = 6; % || Coarse: 6 || Fine: 12 ||

% Patch 1 :
% _________

noKnotsXi1 = ceil(5/9*(iHRef + noHRefStart - 1));
noKnotsEta1 = noKnotsXi1;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface...
    (p1,Xi1,q1,Eta1,CP1,noKnotsXi1,noKnotsEta1,'');

% Patch 2 :
% _________

noKnotsXi2 = ceil(6/9*(iHRef + noHRefStart - 1));
noKnotsEta2 = noKnotsXi2;
[Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface...
    (p2,Xi2,q2,Eta2,CP2,noKnotsXi2,noKnotsEta2,'');

% Patch 3 :
% _________

noKnotsXi3 = ceil(4/9*(iHRef + noHRefStart - 1));
noKnotsEta3 = noKnotsXi3;
[Xi3,Eta3,CP3] = knotRefineUniformlyBSplineSurface...
    (p3,Xi3,q3,Eta3,CP3,noKnotsXi3,noKnotsEta3,'');

% Patch 4 :
% _________

noKnotsXi4 = ceil(7/9*(iHRef + noHRefStart - 1));
noKnotsEta4 = noKnotsXi4;
[Xi4,Eta4,CP4] = knotRefineUniformlyBSplineSurface...
    (p4,Xi4,q4,Eta4,CP4,noKnotsXi4,noKnotsEta4,'');

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous DOFs
homDOFs1 = [];
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

% Parameters for the internal pressure
scaling = 1.0;
loadAmplitude1 = scaling*-43.2;
dirForce1 = 'normal';
isFollower1 = true;

% Parameters for weight load
scaling = 1.0;
gravitationalAcceleration = - 9.81;
loadAmplitude2 = scaling*gravitationalAcceleration*density*thickness; % -1.6995825
dirForce2 = 'z';
isFollower2 = false;

% Patch 1 :
% _________

NBC1.noCnd = 2;
xib1 = [0 1];   etab1 = [0 1];
NBC1.xiLoadExtension = {xib1 xib1};
NBC1.etaLoadExtension = {etab1 etab1};
NBC1.loadAmplitude = {-loadAmplitude1 loadAmplitude2};
NBC1.loadDirection = {dirForce1 dirForce2};
NBC1.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
NBC1.isFollower = [isFollower1; isFollower2];
NBC1.isTimeDependent = [false; false];
        
% Patch 2 :
% _________

NBC2.noCnd = 2;
xib2 = [0 1];   etab2 = [0 1];
NBC2.xiLoadExtension = {xib2 xib2};
NBC2.etaLoadExtension = {etab2 etab2};
NBC2.loadAmplitude = {loadAmplitude1 loadAmplitude2};
NBC2.loadDirection = {dirForce1 dirForce2};
NBC2.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
NBC2.isFollower = [isFollower1; isFollower2];
NBC2.isTimeDependent = [false; false];

% Patch 3 :
% _________

NBC3.noCnd = 2;
xib3 = [0 1];   etab3 = [0 1];
NBC3.xiLoadExtension = {xib3 xib3};
NBC3.etaLoadExtension = {etab3 etab3};
NBC3.loadAmplitude = {-loadAmplitude1 loadAmplitude2};
NBC3.loadDirection = {dirForce1 dirForce2};
NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
NBC3.isFollower = [isFollower1; isFollower2];
NBC3.isTimeDependent = [false; false];

% Patch 4 :
% _________

NBC4.noCnd = 2;
xib4 = [0 1];   etab4 = [0 1];
NBC4.xiLoadExtension = {xib4 xib4};
NBC4.etaLoadExtension = {etab4 etab4};
NBC4.loadAmplitude = {loadAmplitude1 loadAmplitude2};
NBC4.loadDirection = {dirForce1 dirForce2};
NBC4.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
NBC4.isFollower = [isFollower1; isFollower2];
NBC4.isTimeDependent = [false; false];

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
connectionsLM = connections;
connectionsLM.lambda = struct([]);
connectionsALM = connections;
connectionsALM.lambda = struct([]);

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

%% Create Lagrange Multipiers fields for all interfaces
fprintf('**************************************************\n');
fprintf('* Creating interface Lagrange Multipliers fields *\n');
fprintf('**************************************************\n\n');
for iConnections = 1:connections.No
    %% Get the IDs of the patches involved
    idI = connections.xiEtaCoup(iConnections,1);
    idJ = connections.xiEtaCoup(iConnections,2);
    fprintf('Coupling between patches %d and %d \n',idI,idJ);
    fprintf('---------------------------------- \n\n');
    
    %% Create a basic Lagrange Multipliers field
    pLambda = 0;
    XiLambda = [0 1];
    CPLambda(:,4) = [1];
    isNURBSLambda = 0;
    nxiLambda = length(CPLambda(:,1,1));
    for i = 1:nxiLambda
        if CPLambda(i,4)~=1
            isNURBSLambda = 1;
            break;
        end
    end
    
    %% Find the interface parametrization for the involved patches
    xiCoupI = connections.xiEtaCoup(iConnections,3:4);
    etaCoupI = connections.xiEtaCoup(iConnections,5:6);
    if xiCoupI(1,1) ~= xiCoupI(1,2) && etaCoupI(1,1) == etaCoupI(1,2)
        isOnXiI = true;
        fprintf('The coupling interface of patch %d is along xi\n',idI);
    elseif xiCoupI(1,1) == xiCoupI(1,2) && etaCoupI(1,1) ~= etaCoupI(1,2)
        isOnXiI = false;
        fprintf('The coupling interface of patch %d is along eta\n',idI);
    else
        error('Either the xi or the eta direction has to be fixed along the interface for patch %d',idI);
    end
    xiCoupJ = connections.xiEtaCoup(iConnections,7:8);
    etaCoupJ = connections.xiEtaCoup(iConnections,9:10);
    if xiCoupJ(1,1) ~= xiCoupJ(1,2) && etaCoupJ(1,1) == etaCoupJ(1,2)
        isOnXiJ = true;
        fprintf('The coupling interface of patch %d is along xi\n',idJ);
    elseif xiCoupJ(1,1) == xiCoupJ(1,2) && etaCoupJ(1,1) ~= etaCoupJ(1,2)
        isOnXiJ = false;
        fprintf('The coupling interface of patch %d is along eta\n',idJ);
    else
        error('Either the xi or the eta direction has to be fixed along the interface for patch %d',idJ);
    end

    %% Degree elevate the Lagrange Multipliers field
    if isOnXiI
        polOrderI = BSplinePatches{idI}.p;
    else
        polOrderI = BSplinePatches{idI}.q;
    end
    if isOnXiJ
        polOrderJ = BSplinePatches{idJ}.p;
    else
        polOrderJ = BSplinePatches{idJ}.q;
    end
    pLM = min(polOrderI,polOrderJ);
    fprintf('Degree elevating the Lagrange Multipliers field to %d\n',pLM);
    if pLM > 0
        clear pLambda XiLambda CPLambda;
        pLambda = 1;
        XiLambda = [0 0 1 1];
        CPLambda(:,4) = [1 1];
        isNURBSLambda = 0;
        nxiLambda = length(CPLambda(:,1,1));
        for i = 1:nxiLambda
            if CPLambda(i,4) ~= 1
                isNURBSLambda = 1;
                break;
            end
        end
        
        % Perform accordingly a p-refinement
        tpLambda = pLM;
        [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
            (pLambda,XiLambda,CPLambda,tpLambda,'');
    end
    
    %% Perform a knot insertion to the Lagrange Multipliers field
    if isOnXiI
        noKnotsI = length(unique(BSplinePatches{idI}.Xi)) - 2;
    else
        noKnotsI = length(unique(BSplinePatches{idI}.Eta)) - 2;
    end
    if isOnXiJ
        noKnotsJ = length(unique(BSplinePatches{idJ}.Xi)) - 2;
    else
        noKnotsJ = length(unique(BSplinePatches{idJ}.Eta)) - 2;
    end
    scaleLM = 1.0;
    noLambda = ceil(min([noKnotsI noKnotsJ])*scaleLM);
    fprintf('Uniformly inserting %d knots in the Lagrange Multipliers field\n',noLambda);
    [XiLambdaLM,CPLambdaLM] = knotRefineUniformlyBSplineCurve...
        (noLambda,pLambda,XiLambda,CPLambda,'');
    scaleALM = 0.5;
    noLambda = ceil(min([noKnotsI noKnotsJ])*scaleALM);
    fprintf('Uniformly inserting %d knots in the augmented Lagrange Multipliers field\n',noLambda);
    [XiLambdaALM,CPLambdaALM] = knotRefineUniformlyBSplineCurve...
        (noLambda,pLambda,XiLambda,CPLambda,'');
    
    %% Fill up the Lagrange Multipliers patch and add it to the array
    lambdaLM = fillUpLagrangeMultipliers...
        (pLambda,XiLambdaLM,CPLambdaLM,isNURBSLambda);
    connectionsLM.lambda{iConnections} = lambdaLM;
    lambdaALM = fillUpLagrangeMultipliers...
        (pLambda,XiLambdaALM,CPLambdaALM,isNURBSLambda);
    connectionsALM.lambda{iConnections} = lambdaALM;
    fprintf('\n');
end

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-4;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 35;

%% Coupling properties

% % Initialize the rotational coupling parameters
% propCoupling.alphaR = zeros(4,1);
% 
% % General parameters
% penaltyPrmScale = 0;
% 
% for iConnections = 1:connections.No
%     % Get the id's of the patches
%     IDPatchI = connections.xiEtaCoup(iConnections,1);
%     IDPatchJ = connections.xiEtaCoup(iConnections,2);
%     
%     % Get the mean polynomial order between the patches
%     isOnXiI = false;
%     if BSplinePatches{IDPatchI}.weakDBC.etaExtension{1}(1,1) == ...
%             BSplinePatches{IDPatchI}.weakDBC.etaExtension{1}(1,2)
%         isOnXiI = true;
%     end
%     if isOnXiI
%         polOrderI = BSplinePatches{IDPatchI}.p;
%     else
%         polOrderI = BSplinePatches{IDPatchI}.q;
%     end
%     isOnXiJ = false;
%     if BSplinePatches{IDPatchJ}.weakDBC.etaExtension{1}(1,1) == ...
%             BSplinePatches{IDPatchJ}.weakDBC.etaExtension{1}(1,2)
%         isOnXiJ = true;
%     end
%     if isOnXiJ
%         polOrderJ = BSplinePatches{IDPatchJ}.p;
%     else
%         polOrderJ = BSplinePatches{IDPatchJ}.q;
%     end
%     polOrderMean = mean([polOrderI polOrderJ]);
%     
%     % Assign the penalty parameters
% %     propCoupling.alphaR(iConnections,1) = max([norm(eval(['Db' num2str(IDPatchI)])) ...
% %         norm(eval(['Db' num2str(IDPatchJ)]))])*polOrderMean/...
% %         min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea])*...
% %         penaltyPrmScale;
% end

%% Choose a method for the application of weak Dirichlet boundary conditions and the multipatch coupling
method = 'Nitsche';
if ~strcmp(method,'Penalty') && ~strcmp(method,'LagrangeMultipliers') && ...
        ~strcmp(method,'AugmentedLagrangeMultipliers') && ~strcmp(method,'Nitsche')
    error('%s is not a valid method (Nitsche, Penalty, LagrangeMultipliers, AugmentedLagrangeMultipliers)',method);
end

%% Set up the parameters and properties for each method
if strcmp(method,'Penalty')
    % General parameters
    penaltyPrmScale = 1e0;

    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    for iPatches = 1:noPatches
        % Assign the method name
        BSplinePatches{iPatches}.weakDBC.method = 'penalty';

        % Get the polynomial order along the Dirichlet boundary
        isOnXi = false;
        if BSplinePatches{iPatches}.weakDBC.etaExtension{1}(1,1) == ...
                BSplinePatches{iPatches}.weakDBC.etaExtension{1}(1,2)
            isOnXi = true;
        end
        if isOnXi
            polOrder = BSplinePatches{iPatches}.p;
        else
            polOrder = BSplinePatches{iPatches}.q;
        end

        % Assign the penalty parameter
        BSplinePatches{iPatches}.weakDBC.alpha = ...
            norm(eval(['Dm' num2str(iPatches)]))*polOrder/...
            sqrt(BSplinePatches{iPatches}.minElArea)*...
            penaltyPrmScale;
    end

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'penalty';
    propCoupling.intC = intC;
    propCoupling.alphaD = zeros(connections.No,1);
    propCoupling.alphaR = zeros(connections.No,1);
    for iConnections = 1:connections.No
        % Get the id's of the patches
        IDPatchI = connections.xiEtaCoup(iConnections,1);
        IDPatchJ = connections.xiEtaCoup(iConnections,2);

        % Get the mean polynomial order between the patches
        isOnXiI = false;
        if BSplinePatches{IDPatchI}.weakDBC.etaExtension{1}(1,1) == ...
                BSplinePatches{IDPatchI}.weakDBC.etaExtension{1}(1,2)
            isOnXiI = true;
        end
        if isOnXiI
            polOrderI = BSplinePatches{IDPatchI}.p;
        else
            polOrderI = BSplinePatches{IDPatchI}.q;
        end
        isOnXiJ = false;
        if BSplinePatches{IDPatchJ}.weakDBC.etaExtension{1}(1,1) == ...
                BSplinePatches{IDPatchJ}.weakDBC.etaExtension{1}(1,2)
            isOnXiJ = true;
        end
        if isOnXiJ
            polOrderJ = BSplinePatches{IDPatchJ}.p;
        else
            polOrderJ = BSplinePatches{IDPatchJ}.q;
        end
        polOrderMean = mean([polOrderI polOrderJ]);

        % Assign the penalty parameters
        propCoupling.alphaD(iConnections,1) = ...
            max([norm(eval(['Dm' num2str(IDPatchI)])) ...
            norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
            sqrt(min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea]))*...
            penaltyPrmScale;
    end
elseif strcmp(method,'LagrangeMultipliers')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    % Properties for the weak Dirichlet boundary conditions
    for iPatches = 1:noPatches
        BSplinePatches{iPatches}.weakDBC.method = 'lagrangeMultipliers';
        BSplinePatches{iPatches}.weakDBC.alpha = 0;

        % Find along which parametric line the weak Dirichlet 
        % condition is to be imposed
        isOnXi = false;
        if BSplinePatches{iPatches}.weakDBC.xiExtension{1}(1,1) == ...
             BSplinePatches{iPatches}.weakDBC.xiExtension{1}(1,2)
            isOnXi = true;
        end

        % Make a Lagrange Multipliers discretization
        clear pLambda XiLambda CPLambda isNURBSLambda; 

        % Find out up to which polynomial degree the Lagrange
        % Multipliers discretization needs to be increased
        if isOnXi
            polOrderPatch =  BSplinePatches{iPatches}.p;
        else
            polOrderPatch =  BSplinePatches{iPatches}.q;
        end
        pLM = polOrderPatch - 2;

        if pLM <= 0
            pLambda = 0;
            XiLambda = [0 1];
            CPLambda = zeros(1,4);
        else
            pLambda = 1;
            XiLambda = [0 0 1 1];
            CPLambda = zeros(2,4);

            % Peform a p-refinement
            tpLambda = polOrderPatch - pLM;
            [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
                (pLambda,XiLambda,CPLambda,tpLambda,'');
        end
        isNURBSLambda = 0;
        nxiLambda = length(CPLambda(:,1));
        for i = 1:nxiLambda
            if CPLambda(i,4) ~= 1
                isNURBSLambda = 1;
                break;
            end
        end

        % Perform an h-refinement
        percentage = 1.0;
        if isOnXi
            Rxi = unique(BSplinePatches{iPatches}.Xi);
        else
            Rxi = unique(BSplinePatches{iPatches}.Eta);
        end
        noXi = ceil(percentage*(length(Rxi) - 2));
        [XiLambda,CPLambda] = knotRefineUniformlyBSplineCurve...
            (noXi,pLambda,XiLambda,CPLambda,'');

        % Create the corresponding Lagrange Multipliers structure
        BSplinePatches{iPatches}.weakDBC.lambda{1} = fillUpLagrangeMultipliers...
         (pLambda,XiLambda,CPLambda,isNURBSLambda);
    end

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'lagrangeMultipliers';
    connections = connectionsLM;
    propCoupling.alphaD = zeros(connections.No,1);
    propCoupling.alphaR = zeros(connections.No,1);
    propCoupling.intC = intC;
elseif strcmp(method,'AugmentedLagrangeMultipliers')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------
    
    for iPatches = 1:noPatches
        % Assign the method name
        BSplinePatches{iPatches}.weakDBC.method = 'lagrangeMultipliers';

        % Find along which parametric line the weak Dirichlet
        % condition is to be imposed
        isOnXi = false;
        if BSplinePatches{iPatches}.weakDBC.xiExtension{1}(1,1) == ...
             BSplinePatches{iPatches}.weakDBC.xiExtension{1}(1,2)
            isOnXi = true;
        end

        % Find out up to which polynomial degree the Lagrange
        % Multipliers discretization needs to be increased
        if isOnXi
            polOrderPatch = BSplinePatches{iPatches}.p;
        else
            polOrderPatch = BSplinePatches{iPatches}.q;
        end
        pLM = polOrderPatch - 2; %polOrderPatch - 2 has worked well

        % Assign the penalty parameter
        scalePenaltyFactorALM = 1e-2; 
        BSplinePatches{iPatches}.weakDBC.alpha = ...
            norm(eval(['Dm' num2str(iPatches)]))*polOrderPatch/...
            sqrt(BSplinePatches{iPatches}.minElArea)*scalePenaltyFactorALM;

        % Make a Lagrange Multipliers discretization
        clear pLambda XiLambda CPLambda isNURBSLambda; 
        if pLM <= 0
            pLambda = 0;
            XiLambda = [0 1];
            CPLambda = zeros(1,4);
        else
            pLambda = 1;
            XiLambda = [0 0 1 1];
            CPLambda = zeros(2,4);

            % Peform a p-refinement
            tpLambda = polOrderPatch - pLM;
            [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
                (pLambda,XiLambda,CPLambda,tpLambda,'');
        end
        isNURBSLambda = 0;
        nxiLambda = length(CPLambda(:,1));
        for i = 1:nxiLambda
            if CPLambda(i,4) ~= 1
                isNURBSLambda = 1;
                break;
            end
        end

        % Perform an h-refinement
        percentage = .5;
        if isOnXi
            Rxi = unique(BSplinePatches{iPatches}.Xi);
        else
            Rxi = unique(BSplinePatches{iPatches}.Eta);
        end
        noXi = ceil(percentage*(length(Rxi) - 2));
        [XiLambda,CPLambda] = knotRefineUniformlyBSplineCurve...
            (noXi,pLambda,XiLambda,CPLambda,'');

        % Create the corresponding Lagrange Multipliers structure
        BSplinePatches{iPatches}.weakDBC.lambda{1} = fillUpLagrangeMultipliers...
            (pLambda,XiLambda,CPLambda,isNURBSLambda);
    end
    
    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'lagrangeMultipliers';
    connections = connectionsALM;
    propCoupling.intC = intC;
    propCoupling.alphaD = zeros(connections.No,1);
    propCoupling.alphaR = zeros(connections.No,1);
    for iConnections = 1:connections.No
        % Get the Patch IDs
        IDPatchI = connections.xiEtaCoup(iConnections,1);
        IDPatchJ = connections.xiEtaCoup(iConnections,2);

        % Get the mean polynomial order between the patches
        isOnXiI = false;
        if BSplinePatches{IDPatchI}.weakDBC.etaExtension{1}(1,1) == ...
                BSplinePatches{IDPatchI}.weakDBC.etaExtension{1}(1,2)
            isOnXiI = true;
        end
        if isOnXiI
            polOrderI = BSplinePatches{IDPatchI}.p;
        else
            polOrderI = BSplinePatches{IDPatchI}.q;
        end
        isOnXiJ = false;
        if BSplinePatches{IDPatchJ}.weakDBC.etaExtension{1}(1,1) == ...
                BSplinePatches{IDPatchJ}.weakDBC.etaExtension{1}(1,2)
            isOnXiJ = true;
        end
        if isOnXiJ
            polOrderJ = BSplinePatches{IDPatchJ}.p;
        else
            polOrderJ = BSplinePatches{IDPatchJ}.q;
        end
        polOrderMean = mean([polOrderI polOrderJ]);

        % Assign the penalty parameter
        propCoupling.alphaD(iConnections,1) = ...
            max([norm(eval(['Dm' num2str(IDPatchI)])) ...
            norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
            sqrt(min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea]))*...
            scalePenaltyFactorALM;
    end
elseif strcmp(method,'Nitsche')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    % Assing the B-Spline patches for the Penalty method
    if exist('BSplinePatchesNitsche','var');
        clear BSplinePatchesNitsche;
    end

    % Properties for the weak Dirichlet boundary conditions
    for iPatches = 1:noPatches
        BSplinePatches{iPatches}.weakDBC.method = 'nitsche';
        BSplinePatches{iPatches}.weakDBC.estimationStabilPrm = true;
    end

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'nitsche';
    propCoupling.estimationStabilPrm = true;
    propCoupling.gammaTilde = .5;
    propCoupling.intC = intC;
end

%% Perform modal analysis and visualize the chosen eigenmode shape

% Solve the eigenvalue problem
noEig = 400; % 40
[eigenmodeShapes,naturalFrequencies,BSplinePatches,dHat] = ...
    solve_DDMEigenmodeAnalysisIGAMembrane...
    (BSplinePatches,connections,propCoupling,solve_LinearSystem,noEig,...
    propNLinearAnalysis,'outputEnabled');
save(['data_ModalAnalysis' method 'Coarse']);
return;

%% Load a reference solution for the modal analysis
modalAnalysisSP = importdata('./data_ModalAnalysisPool/data_modalAnalysisSemisphereElmnts1024P4Q4FF1e-11.mat');
naturalFrequenciesSP = modalAnalysisSP.naturalFrequencies;

%% Plot the error in the numerical spectra
omegaSP = (2*pi).*naturalFrequenciesSP;
omegaMP = (2*pi).*naturalFrequencies;
omegaRatio = zeros(noEig,1);
for iFreq = 1:noEig
    omegaRatio(iFreq,1) = omegaMP(iFreq,1)/omegaSP(iFreq,1);
end
xAxis = 1:noEig;
plot(xAxis,omegaRatio);
return;

%% Visualize selected mode shapes

% Find the numbering of the DOFs which are purely associated with
% displacements (that is no Lagrange Multiplier DOFs)
noDOFsWTLM = 0;
for iPatches = 1:noPatches
    noDOFsWTLMPatch = 3*BSplinePatches{iPatches}.noCPs;
    noDOFsWTLM = noDOFsWTLM + noDOFsWTLMPatch;
end
EFTPatchesWTLM = zeros(1,noDOFsWTLM);
noDOFsSaved = 0;
for iPatches = 1:noPatches
    noDOFsWTLMPatch = 3*BSplinePatches{iPatches}.noCPs;
    EFTPatchesWTLM(1,noDOFsSaved + 1:noDOFsSaved + noDOFsWTLMPatch) = ...
        BSplinePatches{iPatches}.EFTPatches(1,1:noDOFsWTLMPatch);
    noDOFsSaved = noDOFsSaved + noDOFsWTLMPatch;
end

% Create a new B-Spline array with displaced Control Points according to
% the static equilibrium configuration computed at the static step
BSplinePatchesEq = BSplinePatches;
for iPatches = 1:noPatches
    BSplinePatchesEq{iPatches}.CP = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
    (BSplinePatches{iPatches}.CP,dHat(EFTPatchesWTLM));
end

% Eigenfrequencies used for estimating the damping coefficients
freqI = 41.99;
freqJ = 72.27;

% Find the position of the eigenfrequencies in the frequency array
[~,idI] = min(abs(naturalFrequencies - freqI));
[~,idJ] = min(abs(naturalFrequencies - freqJ));

% Visualization of the eigenmode shapes
for i = idJ % [2 idI idJ]
    scalingFct = .001; % -.001;
    noEigen = i;
    graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
        (BSplinePatchesEq,scalingFct*eigenmodeShapes(EFTPatchesWTLM,noEigen),graph,'outputEnabled');
    az = -310; % -310
    el = 40; % 40
    view(az,el);
    camlight(0,0);
    lighting phong;
    axis off;
    title('');
end

%% END
