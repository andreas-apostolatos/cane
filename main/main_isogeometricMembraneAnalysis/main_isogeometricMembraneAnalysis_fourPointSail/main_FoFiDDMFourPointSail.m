%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Form-finding analysis over a 3-patch four point sail
%
% Date : 10.01.2017
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
        '../../../isogeometricThinStructureAnalysis/penaltyDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../../isogeometricThinStructureAnalysis/output/',...
        '../../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/');

%% CAD modelling via NURBS

% Global variables:
Length = 20;
Width = Length;
Height = Length/2;

% Patch 1 :
% _________

% Polynomial degrees
p1 = 1;
q1 = 2;

% Knot vectors
Xi1 = [0 0 1 1];
Eta1 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP1(:,:,1) = [-Length/2 -Length/2 -Length/2
              -Length/2 -Length/3 Length/2];
         
% y-coordinates
CP1(:,:,2) = [-Width/2 -Width/2 -Width/2
              Width/2  -Width/3 -Width/2];
         
% z-coordinates
CP1(:,:,3) = [0      0                 0
              Height abs(CP1(2,2,1))/2 Height];
       
% Weights
CP1(:,:,4) = [1 1 1
              1 1 1];

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

% Patch 3 :
% _________

% Polynomial degrees
p3 = 1;
q3 = 2;

% Knot vectors
Xi3 = [0 0 1 1];
Eta3 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP3(:,:,1) = [-Length/2 Length/3 Length/2
              Length/2  Length/2 Length/2];
         
% y-coordinates
CP3(:,:,2) = [Width/2 Width/3 -Width/2
              Width/2 Width/2 Width/2];
         
% z-coordinates
CP3(:,:,3) = [Height abs(CP3(1,2,1))/2 Height
              0      0                 0];
       
% Weights
CP3(:,:,4) = [1 1 1
              1 1 1];

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

% Patch 2 :
% _________

% Polynomial degrees
p2 = 1;
q2 = 2;

% Knot vectors
Xi2 = [0 0 1 1];
Eta2 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP2(:,:,1) = [CP1(2,1,1) CP1(2,2,1) CP1(2,3,1)
              CP3(1,1,1) CP3(1,2,1) CP3(1,3,1)];
         
% y-coordinates
CP2(:,:,2) = [CP1(2,1,2) CP1(2,2,2) CP1(2,3,2)
              CP3(1,1,2) CP3(1,2,2) CP3(1,3,2)];
         
% z-coordinates
CP2(:,:,3) = [CP1(2,1,3) CP1(2,2,3) CP1(2,3,3)
              CP3(1,1,3) CP3(1,2,3) CP3(1,3,3)];
       
% Weights
CP2(:,:,4) = [CP1(2,1,4) CP1(2,2,4) CP1(2,3,4)
              CP3(1,1,4) CP3(1,2,4) CP3(1,3,4)];

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

%% Material constants

% general parameters
EYoung = 8e+8;
nue = .4;
thickness = 1e-3;
sigma0 = 3e+3;
prestress.voigtVector = [sigma0/thickness
                         sigma0/thickness
                         0];
density = 8050;

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
parameters3.t = thickness;

% Density of the membrane
parameters3.rho = density;

% Prestress of the membrane
parameters3.prestress = prestress;

% Cable :
% _______

parametersCable.E = 1.6e+11;
parametersCable.radiusCS = 12e-3/2;
parametersCable.areaCS = pi*parametersCable.radiusCS^2;
parametersCable.rho = 8050;
parametersCable.prestress = 6e+4/parametersCable.areaCS;

%% GUI

% Case name
caseName = 'DDM3PatchesFourPointSail';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points

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

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','curvature','force','moment','shearForce'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x','y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

% Define the coupling properties

% Patch 1 :
% _________

Dm1 = parameters1.E*parameters1.t/(1-parameters1.nue^2)*...
      [1              parameters1.nue 0
       parameters1.nue 1              0
       0               0              (1-parameters1.nue)/2];
   
% Patch 2 :
% _________

Dm2 = parameters2.E*parameters2.t/(1-parameters2.nue^2)*...
      [1              parameters2.nue 0
       parameters2.nue 1              0
       0               0              (1-parameters2.nue)/2];

% Patch 3 :
% _________

Dm3 = parameters3.E*parameters3.t/(1-parameters3.nue^2)*...
      [1              parameters3.nue 0
       parameters3.nue 1              0
       0               0              (1-parameters3.nue)/2];

% Assign the penalty factors

% Function handle to writing out the results
% writeOutput = @writeResults4Carat;
% writeOutput = @writeResults4GiD;
writeOutput = 'undefined';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

%% Refinement

% Select p-refinement level
iPRef = 1; % Coarse : 1 || Fine : 2

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = iPRef;
tp1 = a + 1;
tq1 = a;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% Patch 2 :
% _________

b = iPRef + 1;
tp2 = b + 1;
tq2 = b;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');

% Patch 3 :
% _________

c = iPRef;
tp3 = c + 1;
tq3 = c;
[Xi3,Eta3,CP3,p3,q3] = degreeElevateBSplineSurface(p3,q3,Xi3,Eta3,CP3,tp3,tq3,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% General parameter used for preserving the produced elements 
% shaped as close as possible to a rectangle
edgeRatio = ceil((2/3)*Length/Width/(sin(4*pi/9))/2);

% Select h-refinement level
iHRef = 6; % Coarse : 6 || Fine : 12

% Patch 1 :
% _________

noKnotsXi1 = ceil((14/9)*(iHRef + 1));
noKnotsEta1 = edgeRatio*noKnotsXi1;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface...
    (p1,Xi1,q1,Eta1,CP1,noKnotsXi1,noKnotsEta1,'');

% Patch 2 :
% _________

noKnotsXi2 = ceil((iHRef + 1)*(10/15));
noKnotsEta2 = ceil(4*edgeRatio*noKnotsXi2*(2/3));
[Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface...
    (p2,Xi2,q2,Eta2,CP2,noKnotsXi2,noKnotsEta2,'');

% Patch 3 :
% _________

noKnotsXi3 = ceil((9/6)*(iHRef + 1));
noKnotsEta3 = edgeRatio*noKnotsXi3;
[Xi3,Eta3,CP3] = knotRefineUniformlyBSplineSurface...
    (p3,Xi3,q3,Eta3,CP3,noKnotsXi3,noKnotsEta3,'');

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs1 = [];
xisup1 = [0 0];   etasup1 = [0 1];
for dir = 1:3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end
xisup1 = [1 1];   etasup1 = [0 0];
for dir = 1:3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end
xisup1 = [1 1];   etasup1 = [1 1];
for dir = 1:3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC1.noCnd = 0;

% Embedded cables
cables1.No = 2;
cables1.xiExtension = {[0 1] [0 1]};
cables1.etaExtension = {[0 0] [1 1]};
cables1.parameters = {parametersCable parametersCable parametersCable parametersCable};
cables1.int.type = 'default';
% cables1.int.type = 'user';
cables1.int.noGPs = 16;

% Patch 2 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs2 = [];
xisup2 = [0 1];   etasup2 = [0 0];
for dir = 1:3
    homDOFs2 = findDofs3D...
        (homDOFs2,xisup2,etasup2,dir,CP2);
end
xisup2 = [0 1];   etasup2 = [1 1];
for dir = 1:3
    homDOFs2 = findDofs3D...
        (homDOFs2,xisup2,etasup2,dir,CP2);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs2 = [];
valuesInhomDOFs2 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC2.noCnd = 0;

% Embedded cables
cables2.No = 0;

% Patch 3 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs3 = [];
xisup3 = [0 0];   etasup3 = [0 0];
for dir = 1:3
    homDOFs3 = findDofs3D...
        (homDOFs3,xisup3,etasup3,dir,CP3);
end
xisup3 = [0 0];   etasup3 = [1 1];
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

% Weak homogeneous Dirichlet boundary conditions
weakDBC3.noCnd = 0;

% Embedded cables
cables3.No = 2;
cables3.xiExtension = {[0 1] [0 1]};
cables3.etaExtension = {[0 0] [1 1]};
cables3.parameters = {parametersCable parametersCable parametersCable parametersCable};
cables3.int.type = 'default';
% cables3.int.type = 'user';
cables3.int.noGPs = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parameter
loadAmplitude = 0;

% Patch 1 :
% _________

FAmp1 = loadAmplitude;
NBC1.noCnd = 1;
xib1 = [0 1];   etab1 = [0 1];   dirForce1 = 'z';
NBC1.xiLoadExtension = {xib1};
NBC1.etaLoadExtension = {etab1};
NBC1.loadAmplitude = {FAmp1};
NBC1.loadDirection = {dirForce1};
NBC1.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC1.isFollower = false;
NBC1.isTimeDependent = true;

% Patch 2 :
% _________

FAmp2 = - loadAmplitude;
NBC2.noCnd = 1;
xib2 = [0 1];   etab2 = [0 1];   dirForce2 = 'z';
NBC2.xiLoadExtension = {xib2};
NBC2.etaLoadExtension = {etab2};
NBC2.loadAmplitude = {FAmp2};
NBC2.loadDirection = {dirForce2};
NBC2.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC2.isFollower = false;
NBC2.isTimeDependent = true;

% Patch 3 :
% _________

FAmp3 = - loadAmplitude;
NBC3.noCnd = 1;
xib3 = [0 1];   etab3 = [0 1];   dirForce3 = 'z';
NBC3.xiLoadExtension = {xib3};
NBC3.etaLoadExtension = {etab3};
NBC3.loadAmplitude = {FAmp3};
NBC3.loadDirection = {dirForce3};
NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC3.isFollower = false;
NBC3.isTimeDependent = true;

% Collect all the Neumann boundary conditions into an arra<y
NBC = {NBC1 NBC2 NBC3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface parametrizations     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Connection with patch 2:
xicoup12 = [1 1];   etacoup12 = [0 1];

% Collect all interfaces into arrays:
xicoup1 = xicoup12;
etacoup1 = etacoup12;

% Patch 2 :
% _________

% Connection with patch 1:
xicoup21 = [0 0];   etacoup21 = [0 1];

% Connection with patch 3:
xicoup23 = [1 1];   etacoup23 = [0 1];

% Collect all interfaces into arrays:
xicoup2 = [xicoup21 xicoup23];
etacoup2 = [etacoup21 etacoup23];

% Patch 3 :
% _________

% Connection with patch 2:
xicoup32 = [0 0];   etacoup32 = [0 1];

% Collect all interfaces into arrays:
xicoup3 = xicoup32;
etacoup3 = etacoup32;

% Define connections :
% ____________________

% Number of connections
noConnections = 2;

% Define connections by patch numbers
connections.No = noConnections;
connections.xiEtaCoup = zeros(noConnections,10);
connections.xiEtaCoup(:,:) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21
                              2 3 xicoup23 etacoup23 xicoup32 etacoup32];
connectionsLM = connections;
connectionsALM = connections;

%% Fill up the arrays for the patches

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

% 2nd patch :
% ___________

patch3 = fillUpPatch...
    (analysis,p3,Xi3,q3,Eta3,CP3,isNURBS3,parameters3,homDOFs3,...
    inhomDOFs3,valuesInhomDOFs3,weakDBC3,cables3,NBC3,[],[],[],xicoup3,etacoup3,int3);

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch1 patch2 patch3};
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

%% Plot the distribution of the determinant of the geometrical Jacobian for the multipatch geometry
% figure(graph.index)
% for iPatches = 1:noPatches
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
% color = [217 218 219]/255;
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,color,graph);

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
color = [217 218 219]/255;
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
% az = -40;
% el = 40;
% view(az,el);
% axis off;
% title('');
% camlight left;
% lighting phong;

%% Load a Finite Element solution
% nodes = importdata('../../../preComputedData/isogeometricMembraneAnalysis/4ptSail/referenceFEMSolution/nodes');
% elements = importdata('../../../preComputedData/isogeometricMembraneAnalysis/4ptSail/referenceFEMSolution/elements');
% displacement = importdata('../../../preComputedData/isogeometricMembraneAnalysis/4ptSail/referenceFEMSolution/displacements');
% mesh.elements = elements(:,1:3);
% for iNodes = 1:length(nodes(:,1))
%     [index,~] = find(displacement(iNodes,1) == nodes(:,1));
%     nodes(index,2:4) = nodes(index,2:4) + displacement(iNodes,2:4);
% end
% mesh.nodes = nodes;

%% Plot the form-found geometry using classical Finite Elements
% color = [217 218 219]/255;
% labelsEnabled = false;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,mesh,labelsEnabled,color,graph);
% az = 30;
% el = 45;
% view([az el]);
% camlight(20,40);
% % camlight left;
% lighting phong;

%% Perform a form finding analysis to find the shape of the multipatch membrane with the Nitsche method

% % Properties for the form-finding analysis
% % ----------------------------------------
% 
% propFormFinding.tolerance = 1e-9;
% propFormFinding.maxNoIter = 20;
% propFormFinding.minNoIter = 1;
% penaltyPrmScale = 1e6;
% 
% % Assign the parameters for the application of weak DBC
% % -----------------------------------------------------
% 
% % Assing the B-Spline patches for the Penalty method
% if exist('BSplinePatchesNitsche','var');
%     clear BSplinePatchesNitsche;
% end
% BSplinePatchesNitsche = BSplinePatches;
% 
% % Properties for the weak Dirichlet boundary conditions
% for iPatches = 1:noPatches
%     BSplinePatchesNitsche{iPatches}.weakDBC.method = 'nitsche';
%     BSplinePatchesNitsche{iPatches}.weakDBC.estimationStabilPrm = true;
% end
% 
% % Assign the parameters for multipatch coupling
% % ---------------------------------------------
% 
% propCouplingNitsche.method = 'nitsche';
% propCouplingNitsche.estimationStabilPrm = true;
% propCouplingNitsche.gammaTilde = .5;
% propCouplingNitsche.intC = intC;
% 
% % Solve the problem using the Penalty method
% % ------------------------------------------
% 
% [BSplinePatches,CPHistoryMultipatch,propCouplingNitsche,resHistory,...
%     hasConverged,noIter] = ...
%     solve_DDMFormFindingIGAMembrane(BSplinePatchesNitsche,connections,...
%     propCouplingNitsche,propFormFinding,solve_LinearSystem,'outputEnabled');
% 
% % Plot the result
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');
% 
% return;
% save data_FoFiDDMFourPointSailCoarse;

%% Perform a form finding analysis to find the shape of the multipatch membrane with the Penalty method

% Properties for the form-finding analysis
% ----------------------------------------

propFormFinding.tolerance = 1e-3;
propFormFinding.maxNoIter = 50;
propFormFinding.minNoIter = 1;
penaltyPrmScale = 1e0;

% Assign the parameters for the application of weak DBC
% -----------------------------------------------------

% Assing the B-Spline patches for the Penalty method
if exist('BSplinePatchesPenalty','var')
    clear BSplinePatchesPenalty;
end
BSplinePatchesPenalty = BSplinePatches;

% Properties for the weak Dirichlet boundary conditions
for iPatches = 1:noPatches
    if BSplinePatchesPenalty{iPatches}.weakDBC.noCnd > 0
        % Assign the method name
        BSplinePatchesPenalty{iPatches}.weakDBC.method = 'penalty';

        % Get the polynomial order along the Dirichlet boundary
        isOnXi = false;
        if BSplinePatchesPenalty{iPatches}.weakDBC.etaExtension{1}(1,1) == ...
                BSplinePatchesPenalty{iPatches}.weakDBC.etaExtension{1}(1,2)
            isOnXi = true;
        end
        if isOnXi
            polOrder = BSplinePatchesPenalty{iPatches}.p;
        else
            polOrder = BSplinePatchesPenalty{iPatches}.q;
        end

        % Assign the penalty parameter
        BSplinePatchesPenalty{iPatches}.weakDBC.alpha = ...
            norm(eval(['Dm' num2str(iPatches)]))*polOrder/...
            BSplinePatches{iPatches}.minElArea*...
            penaltyPrmScale;
    end
end

% Assign the parameters for multipatch coupling
% ---------------------------------------------

propCouplingPenalty.method = 'penalty';
propCouplingPenalty.intC = intC;
propCouplingPenalty.alphaD = zeros(connections.No,1);
propCouplingPenalty.alphaR = zeros(connections.No,1);
for iConnections = 1:connections.No
    % Get the id's of the patches
    IDPatchI = connections.xiEtaCoup(iConnections,1);
    IDPatchJ = connections.xiEtaCoup(iConnections,2);

    % Get the mean polynomial order between the patches
    isOnXiI = false;
    if connections.xiEtaCoup(iConnections,5) == connections.xiEtaCoup(iConnections,6)
        isOnXiI = true;
    end
    if isOnXiI
        polOrderI = BSplinePatchesPenalty{IDPatchI}.p;
    else
        polOrderI = BSplinePatchesPenalty{IDPatchI}.q;
    end
    isOnXiJ = false;
    if connections.xiEtaCoup(iConnections,9) == connections.xiEtaCoup(iConnections,10)
        isOnXiJ = true;
    end
    if isOnXiJ
        polOrderJ = BSplinePatchesPenalty{IDPatchJ}.p;
    else
        polOrderJ = BSplinePatchesPenalty{IDPatchJ}.q;
    end
    polOrderMean = mean([polOrderI polOrderJ]);

    % Assign the penalty parameters
    propCouplingPenalty.alphaD(iConnections,1) = ...
        max([norm(eval(['Dm' num2str(IDPatchI)])) ...
        norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
        min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea])*...
        penaltyPrmScale;
    propCouplingPenalty.alphaR(iConnections,1) = 0;
end

% Solve the problem using the Penalty method
% ------------------------------------------

[BSplinePatches,CPHistoryMultipatch,propCoupling,resHistory,...
    hasConverged,noIter] = ...
    solve_DDMFormFindingIGAMembrane(BSplinePatchesPenalty,connections,...
    propCouplingPenalty,propFormFinding,solve_LinearSystem,'outputEnabled');
save data_FoFiDDMFourPointSail;

%% Solve the steady-state problem with the mortar method

% % Properties for the form-finding analysis
% % ----------------------------------------
% 
% propFormFinding.tolerance = 1e-6;
% propFormFinding.maxNoIter = 10;
% propFormFinding.minNoIter = 1;
% 
% % Assign the parameters for the application of weak DBC
% % -----------------------------------------------------
% 
% % Assing the B-Spline patches for the Penalty method
% if exist('BSplinePatchesMortar','var');
%     clear BSplinePatchesMortar;
% end
% BSplinePatchesMortar = BSplinePatches;
% 
% % Properties for the weak Dirichlet boundary conditions
% for iPatches = 1:noPatches
%     BSplinePatchesMortar{iPatches}.weakDBC.method = 'nitsche';
%     BSplinePatchesMortar{iPatches}.weakDBC.estimationStabilPrm = true;
% end
% 
% % Assign the parameters for multipatch coupling
% % ---------------------------------------------
% 
% propCouplingMortar.method = 'mortar';
% propCouplingMortar.isSlaveSideCoarser = false;
% propCouplingMortar.computeRearrangedProblemMtrcs = @computeRearrangedProblemMtrcs4MortarIGAMembrane;
% propCouplingMortar.intC = intC;
% 
% % Solve the problem using the mortar method
% % -------------------------------------------------------
% 
% [BSplinePatches,CPHistoryMultipatch,propCoupling,resHistory,...
%     hasConverged,noIter] = solve_DDMFormFindingIGAMembrane...
%     (BSplinePatchesMortar,connections,propCouplingMortar,propFormFinding,...
%     solve_LinearSystem,'outputEnabled');
% % return;
% % save data_FoFiDDMFourPointSailCoarse;

%% Perform a form finding analysis to find the shape of the multipatch membrane with the Lagrange Multipliers method

% % Scale the penalty parameters
% penaltyPrmScale = 0;%1e-2;
% 
% % Properties for the form-finding analysis
% % ----------------------------------------
% 
% propFormFinding.tolerance = 1e-6;
% propFormFinding.maxNoIter = 1e3;
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
%     if BSplinePatchesLM{iPatches}.weakDBC.noCnd > 0
%         BSplinePatchesLM{iPatches}.weakDBC.method = 'lagrangeMultipliers';
% 
%         % Get the polynomial order along the Dirichlet boundary
%         isOnXi = false;
%         if BSplinePatchesLM{iPatches}.weakDBC.etaExtension{1}(1,1) == ...
%                 BSplinePatchesLM{iPatches}.weakDBC.etaExtension{1}(1,2)
%             isOnXi = true;
%         end
%         if isOnXi
%             polOrder = BSplinePatchesLM{iPatches}.p;
%         else
%             polOrder = BSplinePatchesLM{iPatches}.q;
%         end
% 
%         % Assign the penalty parameter
%         BSplinePatchesLM{iPatches}.weakDBC.alpha = ...
%             norm(eval(['Dm' num2str(iPatches)]))*polOrder/...
%             BSplinePatches{iPatches}.minElArea*...
%             penaltyPrmScale;
% 
%         % Find along which parametric line the weak Dirichlet 
%         % condition is to be imposed
%         isOnXi = false;
%         if BSplinePatchesLM{iPatches}.weakDBC.xiExtension{1}(1,1) == ...
%              BSplinePatchesLM{iPatches}.weakDBC.xiExtension{1}(1,2)
%             isOnXi = true;
%         end
% 
%         % Make a Lagrange Multipliers discretization
%         clear pLambda XiLambda CPLambda isNURBSLambda; 
% 
%         % Find out up to which polynomial degree the Lagrange
%         % Multipliers discretization needs to be increased
%         if isOnXi
%             polOrderPatch =  BSplinePatchesLM{iPatches}.p;
%         else
%             polOrderPatch =  BSplinePatchesLM{iPatches}.q;
%         end
%         pLM = polOrderPatch - 1;
% 
%         if pLM <= 0
%             pLambda = 0;
%             XiLambda = [0 1];
%             CPLambda = zeros(1,4);
%         else
%             pLambda = 1;
%             XiLambda = [0 0 1 1];
%             CPLambda = zeros(2,4);
% 
%             % Peform a p-refinement
%             tpLambda = polOrderPatch - pLM;
%             [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
%                 (pLambda,XiLambda,CPLambda,tpLambda,'');
%         end
%         isNURBSLambda = 0;
%         nxiLambda = length(CPLambda(:,1));
%         for i = 1:nxiLambda
%             if CPLambda(i,4) ~= 1
%                 isNURBSLambda = 1;
%                 break;
%             end
%         end
% 
%         % Perform an h-refinement
%         percentage = .6;
%         if isOnXi
%             Rxi = unique(BSplinePatchesLM{iPatches}.Xi);
%         else
%             Rxi = unique(BSplinePatchesLM{iPatches}.Eta);
%         end
%         noXi = ceil(percentage*(length(Rxi) - 2));
%         [XiLambda,CPLambda] = knotRefineUniformlyBSplineCurve...
%             (noXi,pLambda,XiLambda,CPLambda,'');
% 
%         % Create the corresponding Lagrange Multipliers structure
%         BSplinePatchesLM{iPatches}.weakDBC.lambda{1} = fillUpLagrangeMultipliers...
%          (pLambda,XiLambda,CPLambda,isNURBSLambda);
%     end
% end
% 
% % Assign the parameters for multipatch coupling
% % ---------------------------------------------
% 
% propCouplingLM.method = 'lagrangeMultipliers';
% propCouplingLM.alphaD = zeros(connections.No,1);
% propCouplingLM.alphaR = zeros(connections.No,1);
% propCouplingLM.intC = intC;
% for iConnections = 1:connections.No
%     % Get the id's of the patches
%     IDPatchI = connectionsLM.xiEtaCoup(iConnections,1);
%     IDPatchJ = connectionsLM.xiEtaCoup(iConnections,2);
% 
%     % Get the mean polynomial order between the patches
%     isOnXiI = false;
%     if connectionsLM.xiEtaCoup(iConnections,5) == connectionsLM.xiEtaCoup(iConnections,6)
%         isOnXiI = true;
%     end
%     if isOnXiI
%         polOrderI = BSplinePatchesLM{IDPatchI}.p;
%     else
%         polOrderI = BSplinePatchesLM{IDPatchI}.q;
%     end
%     isOnXiJ = false;
%     if connectionsLM.xiEtaCoup(iConnections,9) == connectionsLM.xiEtaCoup(iConnections,10)
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
% return;
% save data_FoFiDDMFourPointSailCoarse;

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

%% Compute the relative error from the reference results
% referenceData = importdata('../../../preComputedData/isogeometricMembraneAnalysis/referenceSolution_fourPointSailp3q3Xi50Eta50.mat');
% propReferenceSolution.referenceBSplinePatch = referenceData.BSplinePatch;
% propNewtonRapshon.eps = 1e-9;
% propNewtonRapshon.maxIt = 10;
% propError.noSamplingPoints = 10;
% propError.tolClose = 1e-4;
% propInt.type = 'default';
% [relGeoDomainErrL2,relGeoInterfaceErrL2] = ....
%     computeDomainAndInterfaceErrorInL2NormMembraneFormFiding...
%     (BSplinePatches,connections,propReferenceSolution,propNewtonRapshon,...
%     propError,propInt,'outputEnabled');

%% Plot the reference configuration for the multipatch geometry after the form-finding analysis
color = [.85098 .8549 .85882];
% color = 'none';
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = 260;
el = 30;
view(az,el);
camlight(30,70);
lighting phong;
% axis off;
title('');
% limits = [-7.5 2.5 -5 2.5 0 10]
% axis(limits);
axis off;

%% END