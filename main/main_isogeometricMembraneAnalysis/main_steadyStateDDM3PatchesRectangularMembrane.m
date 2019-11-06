%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Steady-state nonlinear analysis over a 3-patch four point sail
%        subject to snow load.
%
% Date : 10.01.2017
%
%% Preamble
clear;
clc;

%% Includes 

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Add system solvers
addpath('../../equationSystemSolvers/');

% Add efficient computation functions
addpath('../../efficientComputation/');

% Add transient analysis solvers
addpath('../../transientAnalysis/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the isogeometric Kirchhoff-Love shell formulation
addpath('../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../isogeometricThinStructureAnalysis/graphicsMultipatches/',...
        '../../isogeometricThinStructureAnalysis/loads/',...
        '../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../isogeometricThinStructureAnalysis/solvers/',...
        '../../isogeometricThinStructureAnalysis/metrics/',...
        '../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../isogeometricThinStructureAnalysis/penaltyDecompositionKLShell/',...
        '../../isogeometricThinStructureAnalysis/penaltyDecompositionMembrane/',...
        '../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionKLShell/',...
        '../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionMembrane/',...
        '../../isogeometricThinStructureAnalysis/nitscheDecompositionKLShell/',...
        '../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../isogeometricThinStructureAnalysis/precomputedData/',...
        '../../isogeometricThinStructureAnalysis/output/',...
        '../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../isogeometricThinStructureAnalysis/formFindingAnalysis/');
    
% Add functions related to the visualization of the multipatch geometry
% related to the isogeometric mortar-based mapping
addpath('../../isogeometricMortarBasedMappingAnalysis/graphics/');

% Check if the operating system is Windows and if yes add some functions
% from a temporary repository
if ispc
    addpath('../../functions4Windows/');
end

%% CAD modelling of the geometry via NURBS

% Global variables:
Length = 4;
Width = 4;

% Patch 1 :
% _________

% Polynomial degrees
p1 = 1;
q1 = 1;

% Knot vectors
Xi1 = [0 0 1 1];
Eta1 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP1(:,:,1) = [0   0
              1/3 1/3]*Length;

% y-coordinates
CP1(:,:,2) = [-1/2  1/2
              -1/2  1/2]*Width;

% z-coordinates
CP1(:,:,3) = [0 0
              0 0];

% Weights
CP1(:,:,4) = [1 1
              1 1];
          
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = 0;
nxi1 = length(CP1(:,1,1));
neta2 = length(CP1(1,:,1));
for i= 1:nxi1
    for j=1:neta2
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
p2 = 1;
q2 = 1;

% Knot vectors
Xi2 = [0 0 1 1];
Eta2 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP2(:,:,1) = [1/3 1/3
              2/3 2/3]*Length;

% y-coordinates
CP2(:,:,2) = [-1/2 1/2
              -1/2 1/2]*Width;

% z-coordinates
CP2(:,:,3) = [0 0
              0 0];

% Weights
CP2(:,:,4) = [1 1
              1 1];
          
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS2 = 0;
nxi2 = length(CP2(:,1,1));
neta2 = length(CP2(1,:,1));
for i= 1:nxi2
    for j=1:neta2
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
p3 = 1;
q3 = 1;

% Knot vectors
Xi3 = [0 0 1 1];
Eta3 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP3(:,:,1) = [2/3 2/3
              3/3 3/3]*Length;

% y-coordinates
CP3(:,:,2) = [-1/2 1/2
              -1/2 1/2]*Width;

% z-coordinates
CP3(:,:,3) = [0 0
              0 0];

% Weights
CP3(:,:,4) = [1 1
              1 1];
          
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

%% Material constants

% general parameters
EYoung = 8e+4;
nue = .4;
thickness = 1e-3;
sigma0 = 1e3; % 1e3
prestress.voigtVector = [sigma0
                         sigma0
                         0];
density = 800;

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

%% GUI

% Case name
caseName = 'DDM4PatchesCircularPlateDistributedLoad';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Choose a method for the application of weak Dirichlet boundary conditions and the multipatch coupling
method = 'Mortar';
if ~strcmp(method,'Penalty') && ~strcmp(method,'LagrangeMultipliers') && ...
        ~strcmp(method,'Mortar') && ~strcmp(method,'AugmentedLagrangeMultipliers') && ...
        ~strcmp(method,'Nitsche')
    error('%s is not a valid method (Nitsche, Penalty, LagrangeMultipliers, AugmentedLagrangeMultipliers)',method);
end

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
graph.postprocConfig = 'current';

% Plot strain or stress field
% .resultant: 'displacement','strain','curvature','force','moment','shearForce'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x','y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

% Function handle to writing out the results
% writeOutput = @writeResults4Carat;
% writeOutput = @writeResults4GiD;
writeOutput = 'undefined';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

%% Refinement

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = 1;
tp1 = a;
tq1 = a;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% Patch 2 :
% _________

b = 0;
tp2 = b;
tq2 = b;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');

% Patch 3 :
% _________

b = 1;
tp3 = b;
tq3 = b;
[Xi3,Eta3,CP3,p3,q3] = degreeElevateBSplineSurface(p3,q3,Xi3,Eta3,CP3,tp3,tq3,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

noKnotsXi1 = 7;
noKnotsEta1 = noKnotsXi1;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface...
    (p1,Xi1,q1,Eta1,CP1,noKnotsXi1,noKnotsEta1,'');

% Patch 2 :
% _________

noKnotsXi2 = 3;
noKnotsEta2 = noKnotsXi2;
[Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface...
    (p2,Xi2,q2,Eta2,CP2,noKnotsXi2,noKnotsEta2,'');

% Patch 3 :
% _________

noKnotsXi3 = 9;
noKnotsEta3 = noKnotsXi3;
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
xisup1 = [0 1];   etasup1 = [0 0];
for dir = 1:3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end
xisup1 = [0 1];   etasup1 = [1 1];
for dir = 1:3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end

% Fix all z-coordinates
% xisup1 = [0 1];   etasup1 = [0 1];
% for dir = 3
%     homDOFs1 = findDofs3D...
%         (homDOFs1,xisup1,etasup1,dir,CP1);
% end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC1.noCnd = 0;

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
xisup2 = [0 1];   etasup2 = [1 1];
for dir = 1:3
    homDOFs2 = findDofs3D...
        (homDOFs2,xisup2,etasup2,dir,CP2);
end

% Fix all z-coordinates
% xisup2 = [0 1];   etasup2 = [0 1];
% for dir = 3
%     homDOFs2 = findDofs3D...
%         (homDOFs2,xisup2,etasup2,dir,CP2);
% end

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
xisup3 = [0 1];   etasup3 = [0 0];
for dir = 1:3
    homDOFs3 = findDofs3D...
        (homDOFs3,xisup3,etasup3,dir,CP3);
end
xisup3 = [0 1];   etasup3 = [1 1];
for dir = 1:3
    homDOFs3 = findDofs3D...
        (homDOFs3,xisup3,etasup3,dir,CP3);
end
xisup3 = [1 1];   etasup3 = [0 1];
for dir = 1:3
    homDOFs3 = findDofs3D... 
        (homDOFs3,xisup3,etasup3,dir,CP3);
end

% Fix all z-coordinates
% xisup3 = [0 1];   etasup3 = [0 1];
% for dir = 3
%     homDOFs3 = findDofs3D...
%         (homDOFs3,xisup3,etasup3,dir,CP3);
% end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs3 = [];
valuesInhomDOFs3 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC3.noCnd = 0;

% Embedded cables
cables3.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parameter
scaling = 1e-2; % 5e1
FAmp = -1e3*scaling;
traction1 = FAmp;
traction2 = FAmp;
traction3 = FAmp;

% Patch 1 :
% _________

% FAmp1 = loadAmplitude;
NBC1.noCnd = 1;
xib1 = [0 1];   etab1 = [0 1];   dirForce1 = 'z';
NBC1.xiLoadExtension = {xib1};
NBC1.etaLoadExtension = {etab1};
NBC1.loadAmplitude = {traction1};
NBC1.loadDirection = {dirForce1};
NBC1.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC1.isFollower = false;
NBC1.isTimeDependent = false;

% Patch 2 :
% _________

% FAmp2 = - loadAmplitude;
NBC2.noCnd = 1;
xib2 = [0 1];   etab2 = [0 1];   dirForce2 = 'z';
NBC2.xiLoadExtension = {xib2};
NBC2.etaLoadExtension = {etab2};
NBC2.loadAmplitude = {traction2};
NBC2.loadDirection = {dirForce2};
NBC2.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC2.isFollower = false;
NBC2.isTimeDependent = false;

% Patch 3 :
% _________

% FAmp3 = - loadAmplitude;
NBC3.noCnd = 1;
xib3 = [0 1];   etab3 = [0 1];   dirForce3 = 'z';
NBC3.xiLoadExtension = {xib3};
NBC3.etaLoadExtension = {etab3};
NBC3.loadAmplitude = {traction3};
NBC3.loadDirection = {dirForce3};
NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC3.isFollower = false;
NBC3.isTimeDependent = false;

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
connections.masterSlave = [true
                           true];
connectionsLM = connections;
connectionsALM = connections;

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

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch1 patch2 patch3};
noPatches = length(BSplinePatches);

%% Compute the load vector for the visualization of the reference configuration
% for counterPatches = 1:noPatches
%     BSplinePatches{counterPatches}.FGamma = ...
%         zeros(3*BSplinePatches{counterPatches}.noCPs,1);
% %     for counterNBC = 1:NBC{counterPatches}.noCnd
% %         funcHandle = str2func(NBC{counterPatches}.computeLoadVct{counterNBC});
% %         BSplinePatches{counterPatches}.FGamma = funcHandle...
% %             (BSplinePatches{counterPatches}.FGamma,...
% %             BSplinePatches{counterPatches},...
% %             NBC{counterPatches}.xiLoadExtension{counterNBC},...
% %             NBC{counterPatches}.etaLoadExtension{counterNBC},...
% %             NBC{counterPatches}.loadAmplitude{counterNBC},...
% %             NBC{counterPatches}.loadDirection{counterNBC},...
% %             NBC{counterPatches}.isFollower(counterNBC,1),0,...
% %             BSplinePatches{counterPatches}.int,'outputEnabled');
% %     end
% end

%% Plot reference configuration
% color = [.85098 .8549 .85882];
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');

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

%% Set up the parameters and properties for each method
if strcmp(method,'Penalty')
    % General parameters
    penaltyPrmScale = 1e0;

    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    for iPatches = 1:noPatches
        % Check if weak boundary conditions are to be applied
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
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
        if connections.xiEtaCoup(iConnections,5) == ...
                connections.xiEtaCoup(iConnections,6)
            isOnXiI = true;
        end
        if isOnXiI
            polOrderI = BSplinePatches{IDPatchI}.p;
        else
            polOrderI = BSplinePatches{IDPatchI}.q;
        end
        isOnXiJ = false;
        if connections.xiEtaCoup(iConnections,9) == ...
                connections.xiEtaCoup(iConnections,10)
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
        % Check if weak boundary conditions are to be applied
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
        % Assign the method name
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
elseif strcmp(method,'Mortar')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'mortar';
    propCoupling.isSlaveSideCoarser = true;
    propCoupling.computeRearrangedProblemMtrcs = @computeRearrangedProblemMtrcs4MortarIGAMembrane;
    propCoupling.intC = intC;
elseif strcmp(method,'AugmentedLagrangeMultipliers')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------
    
    % Scaling factor for the augmented Lagrange Multipliers method
    scalePenaltyFactorALM = 1e-2;
    
    for iPatches = 1:noPatches
        % Check if weak boundary conditions are to be applied
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
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
        if connections.xiEtaCoup(iConnections,5) == ...
                connections.xiEtaCoup(iConnections,6)
            isOnXiI = true;
        end
        if isOnXiI
            polOrderI = BSplinePatches{IDPatchI}.p;
        else
            polOrderI = BSplinePatches{IDPatchI}.q;
        end
        isOnXiJ = false;
        if connections.xiEtaCoup(iConnections,9) == ...
                connections.xiEtaCoup(iConnections,10)
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

    % Properties for the weak Dirichlet boundary conditions
    for iPatches = 1:noPatches
        % Check if weak boundary conditions are to be applied
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
        % Assign the method
        BSplinePatches{iPatches}.weakDBC.method = 'nitsche';
        
        % Assign the estimation of the stabilization parameter
        BSplinePatches{iPatches}.weakDBC.estimationStabilPrm = true;
    end

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'nitsche';
    propCoupling.estimationStabilPrm = true;
    propCoupling.gammaTilde = .5;
    propCoupling.intC = intC;
end

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-4;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 100;

%% Solve the nonlinear problem
plot_IGANLinear = 'undefined';
[dHat,CPHistory,resHistory,hasConverged,~,~,~,~,...
    BSplinePatches,propCoupling,~] = ...
 solve_DDMIGAMembraneMultipatchesNLinear...
 (BSplinePatches,connections,propCoupling,...
 propNLinearAnalysis,solve_LinearSystem,plot_IGANLinear,...
 graph,'outputEnabled');
graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
    (BSplinePatches,dHat,graph,'outputEnabled');
% save data_steadyStateDDMLagrangeMultipliersFourPointSail;
return;

%% END OF SCRIPT