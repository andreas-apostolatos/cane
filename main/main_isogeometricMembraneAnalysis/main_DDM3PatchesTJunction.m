%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : A T-junction is modelled is modelled with 3 membrane patches
%        while the one end is subject to boundary root point excitation
%
% Date : 21.05.2017
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

% Add functions related to the efficient computation
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
        '../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../isogeometricThinStructureAnalysis/output/',...
        '../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../isogeometricThinStructureAnalysis/weakDBCMembrane/');

% Add functions related to the visualization of the multipatch geometry
% related to the isogeometric mortar-based mapping
addpath('../../isogeometricMortarBasedMappingAnalysis/graphics/');

%% CAD modelling via NURBS

% Global variables:
Length = 5;
Width = 2;
Height = 2;

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
CP1(:,:,1) = [-Length/2 -Length/2
              0         0];

% y-coordinates
CP1(:,:,2) = [-Width/2 Width/2
              -Width/2 Width/2];

% z-coordinates
CP1(:,:,3) = [Height Height
              0      0];
          
% Weights
CP1(:,:,4) = [1 1
              1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = 0;
nxi1 = length(CP1(:,1,1));
neta1 = length(CP1(1,:,1));
for i= 1:nxi1
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
p2 = 1;
q2 = 1;

% Knot vectors
Xi2 = [0 0 1 1];
Eta2 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP2(:,:,1) = [-Length/2 -Length/2
              0         0];

% y-coordinates
CP2(:,:,2) = [-Width/2 Width/2
              -Width/2 Width/2];

% z-coordinates
CP2(:,:,3) = [-Height -Height
              0       0];
          
% Weights
CP2(:,:,4) = [1 1
              1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS2 = 0;
nxi2 = length(CP2(:,1,1));
neta2 = length(CP2(1,:,1));
for i= 1:nxi2
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
p3 = 1;
q3 = 1;

% Knot vectors
Xi3 = [0 0 1 1];
Eta3 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP3(:,:,1) = [0        0
              Length/2 Length/2];

% y-coordinates
CP3(:,:,2) = [-Width/2 Width/2
              -Width/2 Width/2];

% z-coordinates
CP3(:,:,3) = [0      0
              Height Height];
          
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

% Patch 4 :
% _________

% Polynomial degrees
p4 = 1;
q4 = 1;

% Knot vectors
Xi4 = [0 0 1 1];
Eta4 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP4(:,:,1) = [0        0
              Length/2 Length/2];

% y-coordinates
CP4(:,:,2) = [-Width/2 Width/2
              -Width/2 Width/2];

% z-coordinates
CP4(:,:,3) = [0       0
              -Height -Height];
          
% Weights
CP4(:,:,4) = [1 1
              1 1];

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

%% Material constants

% general parameters
EYoung = 8e+8;
nue = 0.0;
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
int.type = 'default';
if strcmp(int.type,'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.etaNGPForLoad = 6;
    int.nGPForLoad = 6;
    int.nGPError = 12;
end
int1 = int;
int2 = int;
int3 = int;
int4 = int;

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

tp1 = 1;
tq1 = 1;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% Patch 2 :
% _________

tp2 = 1;
tq2 = 1;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');

% Patch 3 :
% _________

tp3 = 1;
tq3 = 1;
[Xi3,Eta3,CP3,p3,q3] = degreeElevateBSplineSurface(p3,q3,Xi3,Eta3,CP3,tp3,tq3,'');

% Patch 4 :
% _________

tp4 = 1;
tq4 = 1;
[Xi4,Eta4,CP4,p4,q4] = degreeElevateBSplineSurface(p4,q4,Xi4,Eta4,CP4,tp4,tq4,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = 8; % 12
nKnotsXi1 = a;
nKnotsEta1 = ceil(a*Width/Length);
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface(p1,Xi1,q1,Eta1,CP1,nKnotsXi1,nKnotsEta1,'');

% Patch 2 :
% _________

b = 8; % 14
nKnotsXi2 = b;
nKnotsEta2 = ceil(b*Width/Length);
[Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface(p2,Xi2,q2,Eta2,CP2,nKnotsXi2,nKnotsEta2,'');

% Patch 3 :
% _________

b = 8; % 16
nKnotsXi3 = b;
nKnotsEta3 = ceil(b*Width/Length);
[Xi3,Eta3,CP3] = knotRefineUniformlyBSplineSurface(p3,Xi3,q3,Eta3,CP3,nKnotsXi3,nKnotsEta3,'');

% Patch 4 :
% _________

b = 8; % 14
nKnotsXi4 = b;
nKnotsEta4 = ceil(b*Width/Length);
[Xi4,Eta4,CP4] = knotRefineUniformlyBSplineSurface(p4,Xi4,q4,Eta4,CP4,nKnotsXi4,nKnotsEta4,'');

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous DOFs
homDOFs1 = [];
xisup1 = [0 0];
etasup1 = [0 1];
for dirSupp = 1:3
    homDOFs1 = findDofs3D(homDOFs1,xisup1,etasup1,dirSupp,CP1);
end

% Inhomogeneous DOFs
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak Dirichlet boundary conditions
weakDBC1 = [];
weakDBC1.noCnd = 0;

% Cables
cables1.No = 0;

% Patch 2 :
% _________

% Homogeneous DOFs
homDOFs2 = [];
xisup2 = [0 0];   
etasup2 = [0 1];
for dirSupp = 1:3
    homDOFs2 = findDofs3D(homDOFs2,xisup2,etasup2,dirSupp,CP2);
end

% Inhomogeneous DOFs
inhomDOFs2 = [];
valuesInhomDOFs2 = [];

% Weak Dirichlet boundary conditions
weakDBC2 = [];
weakDBC2.noCnd = 0;

% Cables
cables2.No = 0;

% Patch 3 :
% _________

% Homogeneous DOFs
homDOFs3 = [];

% Inhomogeneous DOFs
inhomDOFs3 = [];
xiSup3 = [1 1];   etaSup3 = [0 1];
for dirSupp = 1:3
    inhomDOFs3 = findDofs3D(inhomDOFs3,xiSup3,etaSup3,dirSupp,CP3);
end
valuesInhomDOFs3 = cell({});
for iInhomDBC = 1:length(inhomDOFs3)/3
    valuesInhomDOFs3{3*iInhomDBC - 2} = @(t)4;
    valuesInhomDOFs3{3*iInhomDBC - 1} = @(t)2;
    valuesInhomDOFs3{3*iInhomDBC} = @(t)-2;
end

% Weak Dirichlet boundary conditions
weakDBC3 = [];
weakDBC3.noCnd = 0;

% Cables
cables3.No = 0;

% Patch 4 :
% _________

% Homogeneous DOFs
homDOFs4 = [];

% Inhomogeneous DOFs
inhomDOFs4 = [];
xiSup4 = [1 1];   etaSup4 = [0 1];
for dirSupp = 1:3
    inhomDOFs4 = findDofs3D(inhomDOFs4,xiSup4,etaSup4,dirSupp,CP4);
end
valuesInhomDOFs4 = cell({});
for iInhomDBC = 1:length(inhomDOFs4)/3
    valuesInhomDOFs4{3*iInhomDBC - 2} = @(t)4;
    valuesInhomDOFs4{3*iInhomDBC - 1} = @(t)2;
    valuesInhomDOFs4{3*iInhomDBC} = @(t)-2;
end

% Weak Dirichlet boundary conditions
weakDBC4 = [];
weakDBC4.noCnd = 0;

% Cables
cables4.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load amplitude
FAmp = + 1e1;

% Patch 1 :
% _________

% FAmp1 = FAmp*0;
NBC1.noCnd = 0;

% Patch 2 :
% _________

% FAmp2 = FAmp;
NBC2.noCnd = 0;

% Patch 3 :
% _________

% FAmp3 = FAmp;
NBC3.noCnd = 0;

% Patch 4 :
% _________

% FAmp3 = FAmp;
NBC4.noCnd = 0;

% Collect all the Neumann boundary conditions into an arra<y
NBC = {NBC1 NBC2 NBC3 NBC4};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface parametrization     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Connection with patch 2:
xicoup12 = [1 1];
etacoup12 = [0 1];

% Connection with patch 3:
xicoup13 = [1 1];
etacoup13 = [0 1];

% Connection with patch 4:
xicoup14 = [1 1];
etacoup14 = [0 1];

% Collect all interfaces into arrays:
xicoup1 = [xicoup12 xicoup13 xicoup14];
etacoup1 = [etacoup12 etacoup13 etacoup14];

% Patch 2 :
% _________

% Connection with patch 1:
xicoup21 = [1 1];
etacoup21 = [0 1];

% Connection with patch 3:
xicoup23 = [1 1];
etacoup23 = [0 1];

% Connection with patch 4:
xicoup24 = [1 1];
etacoup24 = [0 1];

% Collect all interfaces into arrays:
xicoup2 = [xicoup21 xicoup23 xicoup24];
etacoup2 = [etacoup21 etacoup23 etacoup24];

% Patch 3 :
% _________

% Connection with patch 1:
xicoup31 = [0 0];
etacoup31 = [0 1];

% Connection with patch 3:
xicoup32 = [0 0];
etacoup32 = [0 1];

% Connection with patch 4:
xicoup34 = [0 0];
etacoup34 = [0 1];

% Collect all interfaces into arrays:
xicoup3 = [xicoup31 xicoup32 xicoup34];
etacoup3 = [etacoup31 etacoup32 etacoup34];

% Patch 4 :
% _________

% Connection with patch 1:
xicoup41 = [0 0];
etacoup41 = [0 1];

% Connection with patch 2:
xicoup42 = [0 0];
etacoup42 = [0 1];

% Connection with patch 3:
xicoup43 = [0 0];
etacoup43 = [0 1];

% Collect all interfaces into arrays:
xicoup4 = [xicoup41 xicoup42 xicoup43];
etacoup4 = [etacoup41 etacoup42 etacoup43];

% Define connections by patch numbers
connections.xiEtaCoup(:,:) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21
                              1 3 xicoup13 etacoup13 xicoup31 etacoup31
                              1 4 xicoup14 etacoup14 xicoup41 etacoup41
                              2 3 xicoup23 etacoup23 xicoup32 etacoup32
                              2 4 xicoup24 etacoup24 xicoup42 etacoup42
                              3 4 xicoup34 etacoup34 xicoup43 etacoup43];
connections.masterSlave = [true
                           true
                           true
                           false
                           false
                           false];
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

%% Plot reference configuration
color = [.85098 .8549 .85882];
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');

%% Plot the multipatch geometry together with the parametric coordinates
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,graph);

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
    propCoupling.isSlaveSideCoarser = false;
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