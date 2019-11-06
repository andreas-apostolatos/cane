%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script documentation
% 
% Task : Write a Carat++ case for a 4 patch semisphere with singularities
%        on the equator
%
% Date : 08.12.2017
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

% Define the patches by polynomial orders, knot vectors and control point
% coordinates/weights
for iPatches = 1

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

    % z-coordinates for an exact semisphere
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

    % z-coordinates for an exact semisphere
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

    % z-coordinates for an exact semisphere
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

    % z-coordinates for an exact semisphere
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
end

% Number of patches
noPatches = 4;

% Rotate each patch by 90 degrees with respect to x-y plane
theta = -pi/2;
transformationMatrix = [cos(theta) -sin(theta)
                        sin(theta) cos(theta)];
for iPatches = 1:noPatches
    CP = eval(['CP' num2str(iPatches)]);
    nxi = length(CP(:,1,1));
    neta = length(CP(:,1,1));
    CPRot = zeros(nxi,neta,4);
    CPRot(:,:,[3 4]) = CP(:,:,[3 4]);
    for iEta = 1:neta
        for iXi = 1:nxi
            CPRot(iXi,iEta,1:2) = transformationMatrix*squeeze(CP(iXi,iEta,1:2));
        end
    end
    assignin('base',['CP' num2str(iPatches)],CPRot);
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
                     
% Young's modulus
parameters.E = EYoung;

% Poisson ratio
parameters.nue = poissonRatio;

% Thickness of the plate
parameters.t = thickness;

% Density of the membrane (used only for dynamics)
parameters.rho = density;

% Prestress for the membrane
parameters.prestress = prestress;
                     
% Loop over all patches
for iPatches = 1:noPatches
    assignin('base',['parameters' num2str(iPatches)],parameters);
end

%% GUI

% Case name
caseName = 'DDM4PatchesSemisphere';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points

% Integration rule for the patches
int.type = 'default';
if strcmp(int.type,'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.etaNGPForLoad = 6;
    int.nGPForLoad = 6;
    int.nGPError = 12;
end

% Loop over all patches
for iPatches = 1:noPatches
    assignin('base',['int' num2str(iPatches)],int);
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

% Material matrices
Dm = parameters.E*parameters.t/(1 - parameters.nue^2)*...
      [1              parameters.nue 0
       parameters.nue 1              0
       0              0              (1 - parameters.nue)/2];
Db = parameters.t^2/12*Dm;

% Loop over all patches
for iPatches = 1:noPatches
    assignin('base',['Dm' num2str(iPatches)],Dm);
    assignin('base',['Db' num2str(iPatches)],Db);
end

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

% Refinement parameter
iPRef = 0; % 1

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = iPRef + 1;
tp1 = a;
tq1 = a;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% Patch 2 :
% _________

b = iPRef;
tp2 = b;
tq2 = b;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');

% Patch 3 :
% _________

c = iPRef + 1;
tp3 = c;
tq3 = c;
[Xi3,Eta3,CP3,p3,q3] = degreeElevateBSplineSurface(p3,q3,Xi3,Eta3,CP3,tp3,tq3,'');

% Patch 4 :
% _________

d = iPRef;
tp4 = d;
tq4 = d;
[Xi4,Eta4,CP4,p4,q4] = degreeElevateBSplineSurface(p4,q4,Xi4,Eta4,CP4,tp4,tq4,'');

%%%%%%%%%%%%%%%%%%%
% Knot insertion  %
%%%%%%%%%%%%%%%%%%%

% Refinement parameters
iHRef = 14; % ["4" "10" 14]
noHRefStart = 8; % 8

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

% NBC2.noCnd = 1;
% xib2 = 0;   etab2 = .5;
% NBC2.xiLoadExtension = {xib2};
% NBC2.etaLoadExtension = {etab2};
% NBC2.loadAmplitude = {loadAmplitude2};
% % NBC2.loadAmplitude = {@(x,y,z,t) - (t<=5)*FAmp - (t>5)*FAmp*9e-1};
% NBC2.loadDirection = {dirForce};
% NBC2.computeLoadVct{1} = 'comp1uteLoadVctPointIGAThinStructure';
% NBC2.isFollower(1,1) = isFollower;
% NBC2.isTimeDependent(1,1) = false;

% Patch 3 :
% _________

% NBC3.noCnd = 1;
% xib3 = [0 1];   etab3 = [0 1];
% NBC3.xiLoadExtension = {xib3};
% NBC3.etaLoadExtension = {etab3};
% NBC3.loadAmplitude = {loadAmplitude1};
% % NBC3.loadAmplitude = {@(x,y,z,t) (t<=5)*FAmp + (t>5)*FAmp*9e-1};
% NBC3.loadDirection = {dirForce};
% NBC3.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
% NBC3.isFollower(1,1) = isFollower;
% NBC3.isTimeDependent(1,1) = false;

% NBC3.noCnd = 1;
% xib3 = 0;   etab3 = .5;
% NBC3.xiLoadExtension = {xib3};
% NBC3.etaLoadExtension = {etab3};
% NBC3.loadAmplitude = {loadAmplitude1};
% % NBC3.loadAmplitude = {@(x,y,z,t) (t<=5)*FAmp + (t>5)*FAmp*9e-1};
% NBC3.loadDirection = {dirForce};
% NBC3.computeLoadVct{1} = 'computeLoadVctPointIGAThinStructure';
% NBC3.isFollower(1,1) = isFollower;
% NBC3.isTimeDependent(1,1) = false;

NBC3.noCnd = 1;
xib3 = 0.341081377402109;   etab3 = 0.602802323379784;
NBC3.xiLoadExtension = {xib3};
NBC3.etaLoadExtension = {etab3};
NBC3.loadAmplitude = {loadAmplitude1};
% NBC3.loadAmplitude = {@(x,y,z,t) (t<=5)*FAmp + (t>5)*FAmp*9e-1};
NBC3.loadDirection = {dirForce};
NBC3.computeLoadVct{1} = 'computeLoadVctPointIGAThinStructure';
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

%% Create the patches
BSplinePatches = struct([]);
for iPatches = 1:noPatches
    p = eval(['p' num2str(iPatches)]);
    Xi = eval(['Xi' num2str(iPatches)]);
    q = eval(['q' num2str(iPatches)]);
    Eta = eval(['Eta' num2str(iPatches)]);
    CP = eval(['CP' num2str(iPatches)]);
    isNURBS = eval(['isNURBS' num2str(iPatches)]);
    homDOFs = eval(['homDOFs' num2str(iPatches)]);
    inhomDOFs = eval(['inhomDOFs' num2str(iPatches)]);
    valuesInhomDOFs = eval(['valuesInhomDOFs' num2str(iPatches)]);
    weakDBC = eval(['weakDBC' num2str(iPatches)]);
    cables = eval(['cables' num2str(iPatches)]);
    NBCPatch = eval(['NBC' num2str(iPatches)]);
    xicoup = [];
    etacoup = [];
    int = eval(['int' num2str(iPatches)]);
    patch = fillUpPatch...
        (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,...
        inhomDOFs,valuesInhomDOFs,weakDBC,cables,NBCPatch,[],[],[],...
        xicoup,etacoup,int);
    assignin('base',['patch' num2str(iPatches)],patch);
    BSplinePatches{iPatches} = patch;
end

%% Project the node on the multipatch geometry
% node170 = sqrt(2)*Radius/2*[1.0; 0; 1.0];
theta = pi/4;
phi = pi/4;
x = Radius*cos(theta)*cos(phi);
y = Radius*cos(theta)*sin(phi);
z = Radius*sin(theta);
node170 = [x; y; z];
parametricCoordinatesOnPatch = zeros(2,noPatches);
for iPatches = 1:noPatches
    xi0 = (BSplinePatches{iPatches}.Xi(1) + BSplinePatches{iPatches}.Xi(end))/2;
    eta0 = (BSplinePatches{iPatches}.Eta(1) + BSplinePatches{iPatches}.Eta(end))/2;
    propNewtonRaphson.eps = 1e-15;
    propNewtonRaphson.maxIt = 1e2;
    [xi,eta,Projected,isProjected,noIter] = computeNearestPointProjectionOnBSplineSurface...
        (node170,BSplinePatches{iPatches}.p,BSplinePatches{iPatches}.Xi,BSplinePatches{iPatches}.q,...
        BSplinePatches{iPatches}.Eta,BSplinePatches{iPatches}.CP,BSplinePatches{iPatches}.isNURBS,...
        xi0,eta0,propNewtonRaphson);
    if isProjected && norm(Projected - node170) < 1e-6
        parametricCoordinatesOnPatch(:,iPatches) = [xi; eta];
    end
end

%% Assign the condition enforcement properties using Penalty

% Properties for the weak Dirichlet boundary conditions
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
        BSplinePatches{iPatches}.minElArea;
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
    
    % Get the connection information
    xicoupIJ = connections.xiEtaCoup(iConnections,3:4);
    etacoupIJ = connections.xiEtaCoup(iConnections,5:6);
    xicoupJI = connections.xiEtaCoup(iConnections,7:8);
    etacoupJI = connections.xiEtaCoup(iConnections,9:10);

    % Get the mean polynomial order between the patches
    isOnXiI = false;
    if etacoupIJ(1,1) == etacoupIJ(1,2)
        isOnXiI = true;
    end
    if isOnXiI
        polOrderI = BSplinePatches{IDPatchI}.p;
    else
        polOrderI = BSplinePatches{IDPatchI}.q;
    end
    isOnXiJ = false;
    if etacoupJI(1,1) == etacoupJI(1,2)
        isOnXiJ = true;
    end
    if isOnXiJ
        polOrderJ = BSplinePatches{IDPatchJ}.p;
    else
        polOrderJ = BSplinePatches{IDPatchJ}.q;
    end
    polOrderMean = mean([polOrderI polOrderJ]);  

    % Assign the penalty parameters
    propCouplingPenalty.alphaD(iConnections,1) = ...
        max([norm(eval(['Dm' num2str(IDPatchI)])) ...
        norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
        min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea]);
    propCouplingPenalty.alphaR(iConnections,1) = 0;
end

%% Plot the multipatch geometry with the patch numbering
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% color = [.85098 .8549 .85882];
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,color,graph);

%% Compute the load vector for the visualization of the reference configuration
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.FGamma = zeros(3*BSplinePatches{iPatches}.noCPs,1);
%     for iNBC = 1:BSplinePatches{iPatches}.NBC.noCnd
%         BSplinePatches{iPatches}.FGamma = zeros(3*BSplinePatches{iPatches}.noCPs,1);
%         funcHandle = str2func(NBC{iPatches}.computeLoadVct{iNBC});
%         BSplinePatches{iPatches}.FGamma = ...
%             funcHandle(BSplinePatches{iPatches}.FGamma,BSplinePatches{iPatches},...
%             NBC{iPatches}.xiLoadExtension{iNBC},...
%             NBC{iPatches}.etaLoadExtension{iNBC},...
%             NBC{iPatches}.loadAmplitude{iNBC},...
%             NBC{iPatches}.loadDirection{iNBC},...
%             NBC{iPatches}.isFollower(iNBC,1),...
%             0,int,'outputEnabled');
%     end
end

%% Plot the reference configuration
color = [217 218 219]/255;
% color = 'none';
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = -40; % -310, -40
el = 40;
axis off;
view(az,el);
camlight(0,0);
lighting phong;
title('');
% print('semisphereSetting','-dpng','-r300');

%% Plot the deviation from an exact semisphere
propSphere.radius = .075;
propSphere.center = [0; 0; 0];
[relErrorInL2Domain,relErrorInL2WeakDBC,errorInL2Interface,graph.index] = ...
    plot_deviationFromSphereBSPlineSurface...
    (BSplinePatches,'undefined',connections,propSphere,int,graph,'outputEnabled');

%% Write the case for Carat++
for iPatches = 1:noPatches
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
    weakDBC.alpha = norm(eval('Dm'))*polOrder/...
        BSplinePatches{iPatches}.minElArea;
end
strongDBC = struct([]);
for iPatches = 1:noPatches
    strongDBCPatch.xiExtensionHom = eval(['xiSup' num2str(iPatches)]);
    strongDBCPatch.etaExtensionHom = eval(['etaSup' num2str(iPatches)]);
    strongDBC{iPatches} = strongDBCPatch;
end
pathToOutput = '../../../inputCarat/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4Carat(BSplinePatches,strongDBC,connections,pathToOutput,caseName);
noElements = 0;
for iPatches = 1:noPatches
    noElements = noElements + BSplinePatches{iPatches}.noElmnts;
end
fprintf('Number of elements on the multipach geometry equals %d\n',noElements);
noCPs = 0;
for iPatches = 1:noPatches
    noCPs = noCPs + BSplinePatches{iPatches}.noCPs;
end
fprintf('Number of Control Points on the multipach geometry equals %d\n',noCPs);
return;

%% Compute the displacement field at the middle of the patch using Carat

% Read the displacement field from Carat
% displacementCarat = importdata('../../../outputIBRACarat/semisphereStatic_IGA4Patches__noElmnts219LowPolynomialOrder/displacements');
% displacementCarat = importdata('../../../outputIBRACarat/semisphereStatic_IGA4Patches_noElmnts504LowPolynomialOrder/displacements');
displacementCarat = importdata('../../../outputIBRACarat/semisphereStatic_IGA4Patches_noElmnts729LowPolynomialOrder/displacements');
% displacementCarat = importdata('../../../outputIBRACarat/semisphereStatic_IGA4Patches_noElmnts504HighPolynomialOrder/displacements');

% Initialize reconstructed displacement field
displacementCaratReconstructed = zeros(noCPs,length(displacementCarat(1,:)));

% Reconstruct the displacement field with the strong boundary conditions
ids = displacementCarat(:,1);
displacementCaratReconstructed(ids,:) = displacementCarat;
displacementCaratReconstructed = displacementCaratReconstructed(:,2:4);

% Distribute the displacement field into a contiguous array
noCPs = length(displacementCaratReconstructed(:,1));
noDOFs = 3*noCPs;
dHatCarat = zeros(noDOFs,1);
for iCPs  = 1:noCPs
    dHatCarat(3*iCPs - 2,1) = displacementCaratReconstructed(iCPs,1);
    dHatCarat(3*iCPs - 1,1) = displacementCaratReconstructed(iCPs,2);
    dHatCarat(3*iCPs,1) = displacementCaratReconstructed(iCPs,3);
end

% Plot the current configuration and the resultants
scaling = 1e0;
graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
    (BSplinePatches,scaling*dHatCarat,graph,'outputEnabled');

% Initialize array containing the displacement field of each patch at the
% top
dMiddle = zeros(3,noPatches);
xiPatches = [BSplinePatches{1}.Xi(1)
             BSplinePatches{2}.Xi(1) 
             BSplinePatches{3}.Xi(1) 
             BSplinePatches{4}.Xi(1)];
etaPatches = [BSplinePatches{1}.Eta(end) 
              BSplinePatches{2}.Eta(end) 
              BSplinePatches{3}.Eta(end) 
              BSplinePatches{4}.Eta(end)];
          
% Create a freedom table for each patch
noDOFs = 0;
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.EFTPatches = (1 + noDOFs):(3*BSplinePatches{iPatches}.noCPs + noDOFs);
    noDOFs = noDOFs + 3*BSplinePatches{iPatches}.noCPs;
end    

% Loop over all patches
for iPatches = 1:noPatches
    % Recover patch data
    p = BSplinePatches{iPatches}.p;
    q = BSplinePatches{iPatches}.q;
    Xi = BSplinePatches{iPatches}.Xi;
    Eta = BSplinePatches{iPatches}.Eta;
    CP = BSplinePatches{iPatches}.CP;
    isNURBS = BSplinePatches{iPatches}.isNURBS;
    
    % Recover the displacement field of the patch
    dHatCaratPatch = dHatCarat(BSplinePatches{iPatches}.EFTPatches);
    
    % Make a DOF numbering for the given patch
    mxi = length(Xi);
    meta = length(Eta);
    nxi = length(CP(:,1,1));
    neta = length(CP(1,:,1));
    dHatElem = zeros(mxi - p - 1,meta - q - 1,3*(p + 1)*(q + 1));
    for etaSpan = (q + 1):(meta - q - 1)
        for xiSpan = (p + 1):(mxi - p - 1)
            xiCounter = 1; 
            for c = etaSpan-q-1:etaSpan-1 
                for b = xiSpan-p:xiSpan
                    dHatElem(xiSpan,etaSpan,xiCounter) = dHatCaratPatch(3*(c*nxi + b)-2);
                    dHatElem(xiSpan,etaSpan,xiCounter + 1) = dHatCaratPatch(3*(c*nxi + b)-1);
                    dHatElem(xiSpan,etaSpan,xiCounter + 2) = dHatCaratPatch(3*(c*nxi + b));

                    % Update counter
                    xiCounter = xiCounter + 3;
                end
            end
        end
    end
    
    % Get the displacement vector at the top of the semisphere
    xi = xiPatches(iPatches,1);
    eta = etaPatches(iPatches,1);
    xiSpan = findKnotSpan(xi,Xi,nxi);
    etaSpan = findKnotSpan(eta,Eta,neta);
    dHatActual = squeeze(dHatElem(xiSpan,etaSpan,:));
    R = computeIGABasisFunctionsAndDerivativesForSurface...
        (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
    X = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R);
    dMiddle(:,iPatches) = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,R,dHatActual);
end

% Compute the mean of the vertical displacement at the top of the 
% semisphere
dMiddleMean = mean(dMiddle(3,:))

%% END