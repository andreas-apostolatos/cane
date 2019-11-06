%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Form-finding analysis over the middle tent of the Olympiadach
%        modelled with 5 patches
%
% Date : 10.01.2017
%
%% Preamble
clear;
clc;

%% Includes 

% Add low order basis functions
addpath('../../../basisFunctions/');

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

% Global variables:
Radius = 1;

% Patch 1 :
% _________

% Polynomial degrees
p1 = 4;
q1 = 4;

% Knot vectors
Xi1 = [0 0 0 0 0 1 1 1 1 1];
Eta1 = [0 0 0 0 0 1 1 1 1 1];

% Control Point coordinates
noCPXi1 = length(Xi1) - p1 - 1;
noCPEta1 = length(Eta1) - q1 - 1;
CP1 = zeros(noCPXi1,noCPEta1,4);
noCoord = 3;
CP1(:,:,1) = [4*(1 - sqrt(3))        , - sqrt(2)                , 0                          , sqrt(2)                     , 4*(sqrt(3) - 1)
              sqrt(2)*(sqrt(3) - 4)  , (2 - 3*sqrt(3))/2        , 0                          , (3*sqrt(3) - 2)/2           , sqrt(2)*(4 - sqrt(3))
              4*(1 - 2*sqrt(3))/3    , sqrt(2)*(2*sqrt(3) - 7)/3, 0                          , sqrt(2)*(7 - 2*sqrt(3))/3   , 4*(2*sqrt(3) - 1)/3
              sqrt(2)*(sqrt(3) - 4)  , (2 - 3*sqrt(3))/ 2       , 0                          , (3*sqrt(3) - 2)/ 2          , sqrt(2)*(4 - sqrt(3))
              4*(1 - sqrt(3))        , -sqrt(2)                 , 0                          , sqrt(2)                     , 4*(sqrt(3) - 1)];
CP1(:,:,2) = [4*(1 - sqrt(3))        , sqrt(2)*(sqrt(3) - 4)    , 4*(1 - 2*sqrt(3))/3        , sqrt(2)*(sqrt(3) - 4)       , 4*(1 - sqrt(3))
              -sqrt(2)               , (2 - 3*sqrt(3))/2        , sqrt(2)*(2*sqrt(3) - 7)/3  , (2 - 3*sqrt(3))/2           , -sqrt(2)
              0                      , 0                        , 0                          , 0                           , 0
              sqrt(2)                , (3*sqrt(3) - 2)/2        , sqrt(2)*(7 - 2*sqrt(3))/3  , (3*sqrt(3) - 2)/2           , sqrt(2) 
              4*(sqrt(3) - 1)        , sqrt(2)*(4 - sqrt(3))    , 4*(2*sqrt(3) - 1)/3        , sqrt(2)*(4 - sqrt(3))       , 4*(sqrt(3) - 1)];
CP1(:,:,3) = -[4*(1 - sqrt(3))        , sqrt(2)*(sqrt(3) - 4)    , 4*(1 - 2*sqrt(3))/3        , sqrt(2)*(sqrt(3) - 4)       , 4*(1 - sqrt(3))
               sqrt(2)*(sqrt(3) - 4)  , -(sqrt(3) + 6)/2         , -5*sqrt(6)/3               , -(sqrt(3) + 6)/2            , sqrt(2)*(sqrt(3) - 4)
               4*(1 - 2*sqrt(3))/3    , -5*sqrt(6)/3             , 4*(sqrt(3) - 5)/3          , -5*sqrt(6)/3                , 4*(1 - 2*sqrt(3))/3
               sqrt(2)*(sqrt(3) - 4)  , -(sqrt(3) + 6)/2         , -5*sqrt(6)/3               , -(sqrt(3) + 6)/2            , sqrt(2)*(sqrt(3) - 4)
               4*(1 - sqrt(3))        , sqrt(2)*(sqrt(3) - 4)    , 4*(1 - 2*sqrt(3))/3        , sqrt(2)*(sqrt(3) - 4)       , 4*(1 - sqrt(3))];
CP1(:,:,4) = [4*(3 - sqrt(3))        , sqrt(2)*(3*sqrt(3) - 2)  , 4*(5 - sqrt(3))/3          , sqrt(2)*(3*sqrt(3) - 2)     , 4*(3 - sqrt(3))
              sqrt(2)*(3*sqrt(3) - 2), (sqrt(3) + 6)/2          , sqrt(2)*(sqrt(3) + 6)/3    , (sqrt(3) + 6)/2             , sqrt(2)*(3*sqrt(3) - 2)
              4*(5 - sqrt(3))/3      , sqrt(2)*(sqrt(3) + 6)/3  , 4*(5*sqrt(3) - 1)/9        , sqrt(2)*(sqrt(3) + 6)/3     , 4*(5 - sqrt(3))/3
              sqrt(2)*(3*sqrt(3) - 2), (sqrt(3) + 6)/2          , sqrt(2)*(sqrt(3) + 6)/3     , (sqrt(3) + 6)/2             , sqrt(2)*(3*sqrt(3) - 2)
              4*(3 - sqrt(3))        , sqrt(2)*(3*sqrt(3) - 2)  , 4*(5 - sqrt(3))/3          , sqrt(2)*(3*sqrt(3) - 2)     , 4*(3 - sqrt(3))];
          
          
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = 0;
nxi1 = length(CP1(:,1,1));
neta1 = length(CP1(1,:,1));
for i = 1:nxi1
    for j = 1:neta1
        if CP1(i,j,4) ~= 1
            isNURBS1 = 1;
            break;
        end
    end
    if isNURBS1
        break;
    end
end

% Number of patches
noPatches = 1;

%% Material constants

% General parameters
EYoung = 2.1e11;
thickness = 5e-3;
diameterCables = 14e-3;
radiusCables = diameterCables/2;
areaCables = pi*radiusCables^2;
areaCablesBoundary = 5e-2;
nue = 4e-1;
n0 = 400000;
nTilde0 = 40000;
scalingYoungModulusCables = 1e2;
rhoPlexiglas = 5e2;
rhoSteel = 7.85e3;

% Homogenisation of the membrane density
a75 = 75e-2;
noCablesPerSquare75 = 4;
massPerSquare75 = rhoPlexiglas*thickness*a75^2 + noCablesPerSquare75*rhoSteel*areaCables*a75;
massPerSquare75 = massPerSquare75 + 0.15*massPerSquare75;
rhoMembrane = massPerSquare75/a75^2/thickness;

% Membrane :
% ----------

parameters.E = EYoung;
parameters.nue = nue;
parameters.t = thickness;
parameters.rho = rhoMembrane;
parameters.prestress.voigtVector = [50*n0
                                    50*n0
                                    0];

% Cables :
% --------
                                
% n19021 = nTilde0*1680;
% n19022 = nTilde0*1730;
% n19023 = nTilde0*1595;
% n19024 = nTilde0*1690;
% n19025 = nTilde0*1350;
% n19026 = nTilde0*1250;
% n19027 = nTilde0*1250;
% n19028 = nTilde0*1575;
% n19029 = nTilde0*1550;
% n19030 = nTilde0*1845;
% constNum = 19020;
% noCables = 10;
% for iCables = 1:noCables
%     evalin('base',['parametersCable' num2str(constNum + iCables) '.E = EYoung;']);
%     evalin('base',['parametersCable' num2str(constNum + iCables) '.areaCS = areaCablesBoundary;']);
%     evalin('base',['parametersCable' num2str(constNum + iCables) '.rho = rhoSteel;']);
%     evalin('base',['parametersCable' num2str(constNum + iCables) '.prestress = n' num2str(constNum + iCables) ';']);
% end

%% GUI

% Case name
caseName = 'DDM5PatchesOlympiadach';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Coupling method
method = 'Nitsche'; % 'Nitsche', 'Penalty', 'LagrangeMultipliers', 'AugmentedLagrangeMultipliers'
if ~strcmp(method,'Penalty') && ~strcmp(method,'LagrangeMultipliers') && ...
        ~strcmp(method,'AugmentedLagrangeMultipliers') && ~strcmp(method,'Nitsche')
    error('%s is not a valid method (Nitsche, Penalty, LagrangeMultipliers, AugmentedLagrangeMultipliers)',method);
end

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

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

% Loop over the patches
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

% Material matrices for the patches :
% ___________________________________

Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
     [1              parameters.nue 0
      parameters.nue 1              0
      0               0              (1-parameters.nue)/2];
Db = parameters.t^2/12*Dm;
for iPatches = 1:noPatches
    assignin('base',['Dm' num2str(iPatches)],Dm);
    assignin('base',['Db' num2str(iPatches)],Db);
end

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

% Co-simulation with Empire
isCosimulationWithEmpire = false;
strMatlabXml = 'undefined';

% Axis limits
limits = [-9.896059304703478 95.282059304703466 ...
          1.0e+02*0.579370000000000 1.0e+02 * 1.408920000000000 ...
          16.504349999999999  57.981849999999994];

%% Refinement

% General refinement variable
noStartHRef = 3;

% Define the refinement step
meshSize = 'coarse'; % 'coarse'
if strcmp(meshSize,'coarse')
    iPRef = 1;
    iHRef = 5;
elseif strcmp(meshSize,'fine')
    iPRef = 2;
    iHRef = 10;
else
    error('Select correct mesh size, coarse or fine');
end

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = iPRef - 1;
tp1 = a;
tq1 = a;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

% % Patch 2 :
% % _________
% 
% b = iPRef;
% tp2 = b;
% tq2 = b;
% [Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'');
% 
% % Patch 3 :
% % _________
% 
% c = iPRef - 1;
% tp3 = c;
% tq3 = c;
% [Xi3,Eta3,CP3,p3,q3] = degreeElevateBSplineSurface(p3,q3,Xi3,Eta3,CP3,tp3,tq3,'');
% 
% % Patch 4 :
% % _________
% 
% c = iPRef;
% tp4 = c;
% tq4 = c;
% [Xi4,Eta4,CP4,p4,q4] = degreeElevateBSplineSurface(p4,q4,Xi4,Eta4,CP4,tp4,tq4,'');
% 
% % Patch 5 :
% % _________
% 
% c = iPRef - 1;
% tp5 = c;
% tq5 = c;
% [Xi5,Eta5,CP5,p5,q5] = degreeElevateBSplineSurface(p5,q5,Xi5,Eta5,CP5,tp5,tq5,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = ceil((11/5)*((iHRef + noStartHRef) + 1));
noKnotsXi1 = a;
noKnotsEta1 = a;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface...
    (p1,Xi1,q1,Eta1,CP1,noKnotsXi1,noKnotsEta1,'');

% % Patch 2 :
% % _________
% 
% b = ceil((7/5)*((iHRef + noStartHRef) + 1)); % ceil((4/5)*(iHRef + 1))
% noKnotsXi2 = b;
% noKnotsEta2 = b;
% [Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface...
%     (p2,Xi2,q2,Eta2,CP2,noKnotsXi2,noKnotsEta2,'');
% 
% % Patch 3 :
% % _________
% 
% c = ceil((14/5)*((iHRef + noStartHRef) + 1));
% noKnotsXi3 = c;
% noKnotsEta3 = c;
% [Xi3,Eta3,CP3] = knotRefineUniformlyBSplineSurface...
%     (p3,Xi3,q3,Eta3,CP3,noKnotsXi3,noKnotsEta3,'');
% 
% % Patch 4 :
% % _________
% 
% c = ceil((3/5)*((iHRef + noStartHRef) + 1));
% noKnotsXi4 = c;
% noKnotsEta4 = c;
% [Xi4,Eta4,CP4] = knotRefineUniformlyBSplineSurface...
%     (p4,Xi4,q4,Eta4,CP4,noKnotsXi4,noKnotsEta4,'');
% 
% % Patch 5 :
% % _________
% 
% c = ceil((12/5)*((iHRef + noStartHRef) + 1));
% noKnotsXi5 = c;
% noKnotsEta5 = c;
% [Xi5,Eta5,CP5] = knotRefineUniformlyBSplineSurface...
%     (p5,Xi5,q5,Eta5,CP5,noKnotsXi5,noKnotsEta5,'');

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs1 = [];
xisup1 = {[0 0] [1 1] [1 1] [0 0]};
etasup1 = {[0 0] [0 0] [1 1] [1 1]};
if length(xisup1) ~= length(etasup1)
    error('Xi and eta extensions of the supports are not matching');
end
noSupp1 = length(xisup1);
for iSupp = 1:noSupp1
    for dir = 1:3
        homDOFs1 = findDofs3D(homDOFs1,xisup1{iSupp},etasup1{iSupp},dir,CP1);
    end
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC1.noCnd = 0;

% Embedded cables
cables1.No = 0;

% % Patch 2 :
% % _________
% 
% % Homogeneous Dirichlet boundary conditions
% homDOFs2 = [];
% xisup2 = {[0 0] [1 1] [0 1]};
% etasup2 = {[0 0] [0 0] [1 1]};
% if length(xisup2) ~= length(etasup2)
%     error('Xi and eta extensions of the supports are not matching');
% end
% noSupp2 = length(xisup2);
% for iSupp = 1:noSupp2
%     for dir = 1:3
%         homDOFs2 = findDofs3D(homDOFs2,xisup2{iSupp},etasup2{iSupp},dir,CP2);
%     end
% end
% 
% % Inhomogeneous Dirichlet boundary conditions
% inhomDOFs2 = [];
% valuesInhomDOFs2 = [];
% 
% % Weak homogeneous Dirichlet boundary conditions
% weakDBC2.noCnd = 0;
% 
% % Embedded cables
% cables2.No = 1;
% cables2.xiExtension = {[1 1]};
% cables2.etaExtension = {[0 1]};
% cables2.parameters = {parametersCable19027};
% cables2.int.type = 'default'; % 'user'
% cables2.int.noGPs = 16;
% 
% % Patch 3 :
% % _________
% 
% % Homogeneous Dirichlet boundary conditions
% homDOFs3 = [];
% xisup3 = {[0 0] [1 1] [1 1] [0 0]};
% etasup3 = {[0 0] [0 0] [1 1] [1 1]};
% if length(xisup3) ~= length(etasup3)
%     error('Xi and eta extensions of the supports are not matching');
% end
% noSupp3 = length(xisup3);
% for iSupp = 1:noSupp3
%     for dir = 1:3
%         homDOFs3 = findDofs3D(homDOFs3,xisup3{iSupp},etasup3{iSupp},dir,CP3);
%     end
% end
% 
% % Inhomogeneous Dirichlet boundary conditions
% inhomDOFs3 = [];
% valuesInhomDOFs3 = [];
% 
% % Weak homogeneous Dirichlet boundary conditions
% weakDBC3.noCnd = 0;
% 
% % Embedded cables
% cables3.No = 1;
% cables3.xiExtension = {[0 1]};
% cables3.etaExtension = {[0 0]};
% cables3.parameters = {parametersCable19026};
% cables3.int.type = 'default'; % 'user'
% cables3.int.noGPs = 16;
% 
% % Patch 4 :
% % _________
% 
% % Homogeneous Dirichlet boundary conditions
% homDOFs4 = [];
% xisup4 = {[0 0] [1 1] [0 1]};
% etasup4 = {[0 0] [0 0] [1 1]};
% if length(xisup4) ~= length(etasup4)
%     error('Xi and eta extensions of the supports are not matching');
% end
% noSupp4 = length(xisup4);
% for iSupp = 1:noSupp4
%     for dir = 1:3
%         homDOFs4 = findDofs3D(homDOFs4,xisup4{iSupp},etasup4{iSupp},dir,CP4);
%     end
% end
% 
% % Inhomogeneous Dirichlet boundary conditions
% inhomDOFs4 = [];
% valuesInhomDOFs4 = [];
% 
% % Weak homogeneous Dirichlet boundary conditions
% weakDBC4.noCnd = 0;
% 
% % Embedded cables
% cables4.No = 2;
% cables4.xiExtension = {[0 1] [0 0]};
% cables4.etaExtension = {[0 0] [0 1]};
% cables4.parameters = {parametersCable19021 parametersCable19022};
% cables4.int.type = 'default'; % 'user'
% cables4.int.noGPs = 16;
% 
% % Patch 5 :
% % _________
% 
% % Homogeneous Dirichlet boundary conditions
% homDOFs5 = [];
% xisup5 = {[0 0] [1 1] [1 1] [0 0]};
% etasup5 = {[0 0] [0 0] [1 1] [1 1]};
% if length(xisup5) ~= length(etasup5)
%     error('Xi and eta extensions of the supports are not matching');
% end
% noSupp5 = length(xisup5);
% for iSupp = 1:noSupp5
%     for dir = 1:3
%         homDOFs5 = findDofs3D(homDOFs5,xisup5{iSupp},etasup5{iSupp},dir,CP5);
%     end
% end
% 
% % Inhomogeneous Dirichlet boundary conditions
% inhomDOFs5 = [];
% valuesInhomDOFs5 = [];
% 
% % Weak homogeneous Dirichlet boundary conditions
% weakDBC5.noCnd = 0;
% 
% % Embedded cables
% cables5.No = 3;
% cables5.xiExtension = {[0 1] [1 1] [0 1]};
% cables5.etaExtension = {[0 0] [0 1] [1 1]};
% cables5.parameters = {parametersCable19023 parametersCable19024 parametersCable19025};
% cables5.int.type = 'default'; % 'user'
% cables5.int.noGPs = 16;

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

% FAmp2 = loadAmplitude;
% NBC2.noCnd = 1;
% xib2 = [0 1];   etab2 = [0 1];   dirForce2 = 'z';
% NBC2.xiLoadExtension = {xib2};
% NBC2.etaLoadExtension = {etab2};
% NBC2.loadAmplitude = {FAmp2};
% NBC2.loadDirection = {dirForce2};
% NBC2.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
% NBC2.isFollower = false;
% NBC2.isTimeDependent = true;
% 
% % Patch 3 :
% % _________
% 
% FAmp3 = loadAmplitude;
% NBC3.noCnd = 1;
% xib3 = [0 1];   etab3 = [0 1];   dirForce3 = 'z';
% NBC3.xiLoadExtension = {xib3};
% NBC3.etaLoadExtension = {etab3};
% NBC3.loadAmplitude = {FAmp3};
% NBC3.loadDirection = {dirForce3};
% NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
% NBC3.isFollower = false;
% NBC3.isTimeDependent = true;
% 
% % Patch 4 :
% % _________
% 
% FAmp4 = loadAmplitude;
% NBC4.noCnd = 1;
% xib4 = [0 1];   etab4 = [0 1];   dirForce4 = 'z';
% NBC4.xiLoadExtension = {xib4};
% NBC4.etaLoadExtension = {etab4};
% NBC4.loadAmplitude = {FAmp4};
% NBC4.loadDirection = {dirForce4};
% NBC4.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
% NBC4.isFollower = false;
% NBC4.isTimeDependent = true;
% 
% % Patch 5 :
% % _________
% 
% FAmp5 = loadAmplitude;
% NBC5.noCnd = 1;
% xib5 = [0 1];   etab5 = [0 1];   dirForce5 = 'z';
% NBC5.xiLoadExtension = {xib5};
% NBC5.etaLoadExtension = {etab5};
% NBC5.loadAmplitude = {FAmp5};
% NBC5.loadDirection = {dirForce5};
% NBC5.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
% NBC5.isFollower = false;
% NBC5.isTimeDependent = true;
% 
% % Collect all the Neumann boundary conditions into an arra<y
% NBC = {NBC1 NBC2 NBC3 NBC4 NBC5};
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Interface parametrizations    %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Patch 1 :
% % _________
% 
% % Connection with patch 2:
% xicoup12 = [0 1];   etacoup12 = [1 1];
% 
% % Collect all interfaces into arrays:
% xicoup1 = xicoup12;
% etacoup1 = etacoup12;
% 
% % Patch 2 :
% % _________
% 
% % Connection with patch 1:
% xicoup21 = [0 1];   etacoup21 = [0 0];
% 
% % Connection with patch 3:
% xicoup23 = [0 0];   etacoup23 = [0 1];
% 
% % Collect all interfaces into arrays:
% xicoup2 = [xicoup21 xicoup23];
% etacoup2 = [etacoup21 etacoup23];
% 
% % Patch 3 :
% % _________
% 
% % Connection with patch 2:
% xicoup32 = [1 1];   etacoup32 = [0 1];
% 
% % Connection with patch 4:
% xicoup34 = [0 0];   etacoup34 = [0 1];
% 
% % Connection with patch 4:
% xicoup35 = [0 1];   etacoup35 = [1 1];
% 
% % Collect all interfaces into arrays:
% xicoup3 = [xicoup32 xicoup34 xicoup35];
% etacoup3 = [etacoup32 etacoup34 etacoup35];
% 
% % Patch 4 :
% % _________
% 
% % Connection with patch 3:
% xicoup43 = [1 1];   etacoup43 = [0 1];
% 
% % Collect all interfaces into arrays:
% xicoup4 = xicoup43;
% etacoup4 = etacoup43;
% 
% % Patch 5 :
% % _________
% 
% % Connection with patch 3:
% xicoup53 = [0 0];   etacoup53 = [0 1];
% 
% % Collect all interfaces into arrays:
% xicoup5 = xicoup53;
% etacoup5 = etacoup53;
% 
% % Define connections :
% % ____________________
% 
% % Number of connections
% noConnections = 4;
% 
% % Define connections by patch numbers
% connections.No = noConnections;
% connections.xiEtaCoup = zeros(noConnections,10);
% connections.xiEtaCoup(:,:) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21
%                               2 3 xicoup23 etacoup23 xicoup32 etacoup32
%                               3 4 xicoup34 etacoup34 xicoup43 etacoup43
%                               3 5 xicoup35 etacoup35 xicoup53 etacoup53];
% connectionsLM = connections;
% connectionsALM = connections;

%% Fill up the arrays for the patches
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

%% Compute the radius of the hemisphere
p = BSplinePatches{1}.p;
p = BSplinePatches{1}.q;
Xi = BSplinePatches{1}.Xi;
Eta = BSplinePatches{1}.Eta;
CP = BSplinePatches{1}.CP;
isNURBS = BSplinePatches{1}.isNURBS;
xi = 0.5;
eta = 0.5;
xiSpan = findKnotSpan(xi,Xi,length(CP(:,1,1)));
etaSpan = findKnotSpan(eta,Eta,length(CP(1,:,1)));
R = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
radius = norm(computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R))

dHat = 'undefined';
connections.No = 0;
propSphere.radius = 4;
propSphere.center = [0, 0, 0]';
int.type = 'default';
figure(graph.index)
graph.index = plot_deviationFromSphereBSPlineSurface(BSplinePatches,dHat,connections,propSphere,int,graph,'outputEnabled');
graph.index = graph.index + 1;

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
% return;

%% Plot the multipatch geometry with the patch numbering
% color = [217 218 219]/255;
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,color,graph);

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

%% Plot the reference configuration for the multipatch geometry before the form-finding analysis
% color = [217 218 219]/255;
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');
% axis(limits);
% az = 30;
% el = 45;
% view(az,el);
% camlight(30,60);
% % camlight left; 
% lighting phong;

% axis off;
% title('');
% camlight left;
% lighting phong;

%% Form-finding properties
propFormFinding.tolerance = 1e-4;
propFormFinding.maxNoIter = 2e1;
propFormFinding.minNoIter = 2;

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

    % Assing the B-Spline patches for the Penalty method
    if exist('BSplinePatchesNitsche','var');
        clear BSplinePatchesNitsche;
    end

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
[BSplinePatchesFoFi,CPHistoryMultipatch,propCouplingNitsche,resHistory,...
    hasConverged,noIter] = ...
    solve_DDMFormFindingIGAMembrane(BSplinePatches,connections,...
    propCoupling,propFormFinding,solve_LinearSystem,'outputEnabled');
save(['data_FoFiDDMMiddleSailOlympiadach' '_' meshSize '_' method]);
return;

%% Compute the load vector for the visualization of the reference configuration
for iPatches = 1:noPatches
    BSplinePatchesFoFi{iPatches}.FGamma = ...
        zeros(3*BSplinePatchesFoFi{iPatches}.noCPs,1);
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

%% Load reference results from a classical Finite Element analysis
nodes = importdata('../../../preComputedData/isogeometricMembraneAnalysis/middleSailOlympiadach/referenceFEMSolution_andreasApostolatos/nodes');
elements = importdata('../../../preComputedData/isogeometricMembraneAnalysis/middleSailOlympiadach/referenceFEMSolution_andreasApostolatos/elements');
displacement = importdata('../../../preComputedData/isogeometricMembraneAnalysis/middleSailOlympiadach/referenceFEMSolution_andreasApostolatos/displacements');
mesh.elements = elements(:,2:4);
for iNodes = 1:length(nodes(:,1))
    [index,~] = find(displacement(iNodes,1) == nodes(:,1));
    nodes(index,2:4) = nodes(index,2:4) + displacement(iNodes,2:4);
end
mesh.nodes = nodes;

%% Plot the form-found geometry using classical Finite Elements
color = [217 218 219]/255;
labelsEnabled = false;
graph.index = plot_referenceConfigIGAMortarMapping...
    ({BSplinePatchesFoFi{1},BSplinePatchesFoFi{2},BSplinePatchesFoFi{5}},mesh,labelsEnabled,color,graph);
az = 30;
el = 45;
view([az el]);
camlight(20,40);
% camlight left;
lighting phong;

%% Compute the relative error from precomputed results using Finite Element analysis
propInt.type = 'default';
[relErrorL2GeoDomain,relErrorL2GeoInterface,graph.index] = ...
    computeDomainAndIterfaceErrorInL2AgainstFEMMembraneFormFinding...
    (BSplinePatchesFoFi,connections,mesh,propInt,graph,'outputEnabled');

%% Plot the reference configuration for the multipatch geometry after the form-finding analysis
color = [217 218 219]/255;
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatchesFoFi,connections,color,graph,'outputEnabled');
axis(limits);
az = 30;
el = 45;
view(az,el);
camlight(20,45);
% camlight left; 
lighting phong;

%% END