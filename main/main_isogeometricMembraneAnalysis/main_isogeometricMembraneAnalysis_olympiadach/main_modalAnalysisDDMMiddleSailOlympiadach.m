%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Modal analysis for a form-found olympiadacg, whose geometry is
%        directly read from a file.
%
% Date : 4.12.2017
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

% Add all functions related to the low order basis functions
addpath('../../../basisFunctions/');

% Add all functions related to the classical Finite Element analysis
addpath('../../../FEMPlateInMembraneActionAnalysis/graphics/');

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
% fileName = 'FoFiDDMFourPointSailCoarse';
fileName = 'FoFiDDMMiddleSailOlympiadach';
method = 'Nitsche'; % Penalty, LagrangeMultipliers, AugmentedLagrangeMultipliers, Nitsche
meshSize = 'fine'; % coarse, fine
FoFiGeo = importdata(['./data_FoFiPool/' 'data_' fileName '_' meshSize '_' method '.mat']);

% Number of patches
noPatches = length(FoFiGeo.BSplinePatches);

% Patch 1 :
% _________

% Polynomial orders
p1 = FoFiGeo.BSplinePatchesFoFi{1}.p;
q1 = FoFiGeo.BSplinePatchesFoFi{1}.q;

% Knot vectors
Xi1 = FoFiGeo.BSplinePatchesFoFi{1}.Xi;
Eta1 = FoFiGeo.BSplinePatchesFoFi{1}.Eta;

% Control Point coordinates and weights
CP1 = FoFiGeo.BSplinePatchesFoFi{1}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS1 = FoFiGeo.BSplinePatchesFoFi{1}.isNURBS;

% Patch 2 :
% _________

% Polynomial orders
p2 = FoFiGeo.BSplinePatchesFoFi{2}.p;
q2 = FoFiGeo.BSplinePatchesFoFi{2}.q;

% Knot vectors
Xi2 = FoFiGeo.BSplinePatchesFoFi{2}.Xi;
Eta2 = FoFiGeo.BSplinePatchesFoFi{2}.Eta;

% Control Point coordinates and weights
CP2 = FoFiGeo.BSplinePatchesFoFi{2}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS2 = FoFiGeo.BSplinePatchesFoFi{2}.isNURBS;

% Patch 3 :
% _________

% Polynomial orders
p3 = FoFiGeo.BSplinePatchesFoFi{3}.p;
q3 = FoFiGeo.BSplinePatchesFoFi{3}.q;

% Knot vectors
Xi3 = FoFiGeo.BSplinePatchesFoFi{3}.Xi;
Eta3 = FoFiGeo.BSplinePatchesFoFi{3}.Eta;

% Control Point coordinates and weights
CP3 = FoFiGeo.BSplinePatchesFoFi{3}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS3 = FoFiGeo.BSplinePatchesFoFi{3}.isNURBS;

% Patch 4 :
% _________

% Polynomial orders
p4 = FoFiGeo.BSplinePatchesFoFi{4}.p;
q4 = FoFiGeo.BSplinePatchesFoFi{4}.q;

% Knot vectors
Xi4 = FoFiGeo.BSplinePatchesFoFi{4}.Xi;
Eta4 = FoFiGeo.BSplinePatchesFoFi{4}.Eta;

% Control Point coordinates and weights
CP4 = FoFiGeo.BSplinePatchesFoFi{4}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS4 = FoFiGeo.BSplinePatchesFoFi{4}.isNURBS;

% Patch 5 :
% _________

% Polynomial orders
p5 = FoFiGeo.BSplinePatchesFoFi{5}.p;
q5 = FoFiGeo.BSplinePatchesFoFi{5}.q;

% Knot vectors
Xi5 = FoFiGeo.BSplinePatchesFoFi{5}.Xi;
Eta5 = FoFiGeo.BSplinePatchesFoFi{5}.Eta;

% Control Point coordinates and weights
CP5 = FoFiGeo.BSplinePatchesFoFi{5}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS5 = FoFiGeo.BSplinePatchesFoFi{5}.isNURBS;

%% Material constants

% General parameters
n0 = 4000000;
nTilde0 = 400000;

% Material properties for the plexiglas (assuming not contributing to the 
% stiffness of the membrane)
propPlexiglas.EYoung = 0;
propPlexiglas.nue = 0;
propPlexiglas.rho = 5e2;
propPlexiglas.thickness = 7e-3;

% Material properties for the inner cables of the roof
propSteelCables.EYoung = 210e9;
propSteelCables.nue = 0.3;
propSteelCables.rho = 7.85e3;
propSteelCables.diameter = 14e-3;
propSteelCables.area = pi*(propSteelCables.diameter/2)^2;

% Material properties of the membrane model
propMembrane.nue = 0.0;

% Characteristic length of homogenization
characteristicLength = 75e-2;

% Number of cables per characteristic element
noCables = 4;

% Compute the homogenized material properties of the membrane and the
% boundary cables
[propMembrane,propCables] = computeHomogenizedMaterialPropertiesOlympiadach...
    (propPlexiglas,propSteelCables,propMembrane,characteristicLength,noCables);

% General parameters
% EYoung = 2.1e11;
% thickness = 5e-3;
% diameterCables = 14e-3;
% radiusCables = diameterCables/2;
% areaCables = pi*radiusCables^2;
% areaCablesBoundary = 5e-2;
% nue = 4e-1;
% n0 = 400000;
% nTilde0 = 40000;
% scalingYoungModulusCables = 1e2;
% rhoPlexiglas = 5e2;
% rhoSteel = 7.85e3;

% Homogenisation of the membrane density
% a75 = 75e-2;
% noCablesPerSquare75 = 4;
% massPerSquare75 = rhoPlexiglas*thickness*a75^2 + noCablesPerSquare75*rhoSteel*areaCables*a75;
% massPerSquare75 = massPerSquare75 + 0.15*massPerSquare75;
% rhoMembrane = massPerSquare75/a75^2/thickness;

% Membrane :
% ----------

parameters.E = propMembrane.EYoung;
parameters.nue = propMembrane.nue;
parameters.t = propMembrane.thickness;
parameters.rho = propMembrane.rho;
parameters.prestress.voigtVector = [50*n0
                                    50*n0
                                    0];

% Cables :
% --------
                                
n19021 = nTilde0*1680;
n19022 = nTilde0*1730;
n19023 = nTilde0*1595;
n19024 = nTilde0*1690;
n19025 = nTilde0*1350;
n19026 = nTilde0*1250;
n19027 = nTilde0*1250;
n19028 = nTilde0*1575;
n19029 = nTilde0*1550;
n19030 = nTilde0*1845;
constNum = 19020;
noCables = 10;
for iCables = 1:noCables
    evalin('base',['parametersCable' num2str(constNum + iCables) '.E = propCables.EYoung;']);
    evalin('base',['parametersCable' num2str(constNum + iCables) '.areaCS = propCables.area;']);
    evalin('base',['parametersCable' num2str(constNum + iCables) '.rho = propCables.rho;']);
    evalin('base',['parametersCable' num2str(constNum + iCables) '.prestress = n' num2str(constNum + iCables) ';']);
end

%% GUI

% Case name
caseName = 'DDM5PatchesMiddleSailOlympiadach';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Coupling method
method = 'AugmentedLagrangeMultipliers'; % 'Nitsche', 'Penalty', 'LagrangeMultipliers', 'AugmentedLagrangeMultipliers'
if ~strcmp(method,'Penalty') && ~strcmp(method,'LagrangeMultipliers') && ...
        ~strcmp(method,'AugmentedLagrangeMultipliers') && ~strcmp(method,'Nitsche')
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

Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
     [1              parameters.nue 0
      parameters.nue 1              0
      0               0              (1-parameters.nue)/2];
Db = parameters.t^2/12*Dm;
for iPatches = 1:noPatches
    assignin('base',['Dm' num2str(iPatches)],Dm);
    assignin('base',['Db' num2str(iPatches)],Db);
end

% Assign the penalty factors

% Function handle to writing out the results
% writeOutput = @writeResults4Carat;
% writeOutput = @writeResults4GiD;
writeOutput = 'undefined';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

% Axis limits
limits = [-9.896059304703478 95.282059304703466 ...
          1.0e+02*0.579370000000000 1.0e+02 * 1.408920000000000 ...
          16.504349999999999  57.981849999999994];

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
cables1.No = 3;
cables1.xiExtension = {[0 0] [0 1] [1 1]};
cables1.etaExtension = {[0 1] [0 0] [0 1]};
cables1.parameters = {parametersCable19028 parametersCable19029 parametersCable19030};
cables1.int.type = 'default'; % 'user'
cables1.int.noGPs = 16;

% Patch 2 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs2 = [];
xisup2 = {[0 0] [1 1] [0 1]};
etasup2 = {[0 0] [0 0] [1 1]};
if length(xisup2) ~= length(etasup2)
    error('Xi and eta extensions of the supports are not matching');
end
noSupp2 = length(xisup2);
for iSupp = 1:noSupp2
    for dir = 1:3
        homDOFs2 = findDofs3D(homDOFs2,xisup2{iSupp},etasup2{iSupp},dir,CP2);
    end
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs2 = [];
valuesInhomDOFs2 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC2.noCnd = 0;

% Embedded cables
cables2.No = 1;
cables2.xiExtension = {[1 1]};
cables2.etaExtension = {[0 1]};
cables2.parameters = {parametersCable19027};
cables2.int.type = 'default'; % 'user'
cables2.int.noGPs = 16;

% Patch 3 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs3 = [];
xisup3 = {[0 0] [1 1] [1 1] [0 0]};
etasup3 = {[0 0] [0 0] [1 1] [1 1]};
if length(xisup3) ~= length(etasup3)
    error('Xi and eta extensions of the supports are not matching');
end
noSupp3 = length(xisup3);
for iSupp = 1:noSupp3
    for dir = 1:3
        homDOFs3 = findDofs3D(homDOFs3,xisup3{iSupp},etasup3{iSupp},dir,CP3);
    end
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs3 = [];
valuesInhomDOFs3 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC3.noCnd = 0;

% Embedded cables
cables3.No = 1;
cables3.xiExtension = {[0 1]};
cables3.etaExtension = {[0 0]};
cables3.parameters = {parametersCable19026};
cables3.int.type = 'default'; % 'user'
cables3.int.noGPs = 16;

% Patch 4 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs4 = [];
xisup4 = {[0 0] [1 1] [0 1]};
etasup4 = {[0 0] [0 0] [1 1]};
if length(xisup4) ~= length(etasup4)
    error('Xi and eta extensions of the supports are not matching');
end
noSupp4 = length(xisup4);
for iSupp = 1:noSupp4
    for dir = 1:3
        homDOFs4 = findDofs3D(homDOFs4,xisup4{iSupp},etasup4{iSupp},dir,CP4);
    end
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs4 = [];
valuesInhomDOFs4 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC4.noCnd = 0;

% Embedded cables
cables4.No = 2;
cables4.xiExtension = {[0 1] [0 0]};
cables4.etaExtension = {[0 0] [0 1]};
cables4.parameters = {parametersCable19021 parametersCable19022};
cables4.int.type = 'default'; % 'user'
cables4.int.noGPs = 16;

% Patch 5 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs5 = [];
xisup5 = {[0 0] [1 1] [1 1] [0 0]};
etasup5 = {[0 0] [0 0] [1 1] [1 1]};
if length(xisup5) ~= length(etasup5)
    error('Xi and eta extensions of the supports are not matching');
end
noSupp5 = length(xisup5);
for iSupp = 1:noSupp5
    for dir = 1:3
        homDOFs5 = findDofs3D(homDOFs5,xisup5{iSupp},etasup5{iSupp},dir,CP5);
    end
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs5 = [];
valuesInhomDOFs5 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC5.noCnd = 0;

% Embedded cables
cables5.No = 3;
cables5.xiExtension = {[0 1] [1 1] [0 1]};
cables5.etaExtension = {[0 0] [0 1] [1 1]};
cables5.parameters = {parametersCable19023 parametersCable19024 parametersCable19025};
cables5.int.type = 'default'; % 'user'
cables5.int.noGPs = 16;

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

FAmp2 = loadAmplitude;
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

FAmp3 = loadAmplitude;
NBC3.noCnd = 1;
xib3 = [0 1];   etab3 = [0 1];   dirForce3 = 'z';
NBC3.xiLoadExtension = {xib3};
NBC3.etaLoadExtension = {etab3};
NBC3.loadAmplitude = {FAmp3};
NBC3.loadDirection = {dirForce3};
NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC3.isFollower = false;
NBC3.isTimeDependent = true;

% Patch 4 :
% _________

FAmp4 = loadAmplitude;
NBC4.noCnd = 1;
xib4 = [0 1];   etab4 = [0 1];   dirForce4 = 'z';
NBC4.xiLoadExtension = {xib4};
NBC4.etaLoadExtension = {etab4};
NBC4.loadAmplitude = {FAmp4};
NBC4.loadDirection = {dirForce4};
NBC4.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC4.isFollower = false;
NBC4.isTimeDependent = true;

% Patch 5 :
% _________

FAmp5 = loadAmplitude;
NBC5.noCnd = 1;
xib5 = [0 1];   etab5 = [0 1];   dirForce5 = 'z';
NBC5.xiLoadExtension = {xib5};
NBC5.etaLoadExtension = {etab5};
NBC5.loadAmplitude = {FAmp5};
NBC5.loadDirection = {dirForce5};
NBC5.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC5.isFollower = false;
NBC5.isTimeDependent = true;

% Collect all the Neumann boundary conditions into an array
NBC = {NBC1 NBC2 NBC3 NBC4 NBC5};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface parametrizations     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Connection with patch 2:
xicoup12 = [0 1];   etacoup12 = [1 1];

% Collect all interfaces into arrays:
xicoup1 = xicoup12;
etacoup1 = etacoup12;

% Patch 2 :
% _________

% Connection with patch 1:
xicoup21 = [0 1];   etacoup21 = [0 0];

% Connection with patch 3:
xicoup23 = [0 0];   etacoup23 = [0 1];

% Collect all interfaces into arrays:
xicoup2 = [xicoup21 xicoup23];
etacoup2 = [etacoup21 etacoup23];

% Patch 3 :
% _________

% Connection with patch 2:
xicoup32 = [1 1];   etacoup32 = [0 1];

% Connection with patch 4:
xicoup34 = [0 0];   etacoup34 = [0 1];

% Connection with patch 4:
xicoup35 = [0 1];   etacoup35 = [1 1];

% Collect all interfaces into arrays:
xicoup3 = [xicoup32 xicoup34 xicoup35];
etacoup3 = [etacoup32 etacoup34 etacoup35];

% Patch 4 :
% _________

% Connection with patch 3:
xicoup43 = [1 1];   etacoup43 = [0 1];

% Collect all interfaces into arrays:
xicoup4 = xicoup43;
etacoup4 = etacoup43;

% Patch 5 :
% _________

% Connection with patch 3:
xicoup53 = [0 0];   etacoup53 = [0 1];

% Collect all interfaces into arrays:
xicoup5 = xicoup53;
etacoup5 = etacoup53;

% Define connections :
% ____________________

% Number of connections
noConnections = 4;

% Define connections by patch numbers
connections.No = noConnections;
connections.xiEtaCoup = zeros(noConnections,10);
connections.xiEtaCoup(:,:) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21
                              2 3 xicoup23 etacoup23 xicoup32 etacoup32
                              3 4 xicoup34 etacoup34 xicoup43 etacoup43
                              3 5 xicoup35 etacoup35 xicoup53 etacoup53];
connectionsLM = connections;
connectionsALM = connections;

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
if strcmp(method,'LagrangeMultipliers') || strcmp(method,'AugmentedLagrangeMultipliers')
    fprintf('**************************************************\n');
    fprintf('* Creating interface Lagrange Multipliers fields *\n');
    fprintf('**************************************************\n\n');
    for iConnections = 1:connections.No
        %% Get the IDs of the patches involved
        idI = connections.xiEtaCoup(iConnections,1);
        idJ = connections.xiEtaCoup(iConnections,2);
        fprintf(['\t' 'Coupling between patches %d and %d \n'],idI,idJ);
        fprintf(['\t' '---------------------------------- \n\n']);

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
            fprintf(['\t' 'The coupling interface of patch %d is along xi\n'],idI);
        elseif xiCoupI(1,1) == xiCoupI(1,2) && etaCoupI(1,1) ~= etaCoupI(1,2)
            isOnXiI = false;
            fprintf(['\t' 'The coupling interface of patch %d is along eta\n'],idI);
        else
            error('Either the xi or the eta direction has to be fixed along the interface for patch %d',idI);
        end
        xiCoupJ = connections.xiEtaCoup(iConnections,7:8);
        etaCoupJ = connections.xiEtaCoup(iConnections,9:10);
        if xiCoupJ(1,1) ~= xiCoupJ(1,2) && etaCoupJ(1,1) == etaCoupJ(1,2)
            isOnXiJ = true;
            fprintf(['\t' 'The coupling interface of patch %d is along xi\n'],idJ);
        elseif xiCoupJ(1,1) == xiCoupJ(1,2) && etaCoupJ(1,1) ~= etaCoupJ(1,2)
            isOnXiJ = false;
            fprintf(['\t' 'The coupling interface of patch %d is along eta\n'],idJ);
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
        fprintf(['\t' 'Degree elevating the Lagrange Multipliers field to %d\n'],pLM);
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
        noLambda = ceil(max([noKnotsI noKnotsJ])*scaleLM);
        fprintf(['\t' 'Uniformly inserting %d knots in the Lagrange Multipliers field\n'],noLambda);
        [XiLambdaLM,CPLambdaLM] = knotRefineUniformlyBSplineCurve...
            (noLambda,pLambda,XiLambda,CPLambda,'');
        scaleALM = 0.8;
        noLambda = ceil(max([noKnotsI noKnotsJ])*scaleALM);
        fprintf(['\t' 'Uniformly inserting %d knots in the augmented Lagrange Multipliers field\n'],noLambda);
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

%% Nonlinear analysis properties
propNLinearAnalysis.method = 'newtonRapshon';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 1e-2;
propNLinearAnalysis.maxIter = 35;

%% Perform modal analysis and visualize the chosen eigenmode shape

% Solve the eigenvalue problem
noEig = 200; % 400
[eigenmodeShapes,naturalFrequencies,BSplinePatches,dHat] = ...
    solve_DDMEigenmodeAnalysisIGAMembrane...
    (BSplinePatches,connections,propCoupling,solve_LinearSystem,noEig,...
    propNLinearAnalysis,'outputEnabled');
save(['./data_ModalAnalysisPool/' 'data_ModalAnalysis' '_' meshSize '_' method]);
return;

%% Load a reference solution for the modal analysis
mesh.nodes = importdata('../../../preComputedData/isogeometricMembraneAnalysis/middleSailOlympiadach/referenceFEMSolution_ModalAnalysis/nodes');
elements = importdata('../../../preComputedData/isogeometricMembraneAnalysis/middleSailOlympiadach/referenceFEMSolution_ModalAnalysis/elements');
mesh.elements = elements(:,1:3);
fileNameFEM = '../../../preComputedData/isogeometricMembraneAnalysis/middleSailOlympiadach/referenceFEMSolution_ModalAnalysis/displacementsModeShapes';
naturalFrequenciesFEM = importdata('../../../preComputedData/isogeometricMembraneAnalysis/middleSailOlympiadach/referenceFEMSolution_ModalAnalysis/eigenfrequencies');
naturalFrequenciesFEM = naturalFrequenciesFEM.data(:,2);
fstring = fileread(fileNameFEM);
idEig = 7; % 1, 4, 7
stringSeparator = ['Result "Eigenvector" "Eigenvector nmb." ' num2str(idEig) ' Vector OnNodes'];
block = regexp(fstring,stringSeparator,'split');
block = regexp(block{2},'Values','split');
block(1) = [];
block(2:end) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
displacements = cell2mat(out);

%% Visualize the chosen reference modal shape
scaling = 1.0;
displacementsOrdered = zeros(3*length(displacements(:,1)),1);
for iDisp = 1:length(displacements(:,1))
%     [index,~] = find(displacements(iNodes,1) == mesh.nodes(:,1));
    displacementsOrdered(3*iDisp - 2,1) = scaling*displacements(iDisp,2);
    displacementsOrdered(3*iDisp - 1,1) = scaling*displacements(iDisp,3);
    displacementsOrdered(3*iDisp,1) = scaling*displacements(iDisp,4);
end
resultant = 'displacement';
component = '2norm';
graph.visualization.geometry = 'reference';
mesh.nodes = mesh.nodes(:,2:4);
graph.index = plot_currentConfigurationAndResultants...
    (mesh,[],displacementsOrdered,[],[],resultant,component,graph);
colormap('jet');
h = colorbar;
limitsColorbar = get(h,'Limits');
axis(limits);
az = 30;
el = 45;
view(az,el);

%% Plot the selected mode shape when using classical Finite Elements
% color = [217 218 219]/255;
% labelsEnabled = false;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     ({},mesh,labelsEnabled,color,graph);
% az = 30;
% el = 45;
% view([az el]);
% camlight(20,40);
% % camlight left;
% lighting phong;

%% Plot the error in the numerical spectra
% omegaFEM = (2*pi).*naturalFrequenciesFEM;
% omega = (2*pi).*naturalFrequencies;
% omegaRatio = zeros(noEig,1);
% for iFreq = 1:noEig
%     omegaRatio(iFreq,1) = omega(iFreq,1)/omegaFEM(iFreq,1);
% end
% xAxis = 1:noEig;
% plot(xAxis,omegaRatio);
% % return;

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

% Visualization of the eigenmode shapes
for i = 6 % idEig % 1, 5, 6 % idEig
    scalingFct = 1.0;
    noEigen = i;
    graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
        (BSplinePatchesEq,scalingFct*eigenmodeShapes(EFTPatchesWTLM,noEigen),graph,'outputEnabled');
    colormap('jet');
    colorbar;
%     caxis(limitsColorbar);
    axis(limits);
    az = 30;
    el = 45;
    view(az,el);
%     lighting phong;
%     axis off;
%     title('');
end

%% END
