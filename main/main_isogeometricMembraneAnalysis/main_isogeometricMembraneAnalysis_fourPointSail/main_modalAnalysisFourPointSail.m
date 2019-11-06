%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Modal analysis for a form-found four-point sail, whose geometry is
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

% Include linear equation system solvers
addpath('../../../equationSystemSolvers/');

% Include transient analysis solvers
addpath('../../../transientAnalysis/');

% Include efficient computation functions
addpath('../../../efficientComputation/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../../CAGDKernel/CAGDKernel_graphics/',...
        '../../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the Isogeometric Kirchhoff-Love shell formulation
addpath('../../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../../isogeometricThinStructureAnalysis/graphicsMultipatches/',...
        '../../../isogeometricThinStructureAnalysis/loads/',...
        '../../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../../isogeometricThinStructureAnalysis/solvers/',...
        '../../../isogeometricThinStructureAnalysis/metrics/',...
        '../../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../../isogeometricThinStructureAnalysis/output/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/initialConditions/');

%% Load the NURBS parameters from a FoFi file
% fileName = 'FoFiFourPointSailp1q1Xi10Eta10';
% fileName = 'FoFiFourPointSailp2q2Xi10Eta10';
fileName = 'FoFiFourPointSailp3q3Xi50Eta50';
FoFiGeo = importdata(['./data_FoFiPool/' 'data_' fileName '.mat']);

% Polynomial orders
p = FoFiGeo.BSplinePatches{1}.p;
q = FoFiGeo.BSplinePatches{1}.q;

% Knot vectors
Xi = FoFiGeo.BSplinePatches{1}.Xi;
Eta = FoFiGeo.BSplinePatches{1}.Eta;

% Control Point coordinates and weights
CP = FoFiGeo.BSplinePatches{1}.CP;

% Flag on whether the basis is a B-Spline or a NURBS
isNURBS = FoFiGeo.BSplinePatches{1}.isNURBS;

%% Material constants

% Young's modulus
parameters.E = 8e+8;

% Poisson ratio
parameters.nue = .4;

% Thickness
parameters.t = 1e-3;

% Density (used only for dynamics)
parameters.rho = 800;

% Prestress for the membrane
sigma0 = 3e+3/parameters.t;
parameters.prestress.voigtVector = [sigma0
                                    sigma0
                                    0];

% Cable parameters
parametersCable.E = 1.6e+11;
parametersCable.radiusCS = 12e-3/2;
parametersCable.areaCS = pi*parametersCable.radiusCS^2;
parametersCable.rho = 8300;
parametersCable.prestress = 6e+4/parametersCable.areaCS;
cables.parameters = {parametersCable parametersCable parametersCable parametersCable};

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'fourPointSail';

% Define linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'default';
if strcmp(int.type,'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.xetaNGPForLoad = 6;
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
% .component: 'x', 'y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

% Function handle to writing out the results
% writeOutput = @writeResults4Carat;
writeOutput = @writeResults4GiD;
% writeOutput = 'undefined';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

%% Refinement

% Degree by which to elevate
a = 0;
tp = a;
tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 0;
refXi = scaling;
refEta = scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

% supports (Dirichlet boundary conditions)
homDOFs = [];
xiSup = [0 0];   etaSup = [0 0];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [0 0];   etaSup = [1 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 0];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [1 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC.noCnd = 0;
% weakDBC.method = 'nitsche';
weakDBC.method = 'penalty';
weakDBC.estimationStabilPrm = true;
weakDBC.alpha = 1e3;
weakDBC.xiExtension = {[0 0] [0 0] [1 1] [1 1]};
weakDBC.etaExtension = {[0 0] [1 1] [0 0] [1 1]};
weakDBC.int.type = 'default';
weakDBC.int.noGPs = 16;

% Embedded cables
cables.No = 4;
cables.xiExtension = {[0 0] [0 1] [0 1] [1 1]};
cables.etaExtension = {[0 1] [0 0] [1 1] [0 1]};
cables.int.type = 'default';
% cables.int.type = 'user';
cables.int.noGPs = 16;

% load (Neuman boundary conditions)
FAmp = 0;
xib = [0 1];   etab = [0 1];   dirForce = 'z';
NBC.noCnd = 1;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
windLoad = {@(x,y,z,t) FAmp*sin(2*pi()*x/Length/1e0)*sin(2*pi()*z/Height/1e0)};
NBC.loadAmplitude(1,1) = windLoad;
NBC.loadDirection = {dirForce};
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC.isFollower(1,1) = false;
NBC.isTimeDependent(1,1) = true;

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);
BSplinePatches = {BSplinePatch};

%% Output the initial geometry as a Carat++ input file
% pathToOutput = '../../outputGiD/isogeometricMembraneAnalysis/';
% writeOutMultipatchBSplineSurface4Carat(BSplinePatches,pathToOutput,caseName);

%% Compute the load vector for the visualization of the reference configuration
BSplinePatches{1}.FGamma = [];
for counterNBC = 1:NBC.noCnd
    BSplinePatches{1}.FGamma = zeros(3*BSplinePatches{1}.noCPs,1);
%     funcHandle = str2func(NBC.computeLoadVct{counterNBC});
%     BSplinePatches{1}.FGamma = ...
%         funcHandle(BSplinePatches{1}.FGamma,BSplinePatches{1},...
%         NBC.xiLoadExtension{counterNBC},...
%         NBC.etaLoadExtension{counterNBC},...
%         NBC.loadAmplitude{counterNBC},...
%         NBC.loadDirection{counterNBC},...
%         NBC.isFollower(counterNBC,1),...
%         0,int,'outputEnabled');
end

%% Plot reference configuration
color = [.85098 .8549 .85882];
connections = [];
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = 260;
el = 30;
view(az,el);
camlight(30,70);
lighting phong;
axis off;
title('');
    
%% Output the initial geometry to be read by GiD
pathToOutput = '../../../outputGiD/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4GiD(BSplinePatches,pathToOutput,caseName);

%% Write the geometry for Carat++
% writeOutMultipatchBSplineSurface4Carat(BSplinePatches,pathToOutput,caseName);

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-6;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 10000;
    
%% Perform a modal analysis over the linear isogeometric membrane system
noEig = 1000; % 2180
[eigenmodeShapes,naturalFrequencies,dHat] = solve_eigenmodeAnalysisIGAMembrane...
    (BSplinePatch,solve_LinearSystem,noEig,propNLinearAnalysis,'outputEnabled');
BSplinePatchesEq = BSplinePatches;
BSplinePatchesEq{1}.CP = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
    (BSplinePatches{1}.CP,dHat);
save(['data_modalAnalysis' caseName fileName]);
return; 

%% Plot the chosen eigenmode shapes
for i = 4
    scalingFct = 20;
    noEigen = i;
    graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
       (BSplinePatchesEq,scalingFct*eigenmodeShapes(:,noEigen),graph,'');
    az = 260;
    el = 30;
    view(az,el);
    camlight(30,70);
    lighting phong;
    axis off;
    title('');
end

%% END