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
% FoFiGeo = importdata('./data_FoFiSemisphereElmnts625P4Q4FF1e-4.mat');
% FoFiGeo = importdata('./data_FoFiSemisphereElmnts144P3Q3FF1e-4.mat');
% FoFiGeo = importdata('./data_FoFiSemisphereElmnts144P3Q3FF1e-4.mat');
% FoFiGeo = importdata('./data_FoFiSemisphereElmnts144P3Q3FF1e-10.mat');
FoFiGeo = importdata('./data_FoFiPool/data_FoFiSemisphereElmnts1024p4q4FF1e-11.mat');
% FoFiGeo = importdata('./data_FoFiSemisphere.mat');

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
parameters.E = 7.e5;

% Poisson ratio
parameters.nue = 0.49;

% thickness
parameters.t = 1.65e-4;

% density of the membrane
parameters.rho = 1050;

% Prestress for the membrane
% parameters.prestress.computeParametricCoordinates = @(X) [X(1,1)^2 + X(2,1)^2
%                                                           atan(X(2,1)/X(1,1))];
% parameters.prestress.computeBaseVectors = @(theta1,theta2) [ cos(theta2) -sin(theta2)
%                                                              sin(theta2) cos(theta2)
%                                                              0           0];
parameters.prestress.voigtVector = [7794.5
                                    7794.5
                                    0];

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'transientCircularMembrane';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Homogeneous Dirichlet boundary conditions
homDOFs = [];
xiSup = [0 0];   etaSup = [0 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [0 1];   etaSup = [0 0];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [0 1];   etaSup = [1 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = cell({});

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables.No = 0;

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
loadAmplitude2 = scaling*gravitationalAcceleration*parameters.rho*parameters.t; % -1.6995825
dirForce2 = 'z';
isFollower2 = false;

% load (Neuman boundary conditions)
xib = [0 1];   etab = [0 1];
NBC.noCnd = 2;
NBC.xiLoadExtension = {xib xib};
NBC.etaLoadExtension = {etab etab};
NBC.loadAmplitude = {loadAmplitude1 loadAmplitude2};
NBC.loadDirection = {dirForce1 dirForce2};
NBC.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
NBC.isFollower = [isFollower1; isFollower2];
NBC.isTimeDependent = [false; false];

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);
BSplinePatches = {BSplinePatch};

%% Output the initial geometry as a Carat++ input file
% pathToOutput = '../../outputGiD/isogeometricMembraneAnalysis/';
% writeOutMultipatchBSplineSurface4Carat(BSplinePatches,pathToOutput,caseName);

%% Compute the load vector for the visualization of the reference configuration
% BSplinePatches{1}.FGamma = [];
% for counterNBC = 1:NBC.noCnd
%     BSplinePatches{1}.FGamma = zeros(3*BSplinePatches{1}.noCPs,1);
% %     funcHandle = str2func(NBC.computeLoadVct{counterNBC});
% %     BSplinePatches{1}.FGamma = ...
% %         funcHandle(BSplinePatches{1}.FGamma,BSplinePatches{1},...
% %         NBC.xiLoadExtension{counterNBC},...
% %         NBC.etaLoadExtension{counterNBC},...
% %         NBC.loadAmplitude{counterNBC},...
% %         NBC.loadDirection{counterNBC},...
% %         NBC.isFollower(counterNBC,1),...
% %         0,int,'outputEnabled');
% end

%% Plot reference configuration
% color = [.85098 .8549 .85882];
% connections = [];
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');

%% Output the initial geometry to be read by GiD
% pathToOutput = '../../../outputGiD/isogeometricMembraneAnalysis/';
% writeOutMultipatchBSplineSurface4GiD(BSplinePatches,pathToOutput,caseName);

%% Write the geometry for Carat++
% writeOutMultipatchBSplineSurface4Carat(BSplinePatches,pathToOutput,caseName);

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-9;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 10000;
    
%% Perform a modal analysis over the linear isogeometric membrane system
noEig = 2180; % 40
[eigenmodeShapes,naturalFrequencies,dHat] = solve_eigenmodeAnalysisIGAMembrane...
    (BSplinePatch,solve_LinearSystem,noEig,propNLinearAnalysis,'outputEnabled');
BSplinePatchesEq = BSplinePatches;
BSplinePatchesEq{1}.CP = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
    (BSplinePatches{1}.CP,dHat);
save data_modalAnalysisSemisphereElmnts1024P4Q4FF1e-11;
return;

%% Plot the chosen eigenmode shapes

% Eigenfrequencies used for estimating the damping coefficients
freqI = 41.99;
% freqJ = 72.27;
freqJ = 57.5;

% Find the position of the eigenfrequencies in the frequency array
[~,idI] = min(abs(naturalFrequencies - freqI));
[~,idJ] = min(abs(naturalFrequencies - freqJ));
idJ = 1;

for i = idJ % [2 idI idJ]
    scalingFct = +0.001;
    noEigen = i;
    graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
       (BSplinePatchesEq,scalingFct*eigenmodeShapes(:,noEigen),graph,'');
    az = -310; % -40
    el = 40; % 40
    view(az,el);
    camlight(0,0);
    lighting phong;
    axis off;
    title(sprintf('f_{%d} = %d',idJ,naturalFrequencies(idJ,1)));
end

%% Plot the eigenfrequencies from both the modal analysis considering the tangent load from the follower loads and without considering the tangent load from the follower loads

% Load the data
data_ModalAnalysis = importdata('./data_ModalAnalysisPool/data_modalAnalysisSemisphereElmnts1024P4Q4FF1e-11.mat');
data_ModalAnalysis_withoutLoadTangent = importdata('./data_ModalAnalysisPool/data_modalAnalysisSemisphereElmnts1024P4Q4FF1e-11_noTangentFromLoad.mat');

% Plot the eigenfrequencies from both analyses
figure(graph.index)
plot(1:data_ModalAnalysis.noEig,data_ModalAnalysis.naturalFrequencies,1:data_ModalAnalysis_withoutLoadTangent.noEig,data_ModalAnalysis_withoutLoadTangent.naturalFrequencies)
grid on;
legend('Considering the tangent from the follower load','not considering the tangent from the follower load');
graph.index = graph.index + 1;

% Plot the ratio of the eigenfrequencies from both analysis
figure(graph.index)
plot(1:data_ModalAnalysis.noEig,data_ModalAnalysis.naturalFrequencies./data_ModalAnalysis_withoutLoadTangent.naturalFrequencies)
grid on;
graph.index = graph.index + 1;



%% END