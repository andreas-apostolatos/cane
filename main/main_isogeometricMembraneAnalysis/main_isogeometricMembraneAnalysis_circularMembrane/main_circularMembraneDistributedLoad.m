%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : A circular membrane is subject to distributed load.
%
% Date : 18.06.2015
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
        '../../../isogeometricThinStructureAnalysis/output/',...
        '../../../isogeometricThinStructureAnalysis/loads/',...
        '../../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../../isogeometricThinStructureAnalysis/solvers/',...
        '../../../isogeometricThinStructureAnalysis/metrics/',...
        '../../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/');

%% NURBS parameters

% Global variables:
Radius = .5;
edgeSquare = Radius*sqrt(2);

% Polynomial degrees
p = 2;
q = 2;

% Knot vectors
Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [-edgeSquare/2 0 edgeSquare/2
             -edgeSquare   0 edgeSquare
             -edgeSquare/2 0 edgeSquare/2];
         
% y-coordinates
CP(:,:,2) = [-edgeSquare/2 -edgeSquare   -edgeSquare/2
             0             0             0
             edgeSquare/2  edgeSquare    edgeSquare/2];
         
% z-coordinates
CP(:,:,3) = [0 0 0
             0 0 0
             0 0 0];
       
% Weights
weight = sqrt(2)/2;
CP(:,:,4) = [1      weight 1
             weight 1      weight
             1      weight 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i= 1:nxi
    for j=1:neta
        if CP(i,j,4)~=1
            isNURBS = 1;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% Material constants

% Young's modulus
parameters.E = 1.23e6;

% Poisson ratio
parameters.nue = 0.49;

% Thickness of the membrane
parameters.t = 1e-4;

% Density of the shell (used only for dynamics)
parameters.rho = 1050;

% Pressure Load amplitude
scaling = 1; % 1
FAmp = - 1e4*scaling;

% Prestress for the membrane
% parameters.prestress.computeParametricCoordinates = @(X) [X(1,1)^2 + X(2,1)^2
%                                                           atan(X(2,1)/X(1,1))];
% parameters.prestress.computeBaseVectors = @(theta1,theta2) [cos(theta2) -sin(theta2)
%                                                      sin(theta2) cos(theta2)
%                                                      0           0];
% scaling = 1;
parameters.prestress.voigtVector = [abs(FAmp)*Radius/2/parameters.t*scaling
                                    abs(FAmp)*Radius/2/parameters.t*scaling
                                    0];
% parameters.prestress.voigtVector = [250000*parameters.t
%                                     250000*parameters.t
%                                     0];

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'circularPlate';

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

% Linear system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','curvature','force','moment','shearForce'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x', 'y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

%% Refinement

% Degree by which to elevate
a = 0;
tp = a;
tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 5; % 5
refXi = scaling;
refEta = scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

% supports (Dirichlet boundary conditions)

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

% Weak boundary conditions
weakDBC.noCnd = 0;
weakDBC.method = 'Nitsche';
weakDBC.estimationStabilPrm = true;
weakDBC.alpha = parameters.E/0.003858024691358/1e-2;
weakDBC.xiExtension = {[0 0] [0 1] [1 1] [0 1]};
weakDBC.etaExtension = {[0 1] [0 0] [0 1] [1 1]};
% weakDBC.computeConstMtx = 'computeWeakDBCMtxPenaltyIGAMembrane';
weakDBC.computeTangMtxResVct = 'computeWeakDBCTangMtxResVctNitscheIGAMembrane';
weakDBC.int.type = 'default';
weakDBC.int.noGPs = 16;

% Embedded cables
cables.No = 0;

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% load (Neuman boundary conditions)
xib = [0 1];   etab = [0 1];   dirForce = 'normal';
NBC.noCnd = 1;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC.isFollower(1,1) = true;
NBC.isTimeDependent(1,1) = false;

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,valuesInhomDOFs,...
    weakDBC,cables,NBC,[],[],[],[],[],int);

%% Plot the distribution of the determinant of the geometrical Jacobian
% figure(graph.index)
% plot_postprocBSplineSurfaceGeometricJacobian(p,q,Xi,Eta,CP,isNURBS);
% shading interp;
% colormap('jet');
% camlight left;
% lighting phong;
% colorbar;
% axis equal;
% graph.index = graph.index + 1;

%% Perform a form finding analysis to find the shape of the membrane
% propFormFinding.tolerance = 1e-7;
% propFormFinding.maxNoIter = 1000;
% [BSplinePatch,CPHistory,resHistory,hasConverged,noIter] = ...
%     solve_formFindingIGAMembrane(BSplinePatch,propFormFinding,...
%     solve_LinearSystem,'outputEnabled');
% save data_cicularPlateFormFoundIntoSemisphereSinglePatch;
% return;

%% Compute the load vector for the visualization of the reference configuration
FGamma = [];
for counterNBC = 1:BSplinePatch.NBC.noCnd
    funcHandle = str2func(BSplinePatch.NBC.computeLoadVct{counterNBC});
    FGamma = funcHandle(FGamma,BSplinePatch,BSplinePatch.NBC.xiLoadExtension{counterNBC},...
                    BSplinePatch.NBC.etaLoadExtension{counterNBC},...
                    BSplinePatch.NBC.loadAmplitude{counterNBC},...
                    BSplinePatch.NBC.loadDirection{counterNBC},...
                    BSplinePatch.NBC.isFollower(counterNBC,1),0,...
                    BSplinePatch.int,'outputEnabled');
end

%% Plot reference configuration
figure(graph.index)
plot_referenceConfigurationIGAThinStructure(BSplinePatch.p,BSplinePatch.q,...
    BSplinePatch.Xi,BSplinePatch.Eta,BSplinePatch.CP,BSplinePatch.isNURBS,...
    BSplinePatch.homDOFs,FGamma,'outputEnabled');
title('Reference configuration for an isogeometric Kirchhoff-Love shell');
graph.index = graph.index + 1;

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-9;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 200;

%% Solve the steady-state problem
plot_IGANLinear = '';
[dHat,CPHistory,resHistory,hasConverged,BSplinePatches,minElASize] = ...
    solve_IGAMembraneNLinear(BSplinePatch,propNLinearAnalysis,...
    solve_LinearSystem,plot_IGANLinear,graph,'outputEnabled');
scaling = 1;
graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
    (BSplinePatches,scaling*dHat,graph,'outputEnabled');
return;

%% Output the initial geometry as a Carat++ input file
pathToOutput = '../../../inputCarat/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4Carat(BSplinePatches,pathToOutput,caseName);

%% Solve the system applying nonlinear analysis
% plot_IGANLinear = @plot_postprocIGAMembraneMultipatchesNLinear;
plot_IGANLinear = '';
[dHatNLinear,CPHistory,resHistory,hasConverged,minElArea] = solve_IGAMembraneNLinear...
    (BSplinePatch,propNLinearAnalysis,solve_LinearSystem,plot_IGANLinear,...
    graph,'outputEnabled');

%% Postprocessing

% Visualize the deformed configuration
graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
    (BSplinePatches,dHatNLinear,graph,'outputEnabled');

% Compute the current position of the middle point of the patch
xi = .5;
eta = .5;
xiSpan = findKnotSpan(xi,Xi,length(CP(:,1,1)));
etaSpan = findKnotSpan(eta,Eta,length(CP(1,:,1)));
R = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CPHistory{1}(:,:,:,end),isNURBS,0);
X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CPHistory{1}(:,:,:,end),R)

%% End