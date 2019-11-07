%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Nonlinear analysis over a single patch semisphere which's shape is
%        obtained using form-finding of a circular plate under uniform
%        pressure load.
%
% Date : 09.11.2016
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
% FoFiGeo = importdata('./data_FoFiPool/data_FoFiSemisphereElmnts1024p3q3FF1e-9.mat');
% FoFiGeo = importdata('./data_FoFiPool/data_FoFiSemisphere_Tol_1e-09_noElmnts_144_p_2_q_2.mat');
% FoFiGeo = importdata('./data_FoFiPool/data_FoFiSemisphere_Tol_1e-09_noElmnts_625_p_2_q_2.mat');
FoFiGeo = importdata('./data_FoFiPool/data_FoFiSemisphere_Tol_1e-09_noElmnts_625_p_3_q_3.mat');

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

% Thickness of the shell
parameters.t = 1.65e-4;

% Density of the shell (used only for dynamics)
parameters.rho = 1050;

% Pressure Load amplitude and weight
% FAmpPressure = - 43.2;
FAmpPressure = - 1e9;
scaling = 0.0;
FAmpWeight = - 9.81*parameters.rho*parameters.t*scaling;
scaling = 0.0;
FAmpPoint = 0.032*scaling;

% Prestress for the membrane
% parameters.prestress.computeParametricCoordinates = @(X) [X(1,1)^2 + X(2,1)^2
%                                                           atan(X(2,1)/X(1,1))];
% parameters.prestress.computeBaseVectors = @(theta1,theta2) [ cos(theta2) -sin(theta2)
%                                                              sin(theta2) cos(theta2)
%                                                              0           0];
% parameters.prestress.voigtVector = [7794.5
%                                     7794.5
%                                     0];
Radius = norm(squeeze(CP(end,1,1:3)));
parameters.prestress.voigtVector = [abs(FAmpPressure)*Radius/2/parameters.t
                                    abs(FAmpPressure)*Radius/2/parameters.t
                                    0];
                                
% Define the material matrices
Dm = parameters.E*parameters.t/(1 - parameters.nue^2)*...
      [1              parameters.nue 0
       parameters.nue 1              0
       0              0              (1 - parameters.nue)/2];
Db = parameters.t^2/12*Dm;

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'transientSemisphere';

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
graph.postprocConfig = 'referenceCurrent';

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

%% Get the parametric location of the point load application
xi0 = .5;
eta0 = .5;
P = Radius/sqrt(2)*[1; 1; 1];
propNewtonRaphson.eps = 1e-15;
propNewtonRaphson.maxIt = 1000;
[xiP,etaP,PProjected,isProjected,noIter] = ...
    computeNearestPointProjectionOnBSplineSurface...
    (P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,propNewtonRaphson);
if ~isProjected
    error('The point of the load application could not be found on patch 2');
end

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

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;
weakDBC.xiExtension = {[0 0] [0 1] [1 1] [0 1]};
weakDBC.etaExtension = {[0 1] [0 0] [0 1] [1 1]};
weakDBC.int.type = 'default';
weakDBC.int.noGPs = 16;
weakDBC.method = 'penalty';

% Embedded cables
cables.No = 0;

% load (Neuman boundary conditions)
% FAmp = 0;
xib1 = [0 1];   etab1 = [0 1];   dirForcePressure = 'normal';
xib2 = [0 1];   etab2 = [0 1];   dirForceWeight = 'z';
xibP2 = xiP;    etaP2 = etaP;    dirForcePoint = 'normal';
NBC.noCnd = 3;
NBC.xiLoadExtension = {xib1 xib2 xibP2};
NBC.etaLoadExtension = {etab1 etab2 etaP2};
NBC.loadAmplitude = {FAmpPressure FAmpWeight FAmpPoint};
NBC.loadDirection = {dirForcePressure dirForceWeight dirForcePoint};
NBC.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure' 'computeLoadVctPointIGAThinStructure'};
NBC.isFollower = [true; false; false];
NBC.isTimeDependent = [false; false; false];

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);
BSplinePatches = {BSplinePatch};
noPatches = length(BSplinePatch);

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
% color = [.85098 .8549 .85882];
% connections = [];
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');

%% Plot deviation from exact semisphere
% propSphere.radius = .075;
% propSphere.center = [0; 0; 0];
% [relErrorInL2Domain,relErrorInL2WeakDBC,errorInL2Interface,index] = ... 
%     plot_deviationFromSphereBSPlineSurface...
%     (BSplinePatches,'undefined',connections,propSphere,int,graph,'outputEnabled');

%% Write the case for Carat++
% for iPatches = 1:noPatches
%     isOnXi = false;
%     if BSplinePatches{iPatches}.weakDBC.etaExtension{1}(1,1) == ...
%             BSplinePatches{iPatches}.weakDBC.etaExtension{1}(1,2)
%         isOnXi = true;
%     end
%     if isOnXi
%         polOrder = BSplinePatches{iPatches}.p;
%     else
%         polOrder = BSplinePatches{iPatches}.q;
%     end
%     weakDBC.alpha = norm(eval('Dm'))*polOrder/...
%         BSplinePatches{iPatches}.minElArea;
% end
% connections = [];
% strongDBC = {homDOFs};
% pathToOutput = '../../../inputCarat/isogeometricMembraneAnalysis/';
% writeOutMultipatchBSplineSurface4Carat(BSplinePatches,strongDBC,connections,pathToOutput,caseName);

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-7; % 1e-10

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 10000;

%% Solve the steady-state problem
plot_IGANLinear = '';
[dHat,CPHistory,resHistory,hasConverged,BSplinePatches,minElASize] = ...
    solve_IGAMembraneNLinear(BSplinePatch,propNLinearAnalysis,...
    solve_LinearSystem,plot_IGANLinear,graph,'outputEnabled');
return;

%% Compute the displacement field at the middle of the patch

% Plot the current configuration and the resultants
scaling = 1e0;
graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
    (BSplinePatches,scaling*dHat,graph,'outputEnabled');

% Make a DOF numbering for the given patch
mxi = length(Xi);
meta = length(Eta);
nxi = length(BSplinePatch.CP(:,1,1));
neta = length(BSplinePatch.CP(1,:,1));
dHatElem = zeros(mxi - p - 1,meta - q - 1,3*(p + 1)*(q + 1));
for etaSpan = (q + 1):(meta - q - 1)
    for xiSpan = (p + 1):(mxi - p - 1)
        xiCounter = 1; 
        for c = etaSpan-q-1:etaSpan-1 
            for b = xiSpan-p:xiSpan
                dHatElem(xiSpan,etaSpan,xiCounter) = dHat(3*(c*nxi + b)-2);
                dHatElem(xiSpan,etaSpan,xiCounter + 1) = dHat(3*(c*nxi + b)-1);
                dHatElem(xiSpan,etaSpan,xiCounter + 2) = dHat(3*(c*nxi + b));
                
                % Update counter
                xiCounter = xiCounter + 3;
            end
        end
    end
end

% Get the data from the B-Spline patch
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Get the displacement vector at the top of the semisphere
xi = (Xi(1) + Xi(end))/2;
eta = (Eta(1) + Eta(end))/2;
xiSpan = findKnotSpan(xi,Xi,nxi);
etaSpan = findKnotSpan(eta,Eta,neta);
dHatActual = squeeze(dHatElem(xiSpan,etaSpan,:));
R = computeIGABasisFunctionsAndDerivativesForSurface...
	(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
dMiddle = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,R,dHatActual)

%% Compute the displacement field at the middle of the patch using Carat

% Read the displacement field from Carat
% displacementCarat = importdata('../../../outputIBRACarat/semisphereStatic_IGASinglePatch_1e-9_noElmnts625p3q3/displacements');
% displacementCarat = importdata('../../../outputIBRACarat/semisphereStatic_IGASinglePatch_1e-9_noElmnts625p2q2/displacements');
displacementCarat = importdata('../../../outputIBRACarat/semisphereStatic_IGASinglePatch_1e-9_noElmnts625p3q3/displacements');
displacementCarat = displacementCarat(:,2:4);
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
dHatCarat = zeros(3*nxi*neta,1);
counter = 1;
for iEta = 1:neta
    for iXi = 1:nxi
        dHatCarat(3*counter - 2,1) = displacementCarat(counter,1);
        dHatCarat(3*counter - 1,1) = displacementCarat(counter,2);
        dHatCarat(3*counter,1) = displacementCarat(counter,3);
        counter = counter + 1;
    end
end

% Plot the current configuration and the resultants
scaling = 1e1;
graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
    (BSplinePatches,scaling*dHatCarat,graph,'outputEnabled');

% Make a DOF numbering for the given patch
mxi = length(Xi);
meta = length(Eta);
nxi = length(BSplinePatch.CP(:,1,1));
neta = length(BSplinePatch.CP(1,:,1));
dHatElem = zeros(mxi - p - 1,meta - q - 1,3*(p + 1)*(q + 1));
for etaSpan = (q + 1):(meta - q - 1)
    for xiSpan = (p + 1):(mxi - p - 1)
        xiCounter = 1; 
        for c = etaSpan-q-1:etaSpan-1 
            for b = xiSpan-p:xiSpan
                dHatElem(xiSpan,etaSpan,xiCounter) = dHatCarat(3*(c*nxi + b)-2);
                dHatElem(xiSpan,etaSpan,xiCounter + 1) = dHatCarat(3*(c*nxi + b)-1);
                dHatElem(xiSpan,etaSpan,xiCounter + 2) = dHatCarat(3*(c*nxi + b));
                
                % Update counter
                xiCounter = xiCounter + 3;
            end
        end
    end
end

% Get the data from the B-Spline patch
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Get the displacement vector at the top of the semisphere
xi = (Xi(1) + Xi(end))/2;
eta = (Eta(1) + Eta(end))/2;
xiSpan = findKnotSpan(xi,Xi,nxi);
etaSpan = findKnotSpan(eta,Eta,neta);
dHatActual = squeeze(dHatElem(xiSpan,etaSpan,:));
R = computeIGABasisFunctionsAndDerivativesForSurface...
	(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
dMiddle = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,R,dHatActual);
dMiddle(3,1)

%% End
