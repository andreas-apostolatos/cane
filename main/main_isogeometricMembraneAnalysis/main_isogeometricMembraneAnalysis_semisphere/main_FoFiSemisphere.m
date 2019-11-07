%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : A single patch circular membrane form-found into a sphere and the 
%        results subsequently saved into a file to be read from an analysis 
%        file in order to perform a steady-state or a transient analysis.
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
        '../../../isogeometricThinStructureAnalysis/initialConditions/',....
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/');

%% NURBS parameters

% Global variables:
Radius = .075;
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
CP(:,:,4) = [1      weight        1
             weight weight^2      weight
             1      weight        1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i = 1:nxi
    for j = 1:neta
        if CP(i,j,4) ~= 1
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
parameters.E = 7.e5;

% Poisson ratio
parameters.nue = 0.49;

% Thickness of the shell
parameters.t = 1.65e-4;

% Density of the shell (used only for dynamics)
parameters.rho = 1050;

% Pressure Load amplitude
FAmp = - 19;

% Prestress for the membrane
% parameters.prestress.computeParametricCoordinates = @(X) [X(1,1)^2 + X(2,1)^2
%                                                           atan(X(2,1)/X(1,1))];
% parameters.prestress.computeBaseVectors = @(theta1,theta2) [ cos(theta2) -sin(theta2)
%                                                              sin(theta2) cos(theta2)
%                                                              0           0];
parameters.prestress.voigtVector = [abs(FAmp)*Radius/2/parameters.t
                                    abs(FAmp)*Radius/2/parameters.t
                                    0];

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'FoFiSemisphere';

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
a = 0; % 2
tp = a;
tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 3; % 25
refXi = scaling;
refEta = scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

%%%%%%%%%%%%%
% Dirichlet %
%%%%%%%%%%%%%

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
weakDBC.method = 'nitsche';
weakDBC.estimationStabilPrm = true;
weakDBC.xiExtension = {[0 0] [0 1] [1 1] [0 1]};
weakDBC.etaExtension = {[0 1] [0 0] [0 1] [1 1]};
weakDBC.int.type = 'default';
weakDBC.int.noGPs = 16;

%%%%%%%%%%
% Cables %
%%%%%%%%%%

cables.No = 0;

%%%%%%%%%%%
% Neumann %
%%%%%%%%%%%

xib = [0 1];   etab = [0 1];   dirForce = 'normal';
NBC.noCnd = 1;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC.isFollower(1,1) = true;
NBC.isTimeDependent(1,1) = false;

% xib = 0.499999999999954;   etab = 0.754462222569290;   dirForce = 'normal';
% NBC.noCnd = 1;
% NBC.xiLoadExtension = {xib};
% NBC.etaLoadExtension = {etab};
% NBC.loadAmplitude = {FAmp};
% NBC.loadDirection = {dirForce};
% NBC.computeLoadVct{1} = 'computeLoadVctPointIGAThinStructure';
% NBC.isFollower(1,1) = true;
% NBC.isTimeDependent(1,1) = false;

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);
BSplinePatches = {BSplinePatch};

%% Compute the load vector for the visualization of the reference configuration
BSplinePatches{1}.FGamma = [];
for counterNBC = 1:NBC.noCnd
    BSplinePatches{1}.FGamma = zeros(3*BSplinePatches{1}.noCPs,1);
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    BSplinePatches{1}.FGamma = ...
        funcHandle(BSplinePatches{1}.FGamma,BSplinePatches{1},...
        NBC.xiLoadExtension{counterNBC},...
        NBC.etaLoadExtension{counterNBC},...
        NBC.loadAmplitude{counterNBC},...
        NBC.loadDirection{counterNBC},...
        NBC.isFollower(counterNBC,1),...
        0,int,'outputEnabled');
end

%% Plot reference configuration before form-finding
color = [.85098 .8549 .85882];
% color = 'none';
connections = [];
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = -310; % -40
el = 40;
view(az,el);
% camlight(0,0);
% colorbar;
% lighting phong;
axis equal;
% xlim([-.075 .075]);
% ylim([-.075 .075]);
% zlim([-.001 .001]);

%% Plot the distribution of the Jacobian
figure(graph.index)
for iPatches = 1
    plot_postprocBSplineSurfaceGeometricJacobian(BSplinePatches{iPatches}.p,...
        BSplinePatches{iPatches}.q,BSplinePatches{iPatches}.Xi,...
        BSplinePatches{iPatches}.Eta,BSplinePatches{iPatches}.CP,...
        BSplinePatches{iPatches}.isNURBS);
    hold on;
end
hold off;
shading interp;
colormap('jet');
lighting phong;
az = -310; % -40
el = 40;
view(az,el);
camlight(0,0);
colorbar;
axis equal;
xlim([-.075 .075]);
ylim([-.075 .075]);
zlim([-.001 .001]);
graph.index = graph.index + 1;

%% Project point on the surface
node170 = sqrt(2)*Radius/2*[1.0; 0; 1.0];
xi0 = (BSplinePatches{1}.Xi(1) + BSplinePatches{1}.Xi(end))/2;
eta0 = (BSplinePatches{1}.Eta(1) + BSplinePatches{1}.Eta(end))/2;
propNewtonRaphson.eps = 1e-15;
propNewtonRaphson.maxIt = 1e2;
[xi,eta,Projected,isProjected,noIter] = computeNearestPointProjectionOnBSplineSurface...
    (node170,BSplinePatches{1}.p,BSplinePatches{1}.Xi,BSplinePatches{1}.q,...
    BSplinePatches{1}.Eta,BSplinePatches{1}.CP,BSplinePatches{1}.isNURBS,...
    xi0,eta0,propNewtonRaphson);

%% Perform form-finding analysis to find the shape of the membrane
propFormFinding.tolerance = 1e-11;
propFormFinding.maxNoIter = 1e4;
[BSplinePatch,CPHistory,resHistory,hasConverged,noIter] = ...
    solve_formFindingIGAMembrane(BSplinePatch,propFormFinding,...
    solve_LinearSystem,'outputEnabled');
BSplinePatches = {BSplinePatch};
save(['./data_FoFiPool/' 'data_' caseName '_Tol_' num2str(propFormFinding.tolerance) '_noElmnts_' num2str(BSplinePatches{1}.noElmnts) '_p_' num2str(BSplinePatches{1}.p) '_q_' num2str(BSplinePatches{1}.q)]);
return;

%% Plot deviation from semisphere
propSphere.center = [0; 0; 0];
propSphere.radius = Radius;
dHat = 'undefined';
connections = [];
[relErrorInL2Domain,relErrorInL2WeakDBC,errorInL2Interface,index] = ...
    plot_deviationFromSphereBSPlineSurface... 
    (BSplinePatches,dHat,connections,propSphere,int,graph,'outputEnabled');
az = -310; % -40
el = 40;
view(az,el);

%% Output the initial geometry as a Carat++ input file
pathToOutput = '../../../inputCarat/';
writeOutMultipatchBSplineSurface4Carat(BSplinePatches,pathToOutput,caseName);

%% Plot the reference configuration after form-finding
color = [.85098 .8549 .85882];
connections = [];
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');

%% Read results from a Carat++ result file

% Define the path to the result file and the file name
pathToResultFile = '../../../resultsCarat/isogeometricMembraneAnalysis/';
fileName = 'resultsSemisphere256ElmntsP3Q3';

% Read the results file
fstring = fileread([pathToResultFile fileName '.res']);

% Get the results in terms of the displacement field
block = regexp(fstring,'Rhino Post Results File 1.0\nResult "Displacement" "Load Case" 1 Vector OnNodes\nValues','split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
dispCarat = out(:,2:4);
noDOFsCarat = 3*length(dispCarat(:,1));
if noDOFsCarat ~= 3*BSplinePatch.noCPs
    error('B-Spline patch data and carat results are not in agreement');
end
dHatCarat = zeros(noDOFsCarat,1);
for iCPs = 1:BSplinePatch.noCPs
    dHatCarat(3*iCPs - 2) = dispCarat(iCPs,1);
    dHatCarat(3*iCPs - 1) = dispCarat(iCPs,2);
    dHatCarat(3*iCPs) = dispCarat(iCPs,3);
end

%% Postprocessing

% Plot the current configuration and the resultants
scaling = 1e0;
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
XMiddle = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R)
dMiddle = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,R,dHatActual)

%% End
