%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Form-finding analysis over the four-point sail modelled with a
%        single patch
%
% Date : 11.01.2016
%
%% Preamble
clear;
clc;

%% Includes

% Add general math functions
addpath('../../../generalMath/');

% Add general auxiliary functions
addpath('../../../auxiliary/');

% Add all linear equation system solvers
addpath('../../../equationSystemSolvers/');

% Add all efficient computation functions
addpath('../../../efficientComputation/');

% Add graphics-related functions for the isogeometric mortar-based mapping
addpath('../../../isogeometricMortarBasedMappingAnalysis/graphics/');

% Add functions related to writting out
addpath('../../../outputIBRACarat/');

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
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/output/');

%% NURBS parameters

% Global variables:
Length = 20;
Width = Length;
Height = Length/2;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [-Length/2 -Length/2
             Length/2  Length/2];
         
% y-coordinates
CP(:,:,2) = [-Width/2 Width/2
             -Width/2 Width/2];
         
% z-coordinates
CP(:,:,3) = [0      Height
             Height 0];
       
% Weights
weight = sqrt(2)/2;
CP(:,:,4) = [1 1
             1 1];

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
parameters.E = 8e+8;

% Poisson ratio
parameters.nue = .4;

% Thickness
parameters.t = 1e-3;

% Density (used only for dynamics)
parameters.rho = 8050;

% Prestress for the membrane
sigma0 = 3e+3/parameters.t; % 
parameters.prestress.voigtVector = [sigma0
                                    sigma0
                                    0];

% Cable parameters
parametersCable.E = 1.6e+11;
parametersCable.radiusCS = 12e-3/2;
parametersCable.areaCS = pi*parametersCable.radiusCS^2;
parametersCable.rho = 8050;
parametersCable.prestress = 6e+4/parametersCable.areaCS;
cables.parameters = {parametersCable parametersCable parametersCable parametersCable};

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'FoFiFourPointSail';

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

%% Refinement

% Degree by which to elevate
a = 0; % 2
tp = a;
tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 10; % 50
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
% % weakDBC.method = 'nitsche';
% weakDBC.method = 'penalty';
% weakDBC.estimationStabilPrm = false;
% weakDBC.alpha = 1e11;
% weakDBC.xiExtension = {[0 0] [0 0] [1 1] [1 1]};
% weakDBC.etaExtension = {[0 0] [1 1] [0 0] [1 1]};
% weakDBC.int.type = 'default';
% weakDBC.int.noGPs = 16;

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
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,valuesInhomDOFs,...
    weakDBC,cables,NBC,[],[],[],[],[],int);
BSplinePatches = {BSplinePatch};
noPatches = length(BSplinePatches);
connections = [];

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

%% Compute the load vector for the visualization of the reference configuration
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.FGamma = ...
        zeros(3*BSplinePatches{iPatches}.noCPs,1);
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

%% Plot the reference configuration for the multipatch geometry before the form-finding analysis
color = [217 218 219]/255;
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = -40;
el = 40;
view(az,el);
axis off;
title('');
camlight left;
lighting phong;

%% Nonlinear analysis properties

% % number of load steps for the non-linear analysis
% propNLinearAnalysis.method = 'newtonRapshon';
% 
% % number of load steps for the non-linear analysis
% propNLinearAnalysis.noLoadSteps = 1;
% 
% % Assign a tolerance for the Newton iterations
% propNLinearAnalysis.eps = 1.0e-9;
% 
% % Assign the maximum number of iterations
% propNLinearAnalysis.maxIter = 1000;

%% Solve the system applying nonlinear analysis
% % plot_IGANLinear = @plot_postprocIGAMembraneNLinear;
% plot_IGANLinear = '';
% [dHatNLinear,CPHistory,resHistory,hasConverged,minElArea] = solve_IGAMembraneNLinear...
%     (BSplinePatch,propNLinearAnalysis,solve_LinearSystem,plot_IGANLinear,...
%     graph,'outputEnabled');
% scaling = 1;
% graph.index = plot_postprocIGAMembraneNLinear...
%     (BSplinePatch,scaling*dHatNLinear,graph,'outputEnabled');
% return;

%% Output the initial geometry as a Carat++ input file
% pathToOutput = '../../../inputCarat/';
% writeOutMultipatchBSplineSurface4Carat(BSplinePatches,pathToOutput,caseName);

%% Perform a form finding analysis to find the shape of the membrane
propFormFinding.tolerance = 1e-4;
propFormFinding.maxNoIter = 100; % 100
[BSplinePatch,CPHistory,resHistory,hasConverged,noIter] = ...
    solve_formFindingIGAMembrane(BSplinePatch,propFormFinding,...
    solve_LinearSystem,'outputEnabled');
BSplinePatches = {BSplinePatch};
% save data_FoFiFourPointSailCoarse.mat;
% return;

%% Compute the load vector for the visualization of the reference configuration
FGamma = [];
for iNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{iNBC});
    FGamma = funcHandle(FGamma,BSplinePatch,NBC.xiLoadExtension{iNBC},...
                    NBC.etaLoadExtension{iNBC},...
                    NBC.loadAmplitude{iNBC},...
                    NBC.loadDirection{iNBC},...
                    NBC.isFollower(iNBC,1),0,int,'outputEnabled');
end

%% Compute the load vector for the visualization of the reference configuration
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.FGamma = ...
        zeros(3*BSplinePatches{iPatches}.noCPs,1);
%     for iNBC = 1:NBC{iPatches}.noCnd
%         funcHandle = str2func(NBC{iPatches}.computeLoadVct{iNBC});
%         BSplinePatches{iPatches}.FGamma = funcHandle...
%             (BSplinePatches{iPatches}.FGamma,...
%             BSplinePatches{iPatches},...
%             NBC{iPatches}.xiLoadExtension{iNBC},...
%             NBC{iPatches}.etaLoadExtension{iNBC},...
%             NBC{iPatches}.loadAmplitude{iNBC},...
%             NBC{iPatches}.loadDirection{iNBC},...
%             NBC{iPatches}.isFollower(iNBC,1),0,...
%             BSplinePatches{iPatches}.int,'outputEnabled');
%     end
end

%% Plot the reference configuration for the multipatch geometry after the form-finding analysis
color = [217 218 219]/255;
% color = 'none';
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = 260; % 160
el = 30; % 30
view(az,el);
camlight(30,70); % (30,60)
lighting phong;
axis off;
title('');

%% End
