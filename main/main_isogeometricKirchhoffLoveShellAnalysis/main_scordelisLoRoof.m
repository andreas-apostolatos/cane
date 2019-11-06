%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : The full Scordelis-Lo-Roof benchmark problem is modelled and
%        solved in a single patch geometry.
%
% Date : 12.02.2015
%
%% Preamble
clear;
clc;

%% Includes 

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Include linear equation system solvers
addpath('../../equationSystemSolvers/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the Isogeometric Kirchhoff-Love shell formulation
addpath('../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../isogeometricThinStructureAnalysis/loads/',...
        '../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../isogeometricThinStructureAnalysis/solvers/',...
        '../../isogeometricThinStructureAnalysis/metrics/',...
        '../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../isogeometricThinStructureAnalysis/output/');

%% NURBS parameters

% Global variables
Length = 50;
Radius = 25;

% Polynomial degrees
p = 1;
q = 2;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [-Length/2 -Length/2 -Length/2
             Length/2  Length/2  Length/2];
         
% y-coordinates
CP(:,:,2) = [-Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)
             -Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)];
         
% z-coordinates
CP(:,:,3) = [Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)
             Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)];
       
% Weights
weight = cos(2*pi/9);
CP(:,:,4) = [1 weight 1
             1 weight 1];

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
parameters.E = 4.32e8;

% Poisson ratio
parameters.nue = .0;

% Thickness of the shell
parameters.t = .25;

% Density of the shell (used only for dynamics)
parameters.rho = 7850;

%% GUI

% Analysis type
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points
int.type = 'default';
if strcmp(int.type,'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.xetaNGPForLoad = 6;
end

% Equation system solvers
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

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
graph.component = 'z';

%% Refinement

% Degree by which to elevate
tp = 1; % 9
tq = 1; % 8 
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 1; % 14
edgeRatio = ceil(Length/Radius/(sin(4*pi/9)));
refXi = edgeRatio*scaling;
refEta = ceil(4/3)*scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

% supports (Dirichlet boundary conditions)

% back and front curved edges are a rigid diaphragm
homDOFs = [];
xiSup = [0 0];   etaSup = [0 1];    
for dirSupp = [2 3]
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 1];    
for dirSupp = [2 3]
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Fix the back left corner of the shell to avoid rigid body motions
xiSup = [0 0];   etaSup = [0 0];   dirSupp = 1;
homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables.No = 0;

% load (Neuman boundary conditions)
FAmp = - 9e1;
NBC.noCnd = 1;
xib = [0 1];   etab = [0 1];   dirForce = 'z';
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.isFollower(1,1) = false;
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC.isConservative(1,1) = true;
NBC.isTimeDependent(1,1) = false;

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);

%% Compute the load vectors for each patch (only for the visualization)
FGamma = zeros(3*BSplinePatch.noCPs,1);
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    FGamma = funcHandle...
        (FGamma,BSplinePatch,NBC.xiLoadExtension{counterNBC},...
        NBC.etaLoadExtension{counterNBC},NBC.loadAmplitude{counterNBC},...
        NBC.loadDirection{counterNBC},NBC.isConservative(counterNBC,1),...
        0,BSplinePatch.int,'outputEnabled');
end
    
%% Plot reference configuration
figure(graph.index)
plot_referenceConfigurationIGAThinStructure(p,q,Xi,Eta,CP,isNURBS,homDOFs,FGamma,'outputEnabled');
title('Reference configuration for an isogeometric Kirchhoff-Love shell');
graph.index = graph.index + 1;

%% Write geometry for Carat
BSplinePatches = {BSplinePatch};
strongDBC = [];
connections = [];
pathToOutput = '../../outputIBRACarat/';
caseName = 'scordelisLoRoof';
writeOutMultipatchBSplineSurface4Carat...
    (BSplinePatches,strongDBC,connections,pathToOutput,caseName);

%% Solve the system applying linear analysis
[dHatLinear,F,minElArea] = solve_IGAKirchhoffLoveShellLinear...
    (BSplinePatch,solve_LinearSystem,'outputEnabled');

%% Solve the system applying nonlinear analysis

% Nonlinear analysis method
propNLinearAnalysis.method = 'newtonRaphson';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 20;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 10e-7;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 120;

% On the load conservativeness for each patch
noPatches = 1;
propNLinearAnalysis.conservativeLoad = true(noPatches,1);

%% Solve the nonlinear system
plot_IGANonlinear = 'undefined';
[dHatNLinear,CPHistory,resHistory,hasConverged,minElASize] = solve_IGAKirchhoffLoveShellNLinear...
    (BSplinePatch,propNLinearAnalysis,solve_LinearSystem,plot_IGANonlinear,graph,'outputEnabled');

%% Postprocessing
graph.index = plot_postprocIGAKirchhoffLoveShellLinear...
    (BSplinePatch,dHatLinear,graph,'outputEnabled');
title('Linear analysis');
graph.index = plot_postprocIGAKirchhoffLoveShellNLinear(BSplinePatch,dHatNLinear,graph,'outputEnabled');
title('Nonlinear analysis');

% % Compute the displacement at the middle of the back curved edge:
% 
% % Preliminary parameters
% mxi = length(Xi);
% meta = length(Eta);
% nxi = length(CP(:,1,1));
% neta = length(CP(1,:,1));
% 
% % Make a DOF numbering for the given patch
% dHatElem = zeros(mxi-p-1,meta-q-1,3*(p+1)*(q+1));
% for etaSpan = (q+1):(meta-q-1)
%     for xiSpan = (p+1):(mxi-p-1)
%         xiCounter = 1; 
%         for c = etaSpan-q-1:etaSpan-1 
%             for b = xiSpan-p:xiSpan
%                 dHatElem(xiSpan,etaSpan,xiCounter) = dHatLinear(3*(c*nxi+b)-2);
%                 dHatElem(xiSpan,etaSpan,xiCounter + 1) = dHatLinear(3*(c*nxi+b)-1);
%                 dHatElem(xiSpan,etaSpan,xiCounter + 2) = dHatLinear(3*(c*nxi+b));
%                 
%                 % Update counter
%                 xiCounter = xiCounter + 3;
%             end
%         end
%     end
% end
% 
% % Parametetric location of the middle point of the curves edge
% xi = (Xi(1) + Xi(end))/2;
% % xi = Xi(1);
% % eta = (Eta(1) + Eta(end))/2;
% eta = Eta(1);
% xiSpan = findKnotSpan(xi,Xi,nxi);
% etaSpan = findKnotSpan(eta,Eta,neta);
% dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
% dHatActual(:,1) = dHatElem(xiSpan,etaSpan,:);
% d = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,dR(:,1),dHatActual);
% d(3,1)

%% Save data
% save overkillSolutionScordelisLoRoofKLShellp10q10nXi59nEta39;

%% end