%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
%% Script documentation
% 
% Task : In this script it is modelled a cantilever plate which is force to
%        bend into a circle by allying a constant bending moment at its
%        tip. The applied analysis is non linear.
%
% Date : 10.12.2013
%
%% Preamble
clear;
clc;

%% Includes 

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Add the linear equation solvers
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
        '../../isogeometricThinStructureAnalysis/nonConservativeLoads/');

%% NURBS parameters

% Global variables:
Length = 1;
% Length = 12;
Width = 1e-1;
% Width = 1;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [0  0
             1  1]*Length;
         
% y-coordinates
CP(:,:,2) = [-Width/2 Width/2
             -Width/2 Width/2];
         
% z-coordinates
CP(:,:,3)=[0 0
           0 0];
       
% Weights
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
parameters.E = 2e8;
% parameters.E = 1.2e6;
% parameters.E = 1e2;

% Poisson ratio
parameters.nue = 0.0;

% Thickness of the shell
parameters.t = 1e-2;

% Density of the shell (used only for dynamics)
parameters.rho = 7810;

%% GUI

% Analysis type
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

% Define linear system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'default';
if strcmp(int.type,'manual')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.etaNGPForLoad = 6;
    int.nGPForLoad = 10;
end

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','rotation','curvature','force',
% 'moment','shearForce'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x', 'y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = 'z';

%% Refinement

% Degree by which to elevate
tp = 3;
tq = 3;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
refXi = 1;
refEta = 1;

% scaling = 1;
% refXi = 12*scaling;
% refEta = 0*scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

% Insert an additional knot in Xi knot vector
% newKnot = Xi(length(Xi)) - 1e-2;
% [Xi,Eta,CP] = knotRefineBSplineSurface(p,q,Xi,Eta,CP,[newKnot],[],'');

%% Dirichlet and Neumann boundary conditions 

% supports (Dirichlet boundary conditions)

% Get the number of Control points in the xi-direction
nxi = length(CP(:,1,1));

% Get the number of Control points in the eta-direction
neta = length(CP(1,:,1));

% Homogeneous Dirichlet boundary conditions
homDOFs = [];
xiSup = [0 0]; etaSup = [0 1];
homDOFs = findDOFsClampedSupportIGAKLShell...
    (homDOFs,xiSup,etaSup,p,Xi,q,Eta,CP,isNURBS);
xiSup = [0 1]; etaSup = [0 1];  dir = 2;
homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dir,CP);

% xiSup = [0 1];   etaSup = [0 1];    dirSupp = 2;
% homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Cables
cables.No = 0;

% load (Neuman boundary conditions)
% Radius = Length/2/pi();
% momentAmp = parameters.E*parameters.t^3/12/Radius;
% GShear = parameters.E/2/(1 + parameters.nue);
% theta = 2*pi();
% bendMomentAmp = theta*GShear*It/Length;
% bendMomentAmp = .01;
% twistMomentAmp = .01;
% FAmp = .6e-2;
% pressureAmp = - 1e-3;

% Apply a horizontal tip force towards +x
% xib = 1;   etab = [0 1]; dirForce = 3;
% Fl = computeLoadVctLineIGAKirchhoffLoveShell(Fl,xib,etab,p,q,Xi,Eta,CP,isNURBS,FAmp,dirForce,int,'outputEnabled');

% Apply a tip bending moment
M0 = 25/3;
% bendMomentAmp = - 2*pi()*M0;
loadAmp = 1e-1;
% xib = 1;   etab = [0 1];   dirBendMoment = 1;
xib = 1;   etab = [0 1];   dirForce = 'z';
NBC.noCnd = 1;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
% NBC.loadAmplitude = {bendMomentAmp};
NBC.loadAmplitude = {loadAmp};
% NBC.loadDirection(1,1) = dirBendMoment;
NBC.loadDirection(1,1) = dirForce;
% NBC.computeLoadVct{1} = 'computeLoadVctLineMomentIGAKirchhoffLoveShell';
NBC.computeLoadVct{1} = 'computeLoadVctLineIGAThinStructure';
NBC.isFollower(1,1) = false;

% Apply a tip twisting moment
% xib = 1;   etab = [0 1];   dirBendMoment = 1; dirTwistMoment = 2;
% Fl = computeLoadVctLineMomentIGAKirchhoffLoveShell(Fl,xib,etab,p,q,Xi,Eta,CP,isNURBS,twistMomentAmp,dirTwistMoment,int,'outputEnabled');

% Apply a pressure along the shell's surface towards -x
% xib = [0 1]; etab = [0 1]; dirForce = 3;
% Fl = computeLoadVctAreaIGAKirchhoffLoveSell(Fl,xib,etab,p,q,Xi,Eta,CP,isNURBS,pressureAmp,dirForce,int,'outputEnabled');

%% Create the B-Spline patch array
BSplinePatch = ...
    fillUpPatch ...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);

%% Compute the load vector for the visualization of the reference configuration
FGamma = [];
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    FGamma = funcHandle(FGamma,BSplinePatch,NBC.xiLoadExtension{counterNBC},...
                    NBC.etaLoadExtension{counterNBC},...
                    NBC.loadAmplitude{counterNBC},...
                    NBC.loadDirection(counterNBC,1),...
                    NBC.isFollower(counterNBC),0,int,'outputEnabled');
end

%% Plot reference configuration
figure(graph.index)
plot_referenceConfigurationIGAThinStructure...
    (p,q,Xi,Eta,CP,isNURBS,homDOFs,FGamma,'outputEnabled');
title('Reference configuration for an isogeometric Kirchhoff-Love shell');
graph.index = graph.index + 1;

%% Solve the system applying linear analysis
% [dHat,F,minElArea] = solve_IGAKirchhoffLoveShellLinear(BSplinePatch,'outputEnabled');

%% Solve the system applying nonlinear analysis

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 4;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 10e-9;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 120;

%% Solve the linear system
[dHatLinear,FComplete,minElSize] = solve_IGAKirchhoffLoveShellLinear...
    (BSplinePatch,solve_LinearSystem,'outputEnabled');
minElSize

%% Solve the nonlinear system
% plot_IGANonlinear = '';
% % plot_IGANonlinear = 'plot_postprocIGAKirchhoffLoveShellNLinear';
% [dHat,FComplete,CpHistory,resHistory] = solve_IGAKirchhoffLoveShellNLinear...
%     (BSplinePatch,propNLinearAnalysis,solve_LinearSystem,plot_IGANonlinear,...
%     graph,'outputEnabled');

%% Postprocessing
BSplinePatch.p = p;
BSplinePatch.q = q;
BSplinePatch.Xi = Xi;
BSplinePatch.Eta = Eta;
BSplinePatch.CP = CP;
BSplinePatch.isNURBS = isNURBS;
BSplinePatch.homDOFs = homDOFs;
BSplinePatch.parameters = parameters;
BSplinePatch.FGamma = FGamma;

% Linear analysis
graph.index = plot_postprocIGAKirchhoffLoveShellLinear(BSplinePatch,dHatLinear,graph,'outputEnabled');

% Compute the displacement at the middle of the back curved edge:

% Preliminary parameters
mxi = length(Xi);
meta = length(Eta);
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Make a DOF numbering for the given patch
dHatElem = zeros(mxi-p-1,meta-q-1,3*(p+1)*(q+1));
for etaSpan = (q+1):(meta-q-1)
    for xiSpan = (p+1):(mxi-p-1)
        xiCounter = 1; 
        for c = etaSpan-q-1:etaSpan-1 
            for b = xiSpan-p:xiSpan
                dHatElem(xiSpan,etaSpan,xiCounter) = dHatLinear(3*(c*nxi+b)-2);
                dHatElem(xiSpan,etaSpan,xiCounter + 1) = dHatLinear(3*(c*nxi+b)-1);
                dHatElem(xiSpan,etaSpan,xiCounter + 2) = dHatLinear(3*(c*nxi+b));
                
                % Update counter
                xiCounter = xiCounter + 3;
            end
        end
    end
end

% Parametetric location of the middle point of the curves edge
xi = Xi(end);
eta = (Eta(1) + Eta(end))/2;
xiSpan = findKnotSpan(xi,Xi,nxi);
etaSpan = findKnotSpan(eta,Eta,neta);
dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
dHatActual(:,1) = dHatElem(xiSpan,etaSpan,:);
d = computePostprocDisplacementIGAKirchhoffLoveShell(p,q,dR(:,1),dHatActual);
d(3,1)

dExact = 0.002;

relError = abs(d(3,1) - dExact)/abs(dExact)

ref_minElSize = [0.050000000000000 0.006250000000000 0.002000000018774];
ref_relError = [5.576051773092949e-12 3.045198641860036e-11 9.386831372956705e-09];

loglog(ref_minElSize,ref_relError)

% Nonlinear analysis
% graph.index = plot_postprocIGAKirchhoffLoveShellNLinear(BSplinePatch,dHat,graph,'outputEnabled');

%% end
