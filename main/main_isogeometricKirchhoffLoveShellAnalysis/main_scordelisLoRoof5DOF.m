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
%        solved in a single patch geometry using 5 DOFs per control point.
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
    homDOFs = findDofs5D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 1];    
for dirSupp = [2 3]
    homDOFs = findDofs5D(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Fix the back left corner of the shell to avoid rigid body motions
xiSup = [0 0];   etaSup = [0 0];   dirSupp = 1;
homDOFs = findDofs5D(homDOFs,xiSup,etaSup,dirSupp,CP);

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
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure5DOF';
NBC.isConservative(1,1) = true;
NBC.isTimeDependent(1,1) = false;

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);

% Set the number of DOFs for 5 DOF per control point
BSplinePatch.noDOFs = 5*BSplinePatch.noCPs;

%% Compute the load vectors for each patch (only for the visualization)
FGamma = zeros(5*BSplinePatch.noCPs,1);
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    FGamma = funcHandle...
        (FGamma,BSplinePatch,NBC.xiLoadExtension{counterNBC},...
        NBC.etaLoadExtension{counterNBC},NBC.loadAmplitude{counterNBC},...
        NBC.loadDirection{counterNBC},NBC.isFollower(counterNBC,1),...
        0,BSplinePatch.int,'outputEnabled');
end
    
%% Plot reference configuration
% Filter homDOFs to only include translational DOFs (1,2,3) for plotting
% since the plotting function expects 3 DOFs per control point
homDOFsFiltered = [];
for i = 1:length(homDOFs)
    dofIndex = homDOFs(i);
    % Convert to control point and direction
    p_cp = ceil(dofIndex/5);
    dir = dofIndex - 5*(p_cp-1);
    % Only keep translational DOFs (dir = 1,2,3)
    if dir <= 3
        % Convert back to 3 DOF numbering for plotting
        dof3D = 3*(p_cp-1) + dir;
        homDOFsFiltered = [homDOFsFiltered; dof3D];
    end
end

% Also filter FGamma to match 3 DOF structure for plotting
noCPs = BSplinePatch.noCPs;
FGammaFiltered = zeros(3*noCPs, 1);
for i = 1:noCPs
    % Extract translational components from 5 DOF vector
    FGammaFiltered(3*(i-1)+1) = FGamma(5*(i-1)+1); % x
    FGammaFiltered(3*(i-1)+2) = FGamma(5*(i-1)+2); % y
    FGammaFiltered(3*(i-1)+3) = FGamma(5*(i-1)+3); % z
end

figure(graph.index)
plot_referenceConfigurationIGAThinStructure(p,q,Xi,Eta,CP,isNURBS,homDOFsFiltered,FGammaFiltered,'outputEnabled');
title('Reference configuration for an isogeometric Kirchhoff-Love shell (5 DOF)');
graph.index = graph.index + 1;

%% Write geometry for Carat
% BSplinePatches = {BSplinePatch};
% strongDBC = [];
% connections = [];
% pathToOutput = '../../outputIBRACarat/';
% caseName = 'scordelisLoRoof';
% writeOutMultipatchBSplineSurface4Carat...
%     (BSplinePatches,strongDBC,connections,pathToOutput,caseName);

%% Solve the system applying linear analysis
[dHatLinear,F,minElArea] = solve_IGAKirchhoffLoveShellLinear5DOF...
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
% Note: For nonlinear analysis, you would need to create the 5DOF version of the nonlinear solver as well
 [dHatNLinear, CPHistory, resHistory, isConverged, BSplinePatch, minElASize] = ...
     solve_IGAKirchhoffLoveShellNLinear5DOF ...
     (BSplinePatch, propNLinearAnalysis, solve_LinearSystem, ...
     plot_IGANonlinear, graph, 'outputEnabled');

%% Postprocessing
% graph.index = plot_postprocIGAKirchhoffLoveShellLinear...
%     (BSplinePatch,dHatLinear,graph,'outputEnabled');
% title('Linear analysis (5 DOF)');
% graph.index = plot_postprocIGAKirchhoffLoveShellNLinear(BSplinePatch,dHatNLinear,graph,'outputEnabled');
% title('Nonlinear analysis (5 DOF)');

fprintf('Linear analysis completed with 5 DOFs per control point.\n');
fprintf('Displacement field size: %d x %d\n', size(dHatLinear));
fprintf('Maximum displacement magnitude: %.6e\n', max(abs(dHatLinear)));

% Display the difference in system size compared to 3 DOF version
noCPs = BSplinePatch.noCPs;
fprintf('\nSystem size comparison:\n');
fprintf('3 DOF version would have: %d DOFs\n', 3*noCPs);
fprintf('5 DOF version has: %d DOFs\n', 5*noCPs);
fprintf('Increase in system size: %.1f%%\n', (5*noCPs - 3*noCPs)/(3*noCPs)*100);

%% Save data
% save overkillSolutionScordelisLoRoof5DOFKLShellp10q10nXi59nEta39;

%% end
