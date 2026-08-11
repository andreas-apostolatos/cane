%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Script documentation
% 
% Task : The Scordelis-Lo-Roof benchmark problem is modelled and
%        solved using HSDT (Higher-order Shear Deformation Theory)
%        with 5 DOF per control point: [u, v, w, θx, θy]
%
% Date : [Current Date]
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
    
% Add all functions related to the Isogeometric HSDT shell formulation
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

% Shear correction factor for HSDT (typical values: 5/6, π²/12, 1.0)
parameters.shearCorrection = 5/6;  % Reissner-Mindlin default

%% GUI

% Analysis type
analysis.type = 'isogeometricHSDTShellAnalysis';

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
tp = 1; 
tq = 1; 
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 1; 
edgeRatio = ceil(Length/Radius/(sin(4*pi/9)));
refXi = edgeRatio*scaling;
refEta = ceil(4/3)*scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions for HSDT (5 DOF per CP)

% supports (Dirichlet boundary conditions)
% Note: HSDT has 5 DOF per control point [u, v, w, θx, θy]

% back and front curved edges are a rigid diaphragm
homDOFs = [];
xiSup = [0 0];   etaSup = [0 1];    
for dirSupp = [2 3]  % constrain v and w displacements
    homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 1];    
for dirSupp = [2 3]  % constrain v and w displacements
    homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Fix the back left corner of the shell to avoid rigid body motions
xiSup = [0 0];   etaSup = [0 0];   dirSupp = 1;  % constrain u displacement
homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);

% Additional constraints for rotations if needed (typically at supports)
% Example: Fix rotations at corners to prevent rigid body rotation
% xiSup = [0 0];   etaSup = [0 0];   
% for dirSupp = [4 5]  % constrain θx and θy rotations
%     homDOFs = findDofs5D_HSDT(homDOFs,xiSup,etaSup,dirSupp,CP);
% end

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
xib = [0 1];   etab = [0 1];   dirForce = 'z';  % vertical load
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.isFollower(1,1) = false;
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAHSDTShell';  % Updated for HSDT
NBC.isConservative(1,1) = true;
NBC.isTimeDependent(1,1) = false;

%% Create the B-Spline patch array for HSDT
BSplinePatch = fillUpPatch_HSDT...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);

%% Compute the load vectors for each patch (only for the visualization)
FGamma = zeros(5*BSplinePatch.noCPs,1);  % 5 DOF per control point
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
title('Reference configuration for an isogeometric HSDT shell (5 DOF per CP)');
graph.index = graph.index + 1;

%% Solve the system applying linear HSDT analysis
fprintf('Starting HSDT linear analysis with 5 DOF per control point...\n');
fprintf('DOF ordering: [u, v, w, θx, θy] per control point\n\n');

[dHatLinear_HSDT,F_HSDT,minElArea_HSDT] = solve_IGAHSDTShellLinear...
    (BSplinePatch,solve_LinearSystem,'outputEnabled');

%% Display results summary
fprintf('=== HSDT Analysis Results Summary ===\n');
fprintf('Total DOFs: %d (5 per control point)\n', length(dHatLinear_HSDT));
fprintf('Number of control points: %d x %d = %d\n', nxi, neta, nxi*neta);
fprintf('Minimum element area: %.6e\n', minElArea_HSDT);

% Extract displacement components
numCPs = nxi * neta;
u_HSDT = dHatLinear_HSDT(1:5:end);     % u-displacements
v_HSDT = dHatLinear_HSDT(2:5:end);     % v-displacements  
w_HSDT = dHatLinear_HSDT(3:5:end);     % w-displacements
theta_x_HSDT = dHatLinear_HSDT(4:5:end); % θx-rotations
theta_y_HSDT = dHatLinear_HSDT(5:5:end); % θy-rotations

fprintf('\nDisplacement ranges:\n');
fprintf('u: [%.6e, %.6e]\n', min(u_HSDT), max(u_HSDT));
fprintf('v: [%.6e, %.6e]\n', min(v_HSDT), max(v_HSDT));
fprintf('w: [%.6e, %.6e]\n', min(w_HSDT), max(w_HSDT));
fprintf('θx: [%.6e, %.6e] rad\n', min(theta_x_HSDT), max(theta_x_HSDT));
fprintf('θy: [%.6e, %.6e] rad\n', min(theta_y_HSDT), max(theta_y_HSDT));

% Maximum displacement magnitude
max_w = max(abs(w_HSDT));
fprintf('\nMaximum |w| displacement: %.6e\n', max_w);
fprintf('Maximum |θx| rotation: %.6e rad\n', max(abs(theta_x_HSDT)));
fprintf('Maximum |θy| rotation: %.6e rad\n', max(abs(theta_y_HSDT)));

%% Postprocessing (simplified - full postprocessing functions would need HSDT adaptation)
% graph.index = plot_postprocIGAHSDTShellLinear...
%     (BSplinePatch,dHatLinear_HSDT,graph,'outputEnabled');
% title('HSDT Linear Analysis Results');

%% Save data for comparison
save('scordelisLoRoof_HSDT_Results.mat', 'dHatLinear_HSDT', 'F_HSDT', ...
     'minElArea_HSDT', 'u_HSDT', 'v_HSDT', 'w_HSDT', 'theta_x_HSDT', ...
     'theta_y_HSDT', 'parameters', 'CP', 'Xi', 'Eta', 'p', 'q');

fprintf('\n=== HSDT Analysis Completed Successfully ===\n');
fprintf('Results saved to: scordelisLoRoof_HSDT_Results.mat\n\n');

%% end