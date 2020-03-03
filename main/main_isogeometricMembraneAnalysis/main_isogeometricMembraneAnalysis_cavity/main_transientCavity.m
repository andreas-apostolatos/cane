%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Solve the cavity FSI benchmark using a 1 patch membrane.
%
% Date : 09.10.2017
%
%% Preamble
clear;
clc;
close;

%% Includes

% Add general math functions
addpath('../../../generalMath/');

% Add general auxiliary functions
addpath('../../../auxiliary/');

% Add system solvers
addpath('../../../equationSystemSolvers/');

% Add efficient computation functions
addpath('../../../efficientComputation/');

% Add transient analysis solvers
addpath('../../../transientAnalysis/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../../CAGDKernel/CAGDKernel_graphics/',...
        '../../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the isogeometric Kirchhoff-Love shell formulation
addpath('../../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../../isogeometricThinStructureAnalysis/graphicsMultipatches/',...
        '../../../isogeometricThinStructureAnalysis/loads/',...
        '../../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../../isogeometricThinStructureAnalysis/solvers/',...
        '../../../isogeometricThinStructureAnalysis/metrics/',...
        '../../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../../isogeometricThinStructureAnalysis/nonConservativeLoads/',...
        '../../../isogeometricThinStructureAnalysis/penaltyDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/penaltyDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../../isogeometricThinStructureAnalysis/output/',...
        '../../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../../isogeometricThinStructureAnalysis/penaltyWeakDirichletBoundaryConditionsKLShell/',...
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/');
    
% Add functions related to the visualization of the multipatch geometry
% related to the isogeometric mortar-based mapping
addpath('../../../isogeometricMortarBasedMappingAnalysis/graphics/');

% Add functions related to the co-simulation with Empire
addpath('../../../empireCosimulationAnalysis/');

%% CAD modelling of the geometry via NURBS

% Global variables:
Length = 1;
Width = 0.1;

% Patch 1 :
% _________

% Polynomial degrees
p1 = 1;
q1 = 1;

% Knot vectors
Xi1 = [0 0 1 1];
Eta1 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP1(:,:,1) = [0 0
              1 1]*Length;

% y-coordinates
CP1(:,:,2) = [0 0
              0 0];


% z-coordinates
CP1(:,:,3) = [0 1
              0 1]*Width;

% Weights
CP1(:,:,4) = [1 1
              1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = 0;
nxi1 = length(CP1(:,1,1));
neta1 = length(CP1(1,:,1));
for iPatches= 1:nxi1
    for j=1:neta1         
        if CP1(iPatches,j,4)~=1
            isNURBS1 = 1;
            break;
        end
    end
    if isNURBS1
        break;
    end
end

%% Material constants

% general parameters
EYoung = 250; % 250e+3
nue = 0.0; %0.0
thickness = 0.002; % 0.001
sigma0 = 1.0; % 1e+3
prestress.voigtVector = [sigma0
                         sigma0
                         0];
density = 500; %800

% Young's modulus
parameters1.E = EYoung;

% Poisson ratio
parameters1.nue = nue;

% Thickness of the membrane
parameters1.t = thickness;

% Density of the membrane
parameters1.rho = density;

% Prestress of the membrane
parameters1.prestress = prestress;

%% GUI

% Case name
caseName = 'transientCavity_coSimulation';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver; % solve_LinearSystemGMResWithIncompleteLUPreconditioning

% Co-simulation with EMPIRE
propEmpireCoSimulation.isCoSimulation = true; % Flag on whether co-simulation with EMPIRE is assumed
propEmpireCoSimulation.isInterfaceLayer = false; % Flag on whether the matlab client is used as an interface layer
propEmpireCoSimulation.strMatlabXml = 'empireMatlab'; % Name of the xml file for Matlab

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points

int1.type = 'default';
if strcmp(int1.type,'user')
    int1.xiNGP = 6;
    int1.etaNGP = 6;
    int1.xiNGPForLoad = 6;
    int1.etaNGPForLoad = 6;
    int1.nGPForLoad = 6;
    int1.nGPError = 12;
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
% .component: 'x','y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

% Start and end time of the simulation
TStart = 0.0;
TEnd = 20.0;

%% Refinement

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

a = 1;
tp1 = a;
tq1 = a;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

noKnotsXi1 = 10;
noKnotsEta1 = noKnotsXi1;
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface...
    (p1,Xi1,q1,Eta1,CP1,noKnotsXi1,noKnotsEta1,'');

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs1 = [];

% Fix all Control Points at the left end
xisup1 = [0 0];   etasup1 = [0 1];
for dir = 1:3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end

% Fix all Control Points at the rıght end
xisup1 = [1 1];   etasup1 = [0 1];
for dir = 1:3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end

% Fix the z-displacement of all Control Points
xisup1 = [0 1];   etasup1 = [0 1];
for dir = 3
    homDOFs1 = findDofs3D...
        (homDOFs1,xisup1,etasup1,dir,CP1);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC1.noCnd = 0;

% Embedded cables
cables1.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

traction1 = 0;
NBC1.noCnd = 1;
xib1 = [0 1];   etab1 = [0 1];   dirForce1 = 'y';
NBC1.xiLoadExtension = {xib1};
NBC1.etaLoadExtension = {etab1};
NBC1.loadAmplitude = {traction1};
NBC1.loadDirection = {dirForce1};
NBC1.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC1.isFollower = false;
NBC1.isTimeDependent = false;

% Collect all the Neumann boundary conditions into an arra<y
NBC = {NBC1};

%% Create the patches and the Lagrange Multiplier fields
patch1 = fillUpPatch...
    (analysis,p1,Xi1,q1,Eta1,CP1,isNURBS1,parameters1,homDOFs1,...
    inhomDOFs1,valuesInhomDOFs1,weakDBC1,cables1,NBC1,[],[],[],[],[],int1);
BSplinePatches = {patch1};
noPatches = length(BSplinePatches);

%% Compute the load vector for the visualization of the reference configuration
for counterPatches = 1:noPatches
    BSplinePatches{counterPatches}.FGamma = ...
        zeros(3*BSplinePatches{counterPatches}.noCPs,1);
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

%% Plot reference configuration
% color = [.85098 .8549 .85882];
% connections.No = 0;
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
% (BSplinePatches,connections,color,graph,'outputEnabled');

%% Plot the multipatch geometry with the patch numbering
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,graph);

%% Output the initial geometry to be read by GiD
pathToOutput = '../../outputGiD/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4GiD(BSplinePatches,pathToOutput,[caseName '_' 'SinglePatch']);

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-4;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 100;

%% Determine the Rayleigh damping coefficients
alphaR = 17.4670;
betaR = 1.8875e-04;

%% Transient analysis properties

% Time dependence
propStrDynamics.timeDependence = 'transient';

% Time integration scheme
% propStrDynamics.method = 'explicitEuler';
propStrDynamics.method = 'bossak';
propStrDynamics.alphaB = 0; % -.1
propStrDynamics.betaB = .25; % .5
propStrDynamics.gammaB = .5; % .6
% propStrDynamics.damping.method = 'rayleigh';
% propStrDynamics.damping.computeDampMtx = @computeDampMtxRayleigh;
% propStrDynamics.damping.alpha = alphaR;
% propStrDynamics.damping.beta = betaR;

% Function handle to the computation of the matrices for the defined time
% integration scheme
if strcmp(propStrDynamics.method,'explicitEuler')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsExplicitEuler;
elseif strcmp(propStrDynamics.method,'bossak')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossak;
else
    error('Choose a time integration scheme');
end

% Function handle to the update of the discrete fields for the defined time
% integration scheme
if strcmp(propStrDynamics.method,'explicitEuler')
    propStrDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propStrDynamics.method,'bossak')
    propStrDynamics.computeUpdatedVct = ...
        @computeBossakTransientUpdatedVctAccelerationField;
else
    error('Choose a time integration scheme');
end

% Function handle to the initial conditions computation
computeInitCnds = @computeNullInitCndsIGAThinStructure;
% computeInitCnds = @computeRestartInitCndsIGAGiD;

% Starting time of the simulation
propStrDynamics.TStart = TStart;

% End time of the simulation
propStrDynamics.TEnd = TEnd;

% Number of time steps
timeStepScaling = 100; % 16
propStrDynamics.noTimeSteps = timeStepScaling*(propStrDynamics.TEnd-propStrDynamics.TStart);

% Time step size
propStrDynamics.dt = (propStrDynamics.TEnd - propStrDynamics.TStart)/propStrDynamics.noTimeSteps;
str_timeStep = num2str(propStrDynamics.dt);

%% Solve the transient problem
% writeOutput = 'undefined';
propOutput.writeOutput = @writeResults4GiD;
propOutput.writeFrequency = 5;
propOutput.saveFrequency = 5;
computeInitCnd = @computeNullInitCndsIGAThinStructure;
[dHatHistory,resHistory,BSplinePatches,minElAreaSize] = ...
    solve_IGAMembraneTransient(BSplinePatches,computeInitCnd,...
    propNLinearAnalysis,propStrDynamics,propPostproc,solve_LinearSystem,...
    propOutput,pathToOutput,[caseName '_' 'SinglePatch'],...
    propEmpireCoSimulation,'outputEnabled');
% graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
%     (BSplinePatchesNitsche,dHatHistoryNitsche(:,end),graph,'outputEnabled');
save(['data_' caseName '_' 'SinglePatch' '_dt_' str_timeStep '.mat']);
return;

%% Postprocessing
scaling = 1;
t = 5;
noTimeStep = int64((t - propStrDynamics.TStart)/propStrDynamics.dt + 2);
graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
    (BSplinePatches,scaling*dHatHistory(:,noTimeStep),graph,'outputEnabled');
az = -310; % -310
el = 40; % 40
view(az,el);
camlight(0,0);
lighting phong;
axis off;
title('');
return;

%% END
