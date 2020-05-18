%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Modal analysis and transient simulation of the four-point sail
%
% Date : 20.02.2016
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

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../../CAGDKernel/CAGDKernel_graphics/',...
        '../../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the isogeometric Kirchhoff-Love shell formulation
addpath('../../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../../isogeometricThinStructureAnalysis/loads/',...
        '../../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../../isogeometricThinStructureAnalysis/solvers/',...
        '../../../isogeometricThinStructureAnalysis/metrics/',...
        '../../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../../isogeometricThinStructureAnalysis/output/',...
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../../isogeometricThinStructureAnalysis/graphicsMultipatches/');
    
% Add all functions related to the transient analysis
addpath('../../../transientAnalysis/');
    
%% Load the NURBS parameters from a FoFi file

% Read the data
% FoFiGeo = importdata('./data_FoFiPool/data_FoFiFourPointSailp1q1Xi50Eta50.mat');
% FoFiGeo = importdata('./data_FoFiPool/data_FoFiFourPointSailp3q3Xi50Eta50.mat');
FoFiGeo = importdata('./data_FoFiPool/data_FoFiFourPointSailp1q1Xi10Eta10.mat');

% Global variables:
Length = FoFiGeo.Length;
Width = FoFiGeo.Width;
Height = FoFiGeo.Height;

% Patch 1 :
% _________

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
parameters.E = 8e+8;

% Poisson ratio
parameters.nue = .4;

% Thickness
parameters.t = 1e-3;

% Density (used only for dynamics)
parameters.rho = 800;

% Prestress for the membrane
sigma0 = 3e+3/parameters.t;
parameters.prestress.voigtVector = [sigma0
                                    sigma0
                                    0];

% Cable parameters
parametersCable.E = 1.6e+11;
parametersCable.radiusCS = 12e-3/2;
parametersCable.areaCS = pi*parametersCable.radiusCS^2;
parametersCable.rho = 8300;
parametersCable.prestress = 6e+4/parametersCable.areaCS;
cables.parameters = {parametersCable parametersCable parametersCable parametersCable};

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'transientFourPointSail';

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
% propOutput.writeOutput = @writeResults4Carat;
% propOutput.writeOutput = 'undefined';
propOutput.writeOutput = @writeResults4GiD;
propOutput.writeFrequency = 1;

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

% Co-simulation with Empire
propEmpireCoSimulation.isCoSimulation = false;
strMatlabXml = 'undefined';

% Start and end time of the simulation
TStart = 0.0;
TEnd = 1.0;

%% Refinement

% Degree by which to elevate
a = 0;
tp = a;
tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
scaling = 0;
refXi = ceil((Length/Width)*scaling);
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

% Ground motion
% vibrationPeriod_y = .1;
% vibrationFrequency_y = 1/vibrationPeriod_y;
vibrationPeriod_x = 1/4;
vibrationFrequency_x = 1/vibrationPeriod_x;
% vibrationPeriod_z = .4;
% vibrationFrequency_z = 1/vibrationPeriod_z;
% omegay = 2*pi*vibrationFrequency_y;
omegax = 2*pi*vibrationFrequency_x;
% omegaz = 2*pi*vibrationFrequency_z;
% logA0y = 4.0;
% A0y = 10^logA0y*1e-6;
logA0x = 2.0;
A0 = 1*1e4;
A0x = A0*10^logA0x*1e-6;
% logA0z = 1.0;
% A0z = 10^logA0z*1e-6;

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];
% xiSup = {[0 0] [0 0] [1 1] [1 1]};   etaSup = {[0 0] [1 1] [0 0] [1 1]};
% for iSupp = 1:length(xiSup)
%     for dirSupp = 1:3
%         inhomDOFs = findDofs3D(inhomDOFs,xiSup{iSupp},etaSup{iSupp},dirSupp,CP);
%     end
% end
% valuesInhomDOFs = cell({});
% for iInhomDBC = 1:length(inhomDOFs) % length(inhomDOFs)/3
%     valuesInhomDOFs{iInhomDBC} = @(t)A0x*sin(omegax*t); % 3*iInhomDBC-2
% %     valuesInhomDOFs{3*iInhomDBC-1} = @(t)A0y*sin(omegay*t);
% %     valuesInhomDOFs{3*iInhomDBC} = @(t)A0z*sin(omegaz*t);
% end

% Weak homogeneous Dirichlet boundary conditions
weakDBC.noCnd = 0;
% weakDBC.method = 'penalty'; 'nitsche';
% weakDBC.estimationStabilPrm = true;
% weakDBC.alpha = 1e3;
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
scaling = 1e0; % 1e2
FAmp = -5e2*scaling;
n = 3;
omega = 2*pi*n/(TEnd - TStart);
xib = [0 1];   etab = [0 1];   dirForce = 'z';
NBC.noCnd = 1;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
% snowLoad = {@(x,y,z,t) FAmp*t*(t<=1.0)};
snowLoad = {@(x,y,z,t) FAmp*abs(sin(omega*t))};
NBC.loadAmplitude(1,1) = snowLoad;
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

%% Compute the load vector for the visualization of the reference configuration
% for iPatches = 1:noPatches
%     BSplinePatches{iPatches}.FGamma = ...
%         zeros(3*BSplinePatches{iPatches}.noCPs,1);
% %     NBC = BSplinePatches{iPatches}.NBC;
% %     for iNBC = 1:NBC.noCnd
% %         funcHandle = str2func(NBC.computeLoadVct{iNBC});
% %         BSplinePatches{iPatches}.FGamma = funcHandle...
% %             (BSplinePatches{iPatches}.FGamma,...
% %             BSplinePatches{iPatches},...
% %             NBC.xiLoadExtension{iNBC},...
% %             NBC.etaLoadExtension{iNBC},...
% %             NBC.loadAmplitude{iNBC},...
% %             NBC.loadDirection{iNBC},...
% %             NBC.isFollower(iNBC,1),0,...
% %             BSplinePatches{iPatches}.int,'outputEnabled');
% %     end
% end

%% Plot the reference configuration for the multipatch geometry before the form-finding analysis
% color = [217 218 219]/255;
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');
% % az = -40;
% % el = 40;
% % view(az,el);
% % axis off;
% % title('');
% % camlight left;
% % lighting phong;

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-6;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 1000;

%% Output the initial geometry to be read by GiD
pathToOutput = '../../../outputGiD/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4GiD(BSplinePatches,pathToOutput,caseName);

%% Assign the Rayleigh damping coefficients

% alphaR = 17.4670;
% scaling = 1.0;
% betaR = 1.8875e-04*scaling;

% Assign two significant eigenfrequencies
freqI = 4.656892249680108;
freqJ = 8.953736617642974;

% Assign the corresponding logarithmic decrements
zetaI = .1;
zetaJ = .1;

% Compute the corresponding circular frequencies
omegaI = 2*pi*freqI;
omegaJ = 2*pi*freqJ;

% Compute the Rayleigh damping parameters
detM = omegaJ/omegaI - omegaI/omegaJ;
coef = 2/detM;
alphaR = coef*(omegaJ*zetaI - omegaI*zetaJ);
betaR = coef*(-1/omegaJ*zetaI + 1/omegaI*zetaJ);

%% Transient analysis properties

% Time dependence
propStrDynamics.timeDependence = 'transient';

% Time integration scheme
% propStrDynamics.method = 'explicitEuler';
propStrDynamics.method = 'bossak';
propStrDynamics.alphaB = 0; % -.1
propStrDynamics.betaB = .25; % .5
propStrDynamics.gammaB = .5; % .6
propStrDynamics.damping.method = 'rayleigh';
propStrDynamics.damping.computeDampMtx = @computeDampMtxRayleigh;
propStrDynamics.damping.alpha = alphaR;
propStrDynamics.damping.beta = betaR;

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
timeStepScaling = 500; % 16
propStrDynamics.noTimeSteps = timeStepScaling*(propStrDynamics.TEnd-propStrDynamics.TStart);

% Time step size
propStrDynamics.dt = (propStrDynamics.TEnd - propStrDynamics.TStart)/propStrDynamics.noTimeSteps;

%% Write the time dependent load curve for Carat corresponding to the root point excitation
% dtCarat = 0.010000000000000;
% TStartCarat = TStart;
% TEndCarat = TEnd;
% noTimeStepsCarat = (TEndCarat - TStartCarat)/dtCarat;
% timeInstance = 0.0;
% for iTimeStep = 1:noTimeStepsCarat + 1
%     functionHandle = valuesInhomDOFs{iInhomDBC};
%     fprintf('   TIME = 	%d		VAL = 	%d\n',timeInstance,functionHandle(timeInstance));
%     timeInstance = timeInstance + dtCarat;
% end

%% Write the time dependent load curve for Carat corresponding to the load
% dtCarat = 0.062500000000000;
% TStartCarat = TStart;
% TEndCarat = TEnd;
% noTimeStepsCarat = (TEndCarat - TStartCarat)/dtCarat;
% timeInstance = 0.0;
% for iTimeStep = 1:noTimeStepsCarat + 1
%     functionHandle = NBC.loadAmplitude{1};
%     fprintf('   TIME = 	%d		VAL = 	%d\n',timeInstance,functionHandle{1}(0,0,0,timeInstance));
%     timeInstance = timeInstance + dtCarat;
% end

%% Solve the transient problem
[dHatHistory,resHistory,BSplinePatches,minElAreaSize] = ...
    solve_IGAMembraneTransient ...
    (BSplinePatches, computeInitCnds, propNLinearAnalysis, propStrDynamics, ...
    propPostproc, solve_LinearSystem, propOutput, pathToOutput, caseName,...
    propEmpireCoSimulation, 'outputEnabled');
str_timeStep = num2str(propStrDynamics.dt);
save(['./data_TransientAnalysisPool/' 'data_' caseName '_dt_' str_timeStep '.mat']);
return;

%% Postprocessing

% On the postprocessing
graph.postprocConfig = 'current'; % 'reference','current','referenceCurrent'
graph.resultant = 'force'; % 'displacement','strain','curvature','force','moment','shearForce'
graph.component = '1Principal'; % 'x','y','z','2norm','1','2','12','1Principal','2Principal'

% Visualize the deformed configuration with the selected resultants
% scaling = 1;
% t = .81; % .13;
scaling = 10;
t = .36; % .13;
noTimeStep = int64((t - propStrDynamics.TStart)/propStrDynamics.dt + 2);
graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
    (BSplinePatches,scaling*dHatHistory(:,noTimeStep),graph,'outputEnabled');
az = 260;
el = 30;
view(az,el);
camlight(30,70);
lighting phong;
axis off;
title('');
h = colorbar;
limitsColorbar = get(h,'Limits');

% Write the geometry for GiD
BSplinePatchesMod = BSplinePatches;
BSplinePatchesMod = computeUpdatedGeometryIGAThinStructureMultipatches...
    (BSplinePatches,scaling*dHatHistory(:,noTimeStep));
for iPatches = 1:noPatches
    BSplinePatchesMod{iPatches}.CP = BSplinePatchesMod{iPatches}.CPd;
end
pathToOutput = '../../../outputGiD/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4GiD(BSplinePatchesMod,pathToOutput,[caseName '_' 'SinglePatch']);

%% END