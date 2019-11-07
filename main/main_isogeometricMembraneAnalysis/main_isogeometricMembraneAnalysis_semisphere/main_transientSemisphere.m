%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Transient analysis over a single patch semisphere whose shape was
%        obtained by a form-finding analysis.
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
% FoFiGeo = importdata('./data_FoFiSemisphere.mat');
% FoFiGeo = importdata('./data_FoFiSemisphereP2Q2NoElmnts6.mat');
% FoFiGeo = importdata('./data_FoFiPool/data_FoFiSemisphereElmnts1024p4q4FF1e-11.mat');
FoFiGeo = importdata('./data_FoFiPool/data_FoFiSemisphere_Tol_1e-11_noElmnts_25_p_2_q_2.mat');

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
parameters.E = 7.0e5;

% Poisson ratio
parameters.nue = .49;

% Thickness of the shell
parameters.t = 1.65e-4;

% Density of the shell (used only for dynamics)
parameters.rho = 1050;

% Prestress for the membrane
% parameters.prestress.computeParametricCoordinates = @(X) [X(1,1)^2 + X(2,1)^2
%                                                           atan(X(2,1)/X(1,1))];
% parameters.prestress.computeBaseVectors = @(theta1,theta2) [ cos(theta2) -sin(theta2)
%                                                              sin(theta2) cos(theta2)
%                                                              0           0];
parameters.prestress.voigtVector = [7794.5
                                    7794.5
                                    0];

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'transientSemisphere_debug';

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
graph.component = 'z';

% Function handle to writing out the results
% writeOutput = @writeResults4Carat;
writeOutput = @writeResults4GiD;
% writeOutput = 'undefined';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

% Co-simulation with Empire
propEmpireCoSimulation.isCoSimulation = false;

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
Radius = norm(squeeze(CP(end,1,1:3)));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Amplitude and frequency of the root point excitation
% g20 = 0; % .1;
% omega = 10;

% supports (Dirichlet boundary conditions)

% Homogeneous Dirichlet boundary conditions
homDOFs = [];
xiSup = {[0 0] [0 1] [1 1] [0 1]};
etaSup = {[0 1] [0 0] [0 1] [1 1]};
for iSupp = 1:length(xiSup)
    for dirSupp = 1:3
        homDOFs = findDofs3D(homDOFs,xiSup{iSupp},etaSup{iSupp},dirSupp,CP);
    end
end

% for dirSupp = 1:3
%     homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
% end
% for dirSupp = 1:3
%     homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
% end
% for dirSupp = 1:3
%     homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
% end
% xiSup = [0 0];   etaSup = [0 1];
% xiSup = [0 1];   etaSup = [0 0];
% xiSup = [1 1];   etaSup = [0 1];
% xiSup = [0 1];   etaSup = [1 1];

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];
% xiSup = [0 0];   etaSup = [0 1];
% for dirSupp = 1:3
%     inhomDOFs = findDofs3D(inhomDOFs,xiSup,etaSup,dirSupp,CP);
% end
% xiSup = [0 1];   etaSup = [0 0];
% for dirSupp = 1:3
%     inhomDOFs = findDofs3D(inhomDOFs,xiSup,etaSup,dirSupp,CP);
% end
% xiSup = [1 1];   etaSup = [0 1];
% for dirSupp = 1:3
%     inhomDOFs = findDofs3D(inhomDOFs,xiSup,etaSup,dirSupp,CP);
% end
% xiSup = [0 1];   etaSup = [1 1];
% for dirSupp = 1:3
%     inhomDOFs = findDofs3D(inhomDOFs,xiSup,etaSup,dirSupp,CP);
% end
% valuesInhomDOFs = cell({});
% for iInhomDBC = 1:length(inhomDOFs)/3
%     valuesInhomDOFs{3*iInhomDBC-2} = @(t) g20*sin(omega*t);
%     valuesInhomDOFs{3*iInhomDBC-1} = @(t) 0;
%     valuesInhomDOFs{3*iInhomDBC} = @(t) 0;
% end

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pressure Load amplitude
FAmpPressure = - 43.2;
FAmpWeight = - 9.81*parameters.rho*parameters.t;
scaling = 1e0; % 1e1
FAmpPoint = @(x,y,z,t)scaling*.018*(t<0.001);

% load (Neuman boundary conditions)
% FAmp = 0;
xib1 = [0 1];   etab1 = [0 1];   dirForcePressure = 'normal';
xib2 = [0 1];   etab2 = [0 1];   dirForceWeight = 'z';
xiP2 = xiP;    etaP2 = etaP; dirForcePoint = 'normal';
NBC.noCnd = 3;
NBC.xiLoadExtension = {xib1 xib2 xiP2};
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

%% Array on the strong Dirichlet boundary conditions
strongDBC.xiExtensionHom = xiSup;
strongDBC.etaExtensionHom = etaSup;
strongDBC = {strongDBC};

%% Compute the load vector for the visualization of the reference configuration
BSplinePatches{1}.FGamma = [];
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    BSplinePatches{1}.FGamma = zeros(3*BSplinePatches{1}.noCPs,1);
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
color = [.85098 .8549 .85882];
connections = [];
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');

%% Output the initial geometry to be read by GiD
pathToOutput = '../../../outputGiD/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4GiD(BSplinePatches,pathToOutput,caseName);

%% Write the geometry for Carat++
pathToOutput = '../../../inputCarat/';
writeOutMultipatchBSplineSurface4Carat(BSplinePatches,strongDBC,connections,pathToOutput,caseName);

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-10;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 30;

%% Output the initial geometry to be read by GiD
pathToOutput = '../../../outputGiD/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4GiD(BSplinePatches,pathToOutput,caseName);

%% Transient analysis properties

% The Rayleigh parameters adjusted through experiments for this case
alphaR = 17.4670;
betaR = 1.8875e-04;

% Time dependence
propStrDynamics.timeDependence = 'transient';

% Time integration scheme
propStrDynamics.method = 'bossak';
% propStrDynamics.method = 'explicitEuler';
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

% Starting time of the simulation
TStart = .0;
propStrDynamics.TStart = TStart;

% End time of the simulation
TEnd = .73; % .73, .3
propStrDynamics.TEnd = TEnd;

% Number of time steps
timeStepScaling = 1e+3; % 1e+4
propStrDynamics.noTimeSteps = timeStepScaling*(propStrDynamics.TEnd - propStrDynamics.TStart);

% Time step size
propStrDynamics.dt = (propStrDynamics.TEnd - propStrDynamics.TStart)/propStrDynamics.noTimeSteps;
str_timeStep = num2str(propStrDynamics.dt);

%% Solve the steady-state problem
% plot_IGANLinear = '';
% [dHat,CPHistory,resHistory,hasConverged,BSplinePatches,minElASize] = ...
%     solve_IGAMembraneNLinear(BSplinePatch,propNLinearAnalysis,...
%     solve_LinearSystem,plot_IGANLinear,graph,'outputEnabled');
% scaling = 1;
% graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
%     (BSplinePatches,scaling*dHat,graph,'outputEnabled');
% return;

%% Solve the transient problem
propOutput.writeOutput = @writeResults4GiD;
propOutput.writeFrequency = 1;
computeInitCnds = @computeNullInitCndsIGAThinStructure;
[dHatHistory,resHistory,BSplinePatches,minElAreaSize] = ...
    solve_IGAMembraneTransient...
    (BSplinePatches,computeInitCnds,propNLinearAnalysis,propStrDynamics,...
    propPostproc,solve_LinearSystem,propOutput,pathToOutput,caseName,...
    propEmpireCoSimulation,'outputEnabled');
save(['data_' caseName '_dt_' str_timeStep]);

%% Postprocessing
scaling = 1;
t = .15; % .13;
noTimeStep = int64((t - propStrDynamics.TStart)/propStrDynamics.dt + 2);
% noTimeStep = 1;
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

%% Draw the displacement of the top point against the time

% Get the parametric coordinates of the top point on the form-found
% semisphere
propNewtonRaphson.eps = 1e-9;
propNewtonRaphson.maxIt = 1e3;
% P = [0; 0; 0.075];
P = [-0.00032759; 0.0010616; 0.07498];
xi0 = .5;
eta0 = .5;
p = BSplinePatches{1}.p;
Xi = BSplinePatches{1}.Xi;
q = BSplinePatches{1}.q;
Eta = BSplinePatches{1}.Eta;
CP = BSplinePatches{1}.CP;
isNURBS = BSplinePatches{1}.isNURBS;
mxi = length(Xi);
meta = length(Eta);
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
[xi,eta,Projected,isProjected,noIter] = ...
    computeNearestPointProjectionOnBSplineSurface...
    (P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,propNewtonRaphson);

% Get the knot span indices
xiSpan = findKnotSpan(xi,Xi,nxi);
etaSpan = findKnotSpan(eta,Eta,neta);
R = computeIGABasisFunctionsAndDerivativesForSurface...
    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);

% Get the DOF numbering
DOFNumbering = BSplinePatches{1}.DOFNumbering;

% Get the local number of Control Points and DOFs
noNodesLoc = (p + 1)*(q + 1);
noDOFsLoc = 3*noNodesLoc;

% Get the element freedom table for the element that the point belongs to
EFT = zeros(1,noDOFsLoc);
r = 1;
for cpj = etaSpan - q:etaSpan
    for cpi = xiSpan - p:xiSpan
        EFT(r)   = DOFNumbering(cpi,cpj,1);
        EFT(r + 1) = DOFNumbering(cpi,cpj,2);
        EFT(r + 2) = DOFNumbering(cpi,cpj,3);
        r = r + 3;
    end
end

% Compute the B-operator matrix for the displacement field
BDisplacementsGC = zeros(3,noDOFsLoc);
k = 0;
for c = 0:q
    for b = 0:p
        % Update counter
        k = k + 1;

        % Matrix containing the basis functions
        BDisplacementsGC(1,3*k - 2) = R(k,1);
        BDisplacementsGC(2,3*k - 1) = R(k,1);
        BDisplacementsGC(3,3*k) = R(k,1);
    end
end

% Get the vertical displacement of the top point on the semisphere through
% the time
% dispTopZ = zeros(propStrDynamics.noTimeSteps + 2,1);
dispTopZ = zeros(noTimeStep + 2,1);
for iTimeSteps = 1:noTimeStep + 2 % propStrDynamics.noTimeSteps + 2
    % Get the solution at the current time step
    dHat = dHatHistory(:,iTimeSteps);
    
    % Get the solution at the top point as a vector
    dispVct = BDisplacementsGC*dHat(EFT);
    
    % Save the vertical displacement of the top point into the array
    dispTopZ(iTimeSteps,1) = dispVct(3,1);
end

% Plot the displacement over the time
% plot(1:propStrDynamics.noTimeSteps + 2,dispTopZ);
plot(1:noTimeStep + 2,dispTopZ);
grid on;
xlabel('time');
ylabel('d_z at top');

%% Read results from a Carat++ result file

% Define the path to the result file and the file name
pathToResultFile = '../../../../simulations/semisphereStatic_HSU_TUM/semisphereDynamic_IGA_tol1e-9_144ElmntsP2Q2/';
fileName = 'out.post';

% Read the results file
fstring = fileread([pathToResultFile fileName '.res']);

% Get the number of DOFs
noDOFs = 3*BSplinePatch.noCPs;
dHatCarat = zeros(noDOFs,propStrDynamics.noTimeSteps);

% Write out the geometry
pathToOutput = '../../../outputGiD/isogeometricMembraneAnalysis/';
caseName = 'transientSemisphere_Carat';
writeOutMultipatchBSplineSurface4GiD(BSplinePatches,pathToOutput,caseName);

%% Loop over the time steps
noTimeSteps = 199;
for iTimeSteps = 1:noTimeSteps
    %% Read the results from a carat res file
    block = regexp(fstring,['Result "Displacement" "Load Case" ' num2str(iTimeSteps) ' Vector OnNodes\nValues'],'split');
    block(1) = [];
    out = cell(size(block));
    for k = 1:numel(block)
        out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
        out{k} = horzcat(out{k}{:});
    end
    out = cell2mat(out);
    dispCarat = out(:,2:4);
    noDOFsCarat = 3*length(dispCarat(:,1));
    if noDOFsCarat ~= noDOFs
        error('B-Spline patch data and carat results are not in agreement');
    end
    
    %% Write the results into an array readable from matlab
    for iCPs = 1:BSplinePatch.noCPs
        dHatCarat(3*iCPs - 2,iTimeSteps) = dispCarat(iCPs,1);
        dHatCarat(3*iCPs - 1,iTimeSteps) = dispCarat(iCPs,2);
        dHatCarat(3*iCPs,iTimeSteps) = dispCarat(iCPs,3);
    end
    
    %% Write out the results to a GiD case
    writeResults4GiD(analysis,BSplinePatches,dHatCarat(:,iTimeSteps),...
        'undefined','undefined',propNLinearAnalysis,propStrDynamics,...
        propPostproc,caseName,pathToOutput,'undefined',iTimeSteps - 1);
end

%% Postprocessing

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
        for c = iEtaSpan-q-1:iEtaSpan-1 
            for b = iXiSpan-p:iXiSpan
                dHatElem(iXiSpan,iEtaSpan,xiCounter) = dHatCarat(3*(c*nxi + b)-2);
                dHatElem(iXiSpan,iEtaSpan,xiCounter + 1) = dHatCarat(3*(c*nxi + b)-1);
                dHatElem(iXiSpan,iEtaSpan,xiCounter + 2) = dHatCarat(3*(c*nxi + b));
                
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
