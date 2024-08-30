%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : The benchmark is a cantilever beam subject to uniform pressure 
%        load in 2D analysis. For both the Bernoulli and the Timoshenko 
%        settings there exist an analytical solution in terms of the 
%        displacements (and cross sectional rotations for the Timoshenko 
%        problem), namely:
%
%        Bernoulli :
%             w(x) = p*x^2*(6*L^2 - 4L*x + x^2)/24/E/I
%
%       Timoshenko :
%             w(x) = (p*L*x-p*x^2/2)/G/Aq - (-p*x^4/24+p*L*x^3/6-p*L^2*x^2/4)/E/I
%          beta(x) = - p*x^3/6/E/I
%               Aq = alpha*A
%
% Date : 13.11.2013
%
%% Preamble
clear;
clc;

%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/');
    
% Add all functions related to the isogeometric beam analysis
addpath('../../isogeometricBeamAnalysis/stiffnessMatrices/',...
        '../../isogeometricBeamAnalysis/graphics/',...
        '../../isogeometricBeamAnalysis/loads/',...
        '../../isogeometricBeamAnalysis/postprocessing/',...
        '../../isogeometricBeamAnalysis/solvers/',...
        '../../isogeometricBeamAnalysis/math/',...
        '../../isogeometricBeamAnalysis/auxiliary/',...
        '../../isogeometricBeamAnalysis/errorComputation/');
    
% Add all functions related to the equation system solvers
addpath('../../equationSystemSolvers/');

%% CAD model (NURBS)

% Polynomial degree
p = 1;

% Knot vector
Xi = [0 0 1 1];

% Beam's length
L = 4;

% Control Point coordinates and weights

% x-coordinates
CP(:,1) = [0 1]*L;

% y-coordinates
CP(:,2) = [0 0]*L;      

% z-coordinates
CP(:,3) = [0 0];

% weights
CP(:,4) = [1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = false;
for i = length(CP(:, 1))
    if CP(i, 4) ~= 1
        isNURBS = true;
        break;
    end
end

%% Material constants

% Young's modulus (connected to epsilon_11 = E)
parameters.EYoung = 1e7;

% Poisson ratio
parameters.Nu = .3;

% shear modulus (connected to epsilon_12 = E/(1+nu))
parameters.GShear = parameters.EYoung/(2*(1 + parameters.Nu));

% shear correction factor
parameters.alpha = 5/6;

% Width of the beam
parameters.b = .1;

% height of the beam
parameters.h = .1;

% cross sectional area
parameters.A = parameters.b*parameters.h;

% Moment of inertia (I_z = I for a simple 2D case)
parameters.I = parameters.b*(parameters.h^3)/12;

% shear cross sectional area
parameters.Aq = parameters.alpha*parameters.A;

%% GUI

% Analysis type (Bernoulli or Timoshenko Beam Theory)
analysis.type = 'Bernoulli';
if ~strcmp(analysis.type, 'Bernoulli') && ...
        ~strcmp(analysis.type, 'Timoshenko')
    error('Choose valid analysis type')
end

% Choose linear equation solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration
% .type = 'default' : Default choice of Gauss points
% .type = 'user' : Manual choice of Gauss points
int.type = 'default';
intError.type = 'user';

% Number of Gauss points
int.noGP = 6;
int.noGPLoad = 6;
intError.noGP = 12;

% Postprocessing features
% resultant : 'displacement', 'crossSectionalRotation', 'force', 'moment',
% 'shearForce'
graph.resultant = 'shearForce';
% Component (only for graph.resultant = 'displacement'): 'x', 'y', '2norm'
graph.component = 'y';

% Plot initial and/or deformed geometry
% type 0 : None
% type 1 : Undeformed beam
% type 2 : Deformed beam
% type 3 : Undeformed && deformed beam
graph.postProcVisType = 'deformedAndUndeformedShape';
% initialization of the plot index
graph.index = 1;

% Plot NURBS basis functions and its derivatives after the refinement:
% graph.plotBasisFunctionsAndDerivs :
%
% 0 : No output graph
% 1 : The basis functions themselves
% n : Up to the (n-1)-th derivative 
%
graph.plotBasisFunctionsAndDerivs = 1;

%% Refinement

% Order elevation
tp = 2;
[Xi, CP, p] = degreeElevateBSplineCurve ...
    (p, Xi, CP, tp, 'outputEnabled');

% Knot insertion
n = 10;
[Xi, CP] = knotRefineUniformlyBSplineCurve ...
    (n, p, Xi, CP, 'outputEnabled');
% Rxi = [.25 .5 .5];
% [Xi, CP] = knotRefineBSplineCurve ...
%     (p,Xi,CP,Rxi,'outputEnabled');

%% Plot the shape functions and their derivatives
numEval = 149;
graph.index = plot_IGABasisFunctionsAndDerivativesForCurve ...
    (p, Xi, CP, numEval, isNURBS, graph, 'outputEnabled');

%% Apply boundary conditions

% Dirichlet boundary conditions
homDOFs = [];
    
if strcmp(analysis.type, 'Bernoulli')
    % Clamp the left edge of the beam
    xib = [Xi(1) Xi(p + 1)];
    dir = 1;
    homDOFs = findDofsForBernoulliBeams2D ...
        (homDOFs, xib, dir, CP);
    xib = [Xi(1) Xi(p + 2)];
    dir = 2;
    homDOFs = findDofsForBernoulliBeams2D ...
        (homDOFs, xib, dir, CP);
    
    % Clamp the right edge of the beam
%     xib = [Xi(length(Xi) - p) Xi(end)]; 
%     dir = 1;
%     homDOFs = findDofsForBernoulliBeams2D(homDOFs, xib, dir, CP);
%     xib = [Xi(length(Xi)-p-1) Xi(length(Xi))];
%     dir = 2;
%     homDOFs = findDofsForBernoulliBeams2D(homDOFs, xib, dir, CP);
elseif strcmp(analysis.type,'Timoshenko')
    % Clamp the left edge of the beam (3 DoFs two translations and 1 rotation)
    xib = [Xi(1) Xi(p + 1)];
    dir = 1;
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);
    xib = [Xi(1) Xi(p + 1)];
    dir = 2;
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);
    xib = [Xi(1) Xi(p + 1)];
    dir = 3;
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);

    % Clamp the right edge of the beam (3 DoFs two translations and 1 rotation)
%     xib = [Xi(end - p) Xi(end)];
%     dir = 1;
%     homDOFs = findDofsForTimoshenkoBeams2D(homDOFs, xib, dir, CP);
%     xib = [Xi(end - p) Xi(end)];
%     dir = 2;
%     homDOFs = findDofsForTimoshenkoBeams2D(homDOFs, xib, dir, CP);
%     xib = [Xi(end - p) Xi(end)];
%     dir = 3;
%     homDOFs = findDofsForTimoshenkoBeams2D(homDOFs, xib, dir, CP);
end

% Neumann boundary conditions

% On the application of a pressure load on the beam
NBC.noCnd = 1;
xib = [0 1]; % Distributed load
% xib = 1; % Point load
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {'undefined'};
pLoad = - 1e2; % Distributed load
% pLoad = + 1e2; % Point load
NBC.loadAmplitude = {pLoad};
if strcmp(analysis.type, 'Bernoulli') % Distributed load
    NBC.computeLoadVct = ...
        {'computeLoadVctLinePressureVectorForIGABernoulliBeam2D'};
elseif strcmp(analysis.type, 'Timoshenko')
    NBC.computeLoadVct = {'computeLoadVctLinePressureVectorForIGATimoshenkoBeam2D'};
end
% if strcmp(analysis.type, 'Bernoulli') % Point load
%     NBC.computeLoadVct = ...
%         {'computeLoadPointVectorForIGABernoulliBeam2D'};
% elseif strcmp(analysis.type, 'Timoshenko')
%     NBC.computeLoadVct = {'computeLoadPointVectorForIGATimoshenkoBeam2D'};
% end
loadDir = 2;
NBC.loadDirection = {loadDir};
NBC.isFollower(1, 1) = false;
NBC.isTimeDependent(1, 1) = false;

% On the application of a point load on the beam
% xib = 1;
% dir = 2;
% F = computeLoadPointVectorForIGATimoshenkoBeam2D ...
%     (F, BSplinePatch, xib, etab, pLoad, dir, isFollower, int, 'outputEnabled');

%% Fill up computational patch
BSplinePatch.p = p;
BSplinePatch.Xi = Xi;
BSplinePatch.CP = CP;
BSplinePatch.isNURBS = isNURBS;

%% Compute the load vector only for the visualization of the reference configuration
for iNBC = 1:NBC.noCnd
    F = [];
    funcHandle = str2func(NBC.computeLoadVct{iNBC});
    F = funcHandle ...
        (F, BSplinePatch, NBC.xiLoadExtension{iNBC}, ...
        NBC.etaLoadExtension{iNBC}, NBC.loadAmplitude{iNBC}, ...
        NBC.loadDirection{iNBC}, NBC.isFollower(iNBC,1), 0, int, '');
end

%% Plot reference configuration
graph.index = plot_referenceConfigurationIGABeams ...
    (p, Xi, CP, isNURBS, homDOFs, F, analysis, graph, 'outputEnabled');

%% Solve the problem 
[dHat, F, minElEdgeSize] = solve_IGABeamLinear2D ...
    (analysis, p, Xi, CP, homDOFs, NBC, parameters, isNURBS, ...
    solve_LinearSystem, int, 'outputEnabled');

%% Postprocessing
graph.index = plot_currentConfigurationIGABeamsLinear ...
    (analysis, p, Xi, CP, isNURBS, dHat, homDOFs, parameters, graph, ...
    'outputEnabled');

%% Benchmarking the solution
problemSettings.Length = L;
problemSettings.pressure = pLoad;
problemSettings.EYoung = parameters.EYoung;
problemSettings.I = parameters.I;
problemSettings.GShear = parameters.GShear;
problemSettings.Aq = parameters.Aq;
    
if strcmp(analysis.type,'Bernoulli')
    [relErrL2,minElSize] = computeErrIGABernoulliBeam2D(p,Xi,CP,isNURBS,dHat,...
        @computeExactDispl4BernoulliCantileverBeamInUniformPressure,...
        problemSettings,intError,'outputEnabled');
    relativeErrorL2NormDispBernoulli = relErrL2
elseif strcmp(analysis.type,'Timoshenko')
    [relDisplErrL2,relRotErrL2] = computeErrIGATimoshenkoBeam2D...
            (p,Xi,CP,isNURBS,dHat,...
        @computeExactDispl4TimoshenkoCantileverBeamInUniformPressure,...
        problemSettings,intError,'outputEnabled');
    relativeErrorL2NormDispTimoshenko = relDisplErrL2
    relativeErrorL2NormRotTimoshenko = relRotErrL2
end

if strcmp(analysis.type,'Bernoulli')
    % Compute the reference deflection value on the tip
    x = L;
    wReference = abs(pLoad)*x^2*(6*L^2 - 4*L*x + x^2)/24/parameters.EYoung/parameters.I;
    
    % Compute the numerical deflection value on the tip
    wIGA = abs(dHat(length(dHat)));
    
    % Compute the relative error in the absolute value for the tip
    % deflection
    relativeErrorTipDeflectionBernoulli = abs(wIGA-wReference)/abs(wReference)
elseif strcmp(analysis.type,'Timoshenko')
    % Compute the reference deflection values on the tip
    x = L;
    wReference = (abs(pLoad)*L*x-abs(pLoad)*x^2/2)/parameters.GShear/parameters.Aq -...
        (-abs(pLoad)*x^4/24+abs(pLoad)*L*x^3/6-abs(pLoad)*L^2*x^2/4)/parameters.EYoung/parameters.I;
    betaReference = -abs(pLoad)*x^3/6/parameters.EYoung/parameters.I;
    
    % Compute the numerical values on the tip
    wIGA = abs(dHat(length(dHat)-1));
    betaIGA = dHat(length(dHat));
    
    % Compute the relative error in the absolute value for the tip
    % resultants
    relativeErrorTipDeflectionTimoshenko = abs(wIGA-wReference)/abs(wReference)
    relativeErrorTipCrossSectionalRotationTimoshenko = 100*abs(betaIGA-betaReference)/abs(betaReference)
end

if strcmp(analysis.type,'Timoshenko')
    %             w(x) = (q*L*x-q*x^2/2)/G/Aq - (-q*x^4/24+q*L*x^3/6-q*L^2*x^2/4)/E/I
    %          beta(x) = - q*x^3/6/E/I
    %               Aq = alpha*A
    syms w(x) beta(x) Lsym qsym Gsym Aqsym Esym Isym

    w(x) = (qsym*Lsym*x-qsym*x^2/2)/Gsym/(Aqsym) - (-qsym*x^4/24+qsym*Lsym*x^3/6-qsym*Lsym^2*x^2/4)/Esym/Isym;
    beta(x) = - qsym*x^3/6/Esym/Isym;
    Q(x) = Aqsym*Gsym*(-beta(x) + diff(w, x));
    
    w(x) = subs(w,{Lsym, qsym, Gsym, Esym, Isym, Aqsym}, ...
        {L, abs(pLoad), problemSettings.GShear, problemSettings.EYoung, problemSettings.I, problemSettings.Aq});
    beta(x) = subs(beta,{Lsym, qsym, Gsym, Esym, Isym, Aqsym}, ...
        {L, abs(pLoad), problemSettings.GShear, problemSettings.EYoung, problemSettings.I, problemSettings.Aq});
    Q(x) = subs(Q,{Lsym, qsym, Gsym, Esym, Isym, Aqsym}, ...
        {L, abs(pLoad), problemSettings.GShear, problemSettings.EYoung, problemSettings.I, problemSettings.Aq});

    figure(graph.index)
    tiledlayout(3,1)
    nexttile
    fplot(w(x), [0 L]);
    title("Analytical solution for the vertical deflection field")
    nexttile
    fplot(beta(x), [0 L]);
    title("Analytical solution for the rotation field")
    nexttile
    fplot(Q(x), [0 L]);
    title("Analytical solution for the shear force field")
    graph.index = graph.index + 1;
end

%% End of script
