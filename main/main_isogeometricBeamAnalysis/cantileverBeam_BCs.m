%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Timoshenko and Bernoulli models for cantilever beams
%
% Date : 12.05.2021
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
analysis.type = 'Bernoulli'; % 'Bernoulli', 'Timoshenko'
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
graph.resultant = 'moment';
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
tp = 1;
[Xi, CP, p] = degreeElevateBSplineCurve ...
    (p, Xi, CP, tp, 'outputEnabled');

% Knot insertion
n = 20; % No. elements
[Xi, CP] = knotRefineUniformlyBSplineCurve ...
    (n, p, Xi, CP, 'outputEnabled');
% Rxi = [.25 .5 .5];
% [Xi, CP] = knotRefineBSplineCurve ...
%     (p,Xi,CP,Rxi,'outputEnabled');

%% Plot the shape functions and their derivatives
% numEval = 149;
% graph.index = plot_IGABasisFunctionsAndDerivativesForCurve ...
%     (p, Xi, CP, numEval, isNURBS, graph, 'outputEnabled');

%% Apply boundary conditions

% Dirichlet boundary conditions
homDOFs = [];
if strcmp(analysis.type, 'Bernoulli') % For the Bernoulli beam
    % Fix the displacements at the left edge of the beam
    xib = [Xi(1) Xi(p + 1)];
    dir = 1;
    homDOFs = findDofsForBernoulliBeams2D ...
        (homDOFs, xib, dir, CP);
    
    % Fix the rotations at the left edge of the beam: Note that for
    % clamping to be effective, also the displacements need to be fixed
    xib = [Xi(1) Xi(p + 2)];
    dir = 2;
    homDOFs = findDofsForBernoulliBeams2D ...
        (homDOFs, xib, dir, CP);
    
    % Fix the displacements at the right edge of the beam
    xib = [Xi(length(Xi) - p) Xi(end)]; 
    dir = 1;
    
    % Fix the rotations at the right edge of the beam: Note that for
    % clamping to be effective, also the displacements need to be fixed
    homDOFs = findDofsForBernoulliBeams2D(homDOFs, xib, dir, CP);
    xib = [Xi(length(Xi)-p-1) Xi(length(Xi))];
    dir = 2;
    homDOFs = findDofsForBernoulliBeams2D(homDOFs, xib, dir, CP);
elseif strcmp(analysis.type, 'Timoshenko') % For the Timoshenko beam
    % Fix the displacements at the left edge of the beam
    xib = [Xi(1) Xi(p + 1)];
    dir = 1; % x-component of the displacement
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);
    xib = [Xi(1) Xi(p + 1)];
    dir = 2; % y-component of the displacement
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);
    xib = [Xi(1) Xi(p + 1)];
    
    % Fix the rotations at the left edge of the beam: Note that for
    % clamping to be effective, also the displacements need to be fixed
    dir = 3; % rotation
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);

    % Fix the displacements at the right edge of the beam
    xib = [Xi(end - p) Xi(end)];
    dir = 1; % x-component of the displacement
    homDOFs = findDofsForTimoshenkoBeams2D(homDOFs, xib, dir, CP);
    xib = [Xi(end - p) Xi(end)];
    dir = 2; % y-component of the displacement
    homDOFs = findDofsForTimoshenkoBeams2D(homDOFs, xib, dir, CP);
    xib = [Xi(end - p) Xi(end)];
    
    % Fix the rotations at the right edge of the beam: Note that for
    % clamping to be effective, also the displacements need to be fixed
    dir = 3; % rotation
    homDOFs = findDofsForTimoshenkoBeams2D(homDOFs, xib, dir, CP);
end

% Neumann boundary conditions

% On the application of a pressure load on the beam
NBC.noCnd = 1;
xib = [0 1];
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {'undefined'};
% pLoad = - 1e2;
pLoad = @(xi) - 1e4*(xi - 0.5)^3; % Function handle
NBC.loadAmplitude = {pLoad};
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
    if strcmp(analysis.type, 'Bernoulli')
        NBC.computeLoadVct = ...
            {'computeLoadVctLinePressureVectorForIGABernoulliBeam2D'};
    elseif strcmp(analysis.type, 'Timoshenko')
        NBC.computeLoadVct = {'computeLoadVctLinePressureVectorForIGATimoshenkoBeam2D'};
    end
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
[dHat, Fcomplete, minElEdgeSize] = solve_IGABeamLinear2D ...
    (analysis, p, Xi, CP, homDOFs, NBC, parameters, isNURBS, ...
    solve_LinearSystem, int, 'outputEnabled');

%% Postprocessing
graph.index = plot_currentConfigurationIGABeamsLinear ...
    (analysis, p, Xi, CP, isNURBS, dHat, homDOFs, parameters, graph, ...
    'outputEnabled');

%% End of script