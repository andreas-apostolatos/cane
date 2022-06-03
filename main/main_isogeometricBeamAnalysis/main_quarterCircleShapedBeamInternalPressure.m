%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : The benchmark is a quarter of a circle shaped beam subject to
%        uniform internal pressure. Analytical solution suggest uniform
%        expansion of the beam without any shear deformation for the
%        Timoshenko beam theory
%
% Date : 13.11.2013
%
%% Preamble
clear;
clc;

%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add system solvers
addpath('../../equationSystemSolvers/');

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
        '../../isogeometricBeamAnalysis/auxiliary/');

%% CAD model (NURBS)

% Polynomial degree
p = 2;

% Knot vector
Xi = [0 0 0 1 1 1];

% Beam's length
L = 10;

% Control Point coordinates and weights

% x-coordinates
CP(:, 1) = [0 1 1]*L;

% y-coordinates
CP(:, 2) = [1 1 0]*L;

% z-coordinates
CP(:, 3) = [0 0 0];

% weights
CP(:, 4) = [1 1/sqrt(2) 1];

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
parameters.EYoung = 4e6;

% Poisson ration
parameters.Nu = 0;

% shear modulus (connected to epsilon_12 = E/(1+nu))
parameters.GShear = parameters.EYoung/(2*(1+parameters.Nu));

% shear correction factor
parameters.alpha = 5/6;

% width of the beam
parameters.b = 1;

% height of the beam
parameters.h = 1;

% cross sectional area
parameters.A = parameters.b*parameters.h;

% Moment of inertia (I_z = I for a simple 2D case)
parameters.I = parameters.b*(parameters.h^3)/12;

% shear cross sectional area
parameters.Aq = parameters.alpha*parameters.A;

%% UI

% Analysis type (Bernoulli or Timoshenko Beam Theory)
analysis.type = 'Timoshenko';

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration
% .type = 'default' : Default choice of Gauss points
% .type = 'user' : Manual choice of Gauss points
int.type = 'default';

% Number of Gauss points
int.noGP = 6;
% int.ngaussLoad = ceil(0.5*(p+1));
int.noGPLoad = 6;

% Postprocessing features
% resultant :
% 'displacement'
% 'crossSectionalRotation'
% 'force'
% 'moment'
% 'shearForce'
graph.resultant = 'displacement';
% Component (only for graph.resultant = 'displacement'):
% 'x'
% 'y'
% '2norm'
graph.component = '2norm';

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
% 0 : No outpur graph
% 1 : The basis functions themselves
% n : Up to the (n-1)-th derivative 
%
graph.plotBasisFunctionsAndDerivs = 0;

%% Refinement

% Order elevation
tp = 1;
[Xi, CP, p] = degreeElevateBSplineCurve ...
    (p, Xi, CP, tp, 'outputEnabled');

% Knot insertion
n = 1000;
[Xi, CP] = knotRefineUniformlyBSplineCurve ...
    (n, p, Xi, CP, 'outputEnabled');

%% Plot the shape functions and their derivatives
numEval = 149;
graph.index = plot_IGABasisFunctionsAndDerivativesForCurve ...
    (p, Xi, CP, numEval, isNURBS, graph, 'outputEnabled');
    
%% Apply boundary conditions

% Dirichlet boundary conditions
homDOFs = [];

% Get the number of Control Points
numCPsxi = length(CP(:,1,1));

% Get number of knots
numKnotsxi = numCPsxi + p + 1;
    
if strcmp(analysis.type, 'Bernoulli')
    % Clamp left edge of the beam
    xib = [Xi(1) Xi(1)];
    dir = 1;
    homDOFs = findDofsForBernoulliBeams2D ...
        (homDOFs, xib, dir, CP);
%     homDOFs(end + 1) = 4;
    
    % Clamp left edge of the beam
    xib = [Xi(numKnotsxi) Xi(numKnotsxi)];
    dir = 2;
    homDOFs = findDofsForBernoulliBeams2D ...
        (homDOFs, xib, dir, CP);
%     homDOFs(end + 1) = 2*nxi - 3;
%     homDOFs = sort(homDOFs);
    
elseif strcmp(analysis.type, 'Timoshenko')
    % Fix displacents and rotations accordingly at each end of the beam
    
    % Clamp the left edge
    xib = [Xi(1) Xi(p + 1)];
    dir = 1;
    homDOFs = findDofsForTimoshenkoBeams2D(homDOFs, xib, dir, CP);
    xib = [Xi(1) Xi(p + 1)];
    dir = 3;
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);
    
    % Clamp the right edge
    xib = [Xi(length(Xi) - p) Xi(length(Xi))];
    dir = 2;
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);
    xib = [Xi(length(Xi) - p) Xi(length(Xi))];
    dir = 3;
    homDOFs = findDofsForTimoshenkoBeams2D ...
        (homDOFs, xib, dir, CP);
end

% Neumann boundary conditions

% On the application of a pressure load on the beam
NBC.noCnd = 1;
xib = [0 1];
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {'undefined'};
pLoad = 1e5;
NBC.loadAmplitude = {pLoad};
loadDir = 2;
NBC.loadDirection = {loadDir};
NBC.isFollower(1, 1) = false;
NBC.isTimeDependent(1, 1) = false;

% On the application of a poiint load on the beam
% xi = 1;
% dir = 2;
% F = computeLoadPointVectorBeams2D(F,xi,p,Xi,CP,load,analysis,dir);

%% Fill up computational patch
BSplinePatch.p = p;
BSplinePatch.Xi = Xi;
BSplinePatch.CP = CP;
BSplinePatch.isNURBS = isNURBS;

%% Compute the load vector only for the visualization of the reference configuration
for iNBC = 1:NBC.noCnd
    F = [];
    if strcmp(analysis.type,'Bernoulli')
        NBC.computeLoadVct = {'computeLoadVctLinePressureVectorForIGABernoulliBeam2D'};
    elseif strcmp(analysis.type,'Timoshenko')
        NBC.computeLoadVct = {'computeLoadVctLinePressureVectorForIGATimoshenkoBeam2D'};
    end
    funcHandle = str2func(NBC.computeLoadVct{iNBC});
    F = funcHandle ...
        (F, BSplinePatch, NBC.xiLoadExtension{iNBC}, ...
        NBC.etaLoadExtension{iNBC}, NBC.loadAmplitude{iNBC}, ...
        NBC.loadDirection{iNBC}, NBC.isFollower(iNBC,1), 0, int, '');
end

%% Plot reference configuration
graph.index = plot_referenceConfigurationIGABeams...
    (p,Xi,CP,isNURBS,homDOFs,F,analysis,graph,'outputEnabled');

%% Solve the problem 
[dHat,F,minElEdgeSize] = solve_IGABeamLinear2D...
    (analysis,p,Xi,CP,homDOFs,NBC,parameters,isNURBS,solve_LinearSystem,...
    int,'outputEnabled');

%% Postprocessing
graph.index = plot_currentConfigurationIGABeamsLinear...
    (analysis,p,Xi,CP,isNURBS,dHat,homDOFs,parameters,graph,'outputEnabled');

%% End of script
