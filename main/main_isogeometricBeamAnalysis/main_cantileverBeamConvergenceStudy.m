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

% On the analysis to be performed
addpath('../../equationSystemSolvers/');
    
% Add all functions related to the isogeometric beam analysis
addpath('../../isogeometricBeamAnalysis/stiffnessMatrices/',...
        '../../isogeometricBeamAnalysis/graphics/',...
        '../../isogeometricBeamAnalysis/loads/',...
        '../../isogeometricBeamAnalysis/postprocessing/',...
        '../../isogeometricBeamAnalysis/solvers/',...
        '../../isogeometricBeamAnalysis/math/',...
        '../../isogeometricBeamAnalysis/auxiliary/',...
        '../../isogeometricBeamAnalysis/errorComputation/');

%% CAD model (NURBS)

% Polynomial degree
p = 1;

% Knot vector
Xi = [0 0 1 1];

% Beam's length
L = 10;

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
parameters.EYoung = 4e6;

% Poisson ration
parameters.Nu = 0;

% shear modulus (connected to epsilon_12 = E/(1+nu))
parameters.GShear = parameters.EYoung/(2*(1 + parameters.Nu));

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

%% GUI

% Analysis type : 'Bernoulli', 'Timoshenko'
analysis.type = 'Bernoulli';

% Define linear equation solver
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
% resultant :
% 'displacement'
% 'crossSectionalRotation'
% 'force'
% 'moment'
% 'shearForce'
graph.resultant = 'moment';
% Component (only for graph.resultant = 'displacement'):
% 'x'
% 'y'
% '2norm'
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
% 0 : No outpur graph
% 1 : The basis functions themselves
% n : Up to the (n-1)-th derivative 
%
graph.plotBasisFunctionsAndDerivs = 0;

%% Initial refinement

% Order elevation
a = 1;
if strcmp(analysis.type, 'Bernoulli')
    tp = a;
elseif strcmp(analysis.type, 'Timoshenko')
    tp = a-1;
end
[Xi, CP, p] = degreeElevateBSplineCurve ...
    (p, Xi, CP, tp, 'outputEnabled');

% Knot insertion
n = 0;
[Xi, CP] = knotRefineUniformlyBSplineCurve ...
    (n, p, Xi, CP, 'outputEnabled');

%% Perform a refinement study

% Number of h-refinement steps
numRef = 50;

% Initialize arrays
errorDisplacements = zeros(numRef, 1);
errorRotations = zeros(numRef, 1);
minElSize = zeros(numRef, 1);

fprintf('---------------------------------------- \n');
fprintf('| Loop over all the h-refinement steps | \n');
fprintf('---------------------------------------- \n \n');

for i = 1:numRef
    %% Print message on the current h-refinement
    fprintf('\t Refinement step %d/%d \n', i, numRef);
    fprintf('\t ----------------------- \n \n');
    
    %% Perform an h-refinement
    [XiRef, CPRef] = knotRefineUniformlyBSplineCurve(10*i, p, Xi, CP, '');
    fprintf('\t \t No. Elements = %d \n', length(XiRef));
    if strcmp(analysis.type,'Bernoulli')
        fprintf('\t \t No. DOFs = %d \n', 2*length(CPRef(:, 1)));
    elseif strcmp(analysis.type,'Timoshenko')
        fprintf('\t \t No. DOFs = %d \n', 3*length(CPRef(:, 1)));
    end
    fprintf('\n \n');
    
    %% Fill up computational patch
    BSplinePatch.p = p;
    BSplinePatch.Xi = XiRef;
    BSplinePatch.CP = CPRef;
    BSplinePatch.isNURBS = isNURBS;
    
    %% Compute the load vector
    % On the application of a pressure load on the beam
    NBC.noCnd = 1;
    xib = [0 1];
    NBC.xiLoadExtension = {xib};
    NBC.etaLoadExtension = {'undefined'};
    pLoad = - 1e3;
    NBC.loadAmplitude = {pLoad};
    loadDir = 2;
    NBC.loadDirection = {loadDir};
    NBC.isFollower(1, 1) = false;
    NBC.isTimeDependent(1, 1) = false;
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
            NBC.etaLoadExtension{iNBC}, NBC.loadAmplitude{iNBC},...
            NBC.loadDirection{iNBC}, NBC.isFollower(iNBC,1), 0, int, '');
    end
    
    %% Find the constrained DOFs
    
    % Initialize the array of the constrained DOFs
    homDOFs = [];
    
    if strcmp(analysis.type, 'Bernoulli')
    % Clamp the left edge of the beam
    xib = [XiRef(1) XiRef(p + 1)];
    dir = 1;
    homDOFs = findDofsForBernoulliBeams2D ...
        (homDOFs, xib, dir, CPRef);
    xib = [XiRef(1) XiRef(p + 2)];
    dir = 2;
    homDOFs = findDofsForBernoulliBeams2D ...
        (homDOFs, xib, dir, CPRef);
    
    % Clamp the right edge of the beam
    %     xib = [Xi(end - p) Xi(end)];
    %     dir = 1;
    %     homDOFs = findDofsForBernoulliBeams2D(homDOFs, xib, dir, CP);
    %     xib = [Xi(length(Xi) - p - 1) Xi(length(Xi))];
    %     dir = 2;
    %     homDOFs = findDofsForBernoulliBeams2D(homDOFs, xib, dir, CP);
    elseif strcmp(analysis.type, 'Timoshenko')
        % Clamp the left edge of the beam (3 DoFs two translations and 1 rotation)
        xib = [XiRef(1) XiRef(p + 1)];
        dir = 1;
        homDOFs = findDofsForTimoshenkoBeams2D ...
            (homDOFs, xib, dir, CPRef);
        xib = [XiRef(1) XiRef(p + 1)];
        dir = 2;
        homDOFs = findDofsForTimoshenkoBeams2D ...
            (homDOFs, xib, dir, CPRef);
        xib = [XiRef(1) XiRef(p + 1)];
        dir = 3;
        homDOFs = findDofsForTimoshenkoBeams2D ...
            (homDOFs, xib, dir, CPRef);

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
    
    %% Solve the linear equation system
    [dHat, ~, minElSize(i, 1)] = solve_IGABeamLinear2D ...
        (analysis, p, XiRef, CPRef, homDOFs, NBC, parameters, isNURBS, ...
        solve_LinearSystem, int, '');
    
    %% Compute the relative error
    problemSettings.Length = L;
    problemSettings.pressure = pLoad;
    problemSettings.EYoung = parameters.EYoung;
    problemSettings.I = parameters.I;
    problemSettings.GShear = parameters.GShear;
    problemSettings.Aq = parameters.Aq;
    if strcmp(analysis.type, 'Bernoulli')
        errorDisplacements(i, 1) = ...
            computeErrIGABernoulliBeam2D ...
            (p, XiRef, CPRef, isNURBS, dHat, ...
            @computeExactDispl4BernoulliCantileverBeamInUniformPressure,...
            problemSettings, intError, '');
    elseif strcmp(analysis.type,'Timoshenko')
        [errorDisplacements(i, 1), errorRotations(i, 1)] = ...
            computeErrIGATimoshenkoBeam2D ...
            (p, XiRef, CPRef, isNURBS, dHat, ...
            @computeExactDispl4TimoshenkoCantileverBeamInUniformPressure, ...
            problemSettings, intError, '');
    end
end

%% Plot the convergence graphs
[minElSize, pid] = sort(minElSize) ;
errorDisplacements = errorDisplacements(pid);
figure(graph.index)
loglog(minElSize, errorDisplacements);
title('Refinement study over the L2-norm of the displacement');
grid on;
graph.index = graph.index + 1;
if strcmp(analysis.type, 'Timoshenko')
    errorRotations = errorRotations(pid);
    figure(graph.index)
    loglog(minElSize, errorRotations);
    title('Refinement study over the L2-norm of the rotation');
end
grid on;
graph.index = graph.index + 1;

%% End of script
