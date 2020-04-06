%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Solve a 2D incompressible Stokes flow for the cavity benchmark
%        case for which there is an analytical solution.
%
% date : 05.04.2020
%
%% Preamble
clc;
clear;

%% Includes 

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add general auxiliary functions
addpath('../../auxiliary/');

% Nonlinear solvers
addpath('../../equationSystemSolvers/');

% Transient analysis
addpath('../../transientAnalysis/');
    
% Add all functions related to the isogeometric Computational Fluid
% Dynamics problems
addpath('../../isogeometricComputationalFluidDynamicsAnalysis/solutionMatricesAndVectors/',...
        '../../isogeometricComputationalFluidDynamicsAnalysis/solvers/',...
        '../../isogeometricComputationalFluidDynamicsAnalysis/neumannBoundaryConditions/',...
        '../../isogeometricComputationalFluidDynamicsAnalysis/graphics/',...
        '../../isogeometricComputationalFluidDynamicsAnalysis/postProcessing/',...
        '../../isogeometricComputationalFluidDynamicsAnalysis/inhomogeneousDirichletBoundaryConditions/',...
        '../../isogeometricComputationalFluidDynamicsAnalysis/initialConditions/',...
        '../../isogeometricComputationalFluidDynamicsAnalysis/transientAnalysis/',...
        '../../isogeometricComputationalFluidDynamicsAnalysis/errorComputation/');
    
%% NURBS parameters

% Geometrical parameters
channelLength = 1;
channelHeight = 1;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinate
CP(:,:,1) = [0 0
             channelLength channelLength];
        
% y-coordinate
CP(:,:,2) = [0 channelHeight
             0 channelHeight];
        
% z-coordinate
CP(:,:,3) = [0 0
             0 0];
         
% Control Point weights
CP(:,:,4) = [1 1
             1 1];
         
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = false;
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));
for i= 1:nxi
    for j = 1:neta
        if CP(i, j, 4) ~= 1
            isNURBS = true;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% Material constants 

% Kinematic viscosity
parameters.nue = 5e-2;

% Source vector
amplification = 0;
computeBodyForces  = amplification*[1 0]';
% computeBodyForces = @bodyForcesForAnalyticalSolutionToStokesProblemInUnitSquare;

%% UI

% Analysis type
analysis.type = 'isogeometricIncompressibleFlowAnalysis';

% Function handle to the linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Integration parameters

% type: 'default', 'user'
propInt.type = 'default';
if strcmp(propInt.type,'manual')
    propInt.xiNGP = 6;
    propInt.etaNGP = 3;
    propInt.xiNGPForLoad = 6;
    propInt.etaNGPForLoad = 3;
    propInt.nGPForLoad = 6;
end
propIntError.type = 'user';
propIntError.xiNGP = 10;
propIntError.etaNGP = 10;

% Graphics

% Initialize graph index
propGraph.index = 1;

% postProcComponent: 'xVelocity', 'yVelocity', 'pressure', '2normVelocity',
% 'velocityVectorPlot'
propGraph.postProcComponent = '2normVelocity';

%% Refinement

% Degree elevation
tp = 0;
tq = 0;
[Xi, Eta, CP, p, q] = degreeElevateBSplineSurface ...
    (p, q, Xi, Eta, CP, tp, tq, 'outputEnabled');

% Knot insertion
xiRef = 10;
etaRef = 10;
[Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, xiRef, etaRef, 'outputEnabled');

%% Dirichlet boundary conditions 

% Homogeneous Dirichlet Boundary Conditions
homDOFs = [];

% No Slip condition at the lower wall of the Channel
xiSupp = [0 1];
etaSupp = [0 0];
dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);
xiSupp = [0 1];
etaSupp = [0 0];
dirSupp = 2;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% We have to constraint the pressure at one location so that we ensure
% uniqueness to the pressure space due to the existence of the grad.
xiSupp = [0 0];
etaSupp = [0 0];
dirSupp = 3;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% No Slip condition at the upper wall of the Channel
xiSupp = [0 1];
etaSupp = [1 1];
dirSupp = 2;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% No Slip condition at the left wall of the Channel
xiSupp = [0 0];
etaSupp = [0 1];
dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);
xiSupp = [0 0];
etaSupp = [0 1];
dirSupp = 2;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% No Slip condition at the right wall of the Channel
xiSupp = [1 1];
etaSupp = [0 1];
dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);
xiSupp = [1 1];
etaSupp = [0 1];
dirSupp = 2;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% Inhomogeneous Dirichlet Boundary Conditions
inhomDOFs = [];

propIDBC.numCnd = 1;
propIDBC.xiSpan = zeros(propIDBC.numCnd, 2);
propIDBC.etaSpan = zeros(propIDBC.numCnd, 2);
propIDBC.prescribedDirection = zeros(propIDBC.numCnd, 1);
propIDBC.isUniqueOnBoundary = zeros(propIDBC.numCnd, 1);

% Upper boundary is a moving wall
propIDBC.xiExtension(1, :) = [0 1];
propIDBC.etaExtension(1, :) = [1 1];
propIDBC.prescribedDirection(1) = 1;
propIDBC.isUniqueOnBoundary(1) = 1;
propIDBC.irb(1, :) = findDofs3D ...
    (inhomDOFs, propIDBC.xiExtension(1, :), propIDBC.etaExtension(1, :), ...
        propIDBC.prescribedDirection(1), CP);
propIDBC.prescribedValue = {1};

% Flag on the dominance of the inhomogeneous bc's to the homogeneous
propIDBC.isDominant = 1;

% Find the DOFs where inhomogeneous Dirichlet boundary conditions are
% applied
for i = 1:propIDBC.numCnd
    inhomDOFs = mergesorted(inhomDOFs,propIDBC.irb(i, :));
end
inhomDOFs = unique(inhomDOFs);

%% Neumann boundary conditions

% Transient Neumann boundary conditions

% Initialize the boundary conditions
propNBC.noCnd = 0; 
propNBC.xiSpan = zeros(propNBC.noCnd,2);
propNBC.etaSpan = zeros(propNBC.noCnd,2);
propNBC.loadAmplitude = zeros(propNBC.noCnd,1);
propNBC.loadDirection = zeros(propNBC.noCnd,1);

% Iterate over all the boundary conditions and assign their values
propNBC.xiSpan(1, :) = [1 1];
propNBC.etaSpan(1, :) = [0 1];
propNBC.loadAmplitude(1) = 1000;
propNBC.loadDirection(1) = 1;

% Assign the pointers to the load vector function computations
propNBC.loadVctComputation = {@computeLoadVctLineIGAIncompressibleNavierStokesFlow};

%% Fill up patch
BSplinePatch = fillUpPatch ...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, ...
    inhomDOFs, [], [], [], propNBC, [], [], [], [], [], propInt);

%% Nonlinear analysis parameters
propNLinearAnalysis.method = 'Newton';
propNLinearAnalysis.eps = 1e-9;
propNLinearAnalysis.maxIter = 50;

%% Plot reference configuration
t = 0;
propGraph.index = plot_referenceConfiguration4IGAIncompressibleFlow2D ...
    (p, q, Xi, Eta, CP, isNURBS, homDOFs, inhomDOFs, propNBC, t, ...
    propInt, propGraph, 'outputEnabled');

%% Solve the steady-state Stokes transient Problem
[up, F, ~] = solve_IGAVMSStabSteadyStateStokesE2D...
    (analysis, BSplinePatch, computeBodyForces, solve_LinearSystem, ...
    propIDBC, propNBC, propInt, '');

%% Postprocessing

% Visualize the resultant throughout the domain
plot_postprocIGAIncompressibleFlow2D ...
    (BSplinePatch, up, homDOFs, inhomDOFs, ... 
    F, propGraph, 'outputEnabled');

%% END OF SCRIPT