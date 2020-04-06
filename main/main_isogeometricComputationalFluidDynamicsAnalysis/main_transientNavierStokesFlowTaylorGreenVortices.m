%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Benchmark of Taylor-Green vortices. Problem for which there exist
%        analytical solution to the Navier-Stokes equation. The problem is
%        confronted using the fully nonlinear equation system.
%
% date : 04.04.2020
%
%% Preamble
clc;
clear;

%% Includes 

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Nonlinear solvers
addpath('../../equationSystemSolvers/');

% Transient analysis
addpath('../../transientAnalysis/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
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
channelLength = 2*pi;
channelHeight = 2*pi;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

factor = 1/2;

% Control Point coordinates

% x-coordinate
CP(:,:,1) = [-factor*channelLength -factor*channelLength
             factor*channelLength factor*channelLength];
        
% y-coordinate
CP(:,:,2) = [-factor*channelHeight factor*channelHeight
             -factor*channelHeight factor*channelHeight];
        
% z-coordinate
CP(:,:,3) = [0 0
             0 0];
         
% Control Point weights
CP(:,:,4) = [1 1
             1 1];
         
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = false;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i= 1:nxi
    for j=1:neta
        if CP(i,j,4)~=1
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
parameters.nue = 1e-3;

%% UI

% Analysis type
analysis.type = 'isogeometricIncompressibleFlowAnalysis';

% Function handle to the linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Define the body-force vector
amplification = 0;
computeBodyForces = amplification*[1 0]';

% Function handle to the computation of the initial conditions
computeInitCnds = @computeInitialConditionsForTaylorGreenVorticesIGA4NSE2D;

% Integration parameters

% type: 'default', 'user'
propInt.type = 'default';
if strcmp(propInt.type,'user')
    propInt.xiNGP = 6;
    propInt.etaNGP = 3;
    propInt.xiNGPForLoad = 6;
    propInt.etaNGPForLoad = 3;
    propInt.nGPForLoad = 6;
end

% Graphics :
% __________

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
xiRef = 4;
etaRef = 4;
[Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, xiRef, etaRef, 'outputEnabled');

%% Dirichlet boundary conditions 

% Homogeneous Dirichlet Boundary Conditions
homDOFs = [];

% Inhomogeneous Dirichlet Boundary Conditions
inhomDOFs = [];
propIDBC.numCnd = 12;
propIDBC.xiExtension = zeros(propIDBC.numCnd,2);
propIDBC.etaExtension = zeros(propIDBC.numCnd,2);
propIDBC.prescribedDirection = zeros(propIDBC.numCnd,1);
propIDBC.isUniqueOnBoundary = zeros(propIDBC.numCnd,1);

% Iterate over all the boundary conditions and assign their values
for i = 1:4
    % Find the correct span for the application of the inhomogeneous DBC
    if i == 1
        xiExtension = [0 1];
        etaExtension = [0 0];
    elseif i == 2
        xiExtension = [0 1];
        etaExtension = [1 1];
    elseif i == 3
        xiExtension = [0 0];
        etaExtension = [0 1];
    elseif i == 4
        xiExtension = [1 1];
        etaExtension = [0 1]; 
    end
        
    % x-component of the velocity field
    propIDBC.xiExtension(3*i-2,:) = xiExtension;
    propIDBC.etaExtension(3*i-2,:) = etaExtension;
    propIDBC.prescribedDirection(3*i-2) = 1;
    propIDBC.isUniqueOnBoundary(3*i-2) = true;
    propIDBC.irb(3*i-2,:) = ...
        findDofs3D(inhomDOFs,propIDBC.xiExtension(3*i-2,:), ...
        propIDBC.etaExtension(3*i-2,:), propIDBC.prescribedDirection(3*i-2),CP);
    
    % y-component of the velocity field
    propIDBC.xiExtension(3*i-1,:) = xiExtension;
    propIDBC.etaExtension(3*i-1,:) = etaExtension;
    propIDBC.prescribedDirection(3*i-1) = 2;
    propIDBC.isUniqueOnBoundary(3*i-1) = true;
    propIDBC.irb(3*i-1,:) = ...
        findDofs3D(inhomDOFs,propIDBC.xiExtension(3*i-1,:), ...
        propIDBC.etaExtension(3*i-1,:), propIDBC.prescribedDirection(3*i-1),CP);
    
    % pressure field
    propIDBC.xiExtension(3*i,:) = xiExtension;
    propIDBC.etaExtension(3*i,:) = etaExtension;
    propIDBC.prescribedDirection(3*i) = 3;
    propIDBC.isUniqueOnBoundary(3*i) = true;
    propIDBC.irb(3*i,:) = ...
        findDofs3D(inhomDOFs,propIDBC.xiExtension(3*i,:), ...
        propIDBC.etaExtension(3*i,:), propIDBC.prescribedDirection(3*i),CP);
end

% The prescribed values by function pointers
propIDBC.prescribedValue = ...
    {@computeXVelocityComponentForTaylorGreenVortices2D,...
    @computeYVelocityComponentForTaylorGreenVortices2D,...
    @computePressureFieldForTaylorGreenVortices2D,...
    @computeXVelocityComponentForTaylorGreenVortices2D,...
    @computeYVelocityComponentForTaylorGreenVortices2D,...
    @computePressureFieldForTaylorGreenVortices2D,...
    @computeXVelocityComponentForTaylorGreenVortices2D,...
    @computeYVelocityComponentForTaylorGreenVortices2D,...
    @computePressureFieldForTaylorGreenVortices2D,...
    @computeXVelocityComponentForTaylorGreenVortices2D,...
    @computeYVelocityComponentForTaylorGreenVortices2D,...
    @computePressureFieldForTaylorGreenVortices2D};

% Flag on the dominance of the inhomogeneous bc's to the homogeneous
propIDBC.isDominant = false;

% Find the DOFs where inhomogeneous Dirichlet boundary conditions are
% applied
for i = 1:propIDBC.numCnd
    inhomDOFs = mergesorted(inhomDOFs,propIDBC.irb(i,:));
end
inhomDOFs = unique(inhomDOFs);

%% Neumann boundary conditions
propNBC.noCnd = 2;
propNBC.xiLoadExtension = {[1 1] [1 1]};
propNBC.etaLoadExtension = {[0 1] [0 1]};
propNBC.loadAmplitude = {0 0};
propNBC.loadDirection = {1 2};
propNBC.isFollower = [false
                      false];
propNBC.computeLoadVct = ...
    {'computeLoadVctLineIGAIncompressibleNavierStokesFlow' ...
     'computeLoadVctLineIGAIncompressibleNavierStokesFlow'};
 
%% Fill up patch
BSplinePatch = fillUpPatch ...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, ...
    inhomDOFs, [], [], [], propNBC, [], [], [], [], [], propInt);

%% Transient analysis parameters

% Select time integration scheme
% method: 'explicitEuler', 'Bossak'
propFldDynamics.method = 'Bossak';

% Computation of transient problem matrices
propFldDynamics.timeDependence = 'transient';
propFldDynamics.computeProblemMtrcsTransient = @computeProblemMtrcsBossakIGA4NSE;
propFldDynamics.computeUpdatedVct = @computeBossakTIUpdatedVctAccelerationFieldIGA4NSE;

% Parameters of the selected scheme
% alphaBeta <= .5 (for unconditional stability)
propFldDynamics.alphaBeta = -.1;

% alphaBeta + gamma >= .25 (for unconditional stability)
propFldDynamics.gamma = .5 - propFldDynamics.alphaBeta;

% The start and the end time of the simulation
propFldDynamics.TStart = 0;
propFldDynamics.TEnd = 1;

% The number of time steps
propFldDynamics.noTimeSteps = 2; % 1e1, 1e3

% The time step
propFldDynamics.dt = ...
    (propFldDynamics.TEnd - propFldDynamics.TStart)/propFldDynamics.noTimeSteps;
propFldDynamics.isAdaptive = true;

%% Nonlinear analysis parameters
propNLinearAnalysis.method = 'Newton';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 1e-9;
propNLinearAnalysis.maxIter = 50;

%% Plot reference configuration
% t = propFldDynamics.TStart;
% propGraph.index = plot_referenceConfiguration4IGAIncompressibleFlow2D ...
%     (p, q, Xi, Eta, CP, isNURBS, homDOFs, inhomDOFs, propNBC, t, ...
%     propInt, propGraph, 'outputEnabled');

%% Solve the transient Stokes Problem using the Boosak scheme
[upHistoryStokes, resHistStokes, minElSize] = solve_IGATransientFlow ...
    (analysis, BSplinePatch, computeBodyForces, computeInitCnds, ...
    @computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4StokesE2D, ...
    solve_LinearSystem, @solve_IGALinearSystem, propIDBC, ...
    propFldDynamics, propNLinearAnalysis,'outputEnabled');

%% Solve the transient Navier-Stokes Problem using the Bossak scheme
[upHistoryNavierStokes, resHistNavierStokes, minElSize] = solve_IGATransientFlow ...
    (analysis, BSplinePatch, computeBodyForces, computeInitCnds, ...
    @computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4NSE2D, ...
    solve_LinearSystem, @solve_IGANLinearSystem, propIDBC, ...
    propFldDynamics, propNLinearAnalysis,'outputEnabled');

%% Postprocessing

% Plot the analytical solution versus time
% plot_transientTaylorGreenVortices2D ...
%     (p, Xi, q, Eta, CP, isNURBS, parameters, propFldDynamics, propGraph);

% Display resultant at point over time
xi = .7;
eta = .3;
propGraph.index = plot_resultantAtPointOverTimeForTaylorGreenVorticesProblem...
    (xi, BSplinePatch.p, BSplinePatch.Xi, eta, ...
    BSplinePatch.q, BSplinePatch.Eta, BSplinePatch.CP, BSplinePatch.isNURBS, ...
    parameters, upHistoryStokes, upHistoryNavierStokes, propFldDynamics, propGraph);
legend('Analytical', 'Stokes', 'Navier-Stokes');

% Visualize the solution at the end time when using the Stokes equations
propGraph.index = plot_postprocIGAIncompressibleFlow2D ...
    (BSplinePatch, upHistoryStokes(:, end), homDOFs, inhomDOFs, ...
    zeros(length(upHistoryStokes(:, end)), 1), propGraph, ...
    'outputEnabled');
title('End solution of the Stokes equations');

% Visualize the solution at the end time when using the Navier- Stokes 
% equations
propGraph.index = plot_postprocIGAIncompressibleFlow2D ...
    (BSplinePatch, upHistoryNavierStokes(:, end), homDOFs, inhomDOFs, ...
    zeros(length(upHistoryNavierStokes(:, end)), 1), propGraph, ...
    'outputEnabled');
title('End solution of the Navier-Stokes equations');

%% END OF SCRIPT