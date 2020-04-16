%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Transient nonlinear isogeometric Navier-Stokes flow in L-shaped
%        channel.
%
% date : 05.04.2020
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
addpath('../../functionArchive/');

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

% General geometrical parameters
H = 1;
l1 = 2;
h1 = 6;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 .5 1 1];
Eta = [0 0 1 1];

factor = 1/2;

% Control Point coordinates

% x-coordinate
CP(:,:,1) = [0  0
             l1 l1 + H
             l1 l1 + H];
        
% y-coordinate
CP(:,:,2) = [h1 h1 + H
             h1 h1 + H
             0  0];
        
% z-coordinate
CP(:,:,3) = [0 0
             0 0
             0 0];
         
% Control Point weights
CP(:,:,4) = [1 1
             1 1
             1 1];
         
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = false;
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));
for i = 1:nxi
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
parameters.nue = 1e-4;

% Source vector
amplification = 0;
computeBodyForces = amplification*[1 0]';

%% UI

% Analysis type
analysis.type = 'isogeometricIncompressibleFlowAnalysis';

% Function handle to the linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Function handle to the computation of the initial conditions
computeInitCnds = @computeNullInitialConditionsIGA4NSE2D;

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
tp = 0;  tq = 0;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Knot insertion
Rxi = [.1 .2 .3 .4 .41 .42 .43 .44 .45 .46 .47 .48 .49 .51 .52 .53 .54 ...
    .55 .56 .57 .58 .59 .6 .61 .62 .63 .64 .65 .66 .67 .68 .69 .7 .71 .72 ...
    .73 .74 .75 .76 .77 .78 .79 .8 .81 .82 .83 .84 .85 .86 .88 .9 .92 .94 .97];
Reta = [.1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7 .75 .8 .85 .9 .95];
[Xi,Eta,CP] = knotRefineBSplineSurface(p,Xi,q,Eta,CP,Rxi,Reta,'outputEnabled');

%% Dirichlet boundary conditions 

% Homogeneous Dirichlet Boundary Conditions
homDOFs = [];

% Lower wall is of no-slip condition
xiSupp = [0 1]; etaSupp = [0 0]; dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp,CP);
xiSupp = [0 1]; etaSupp = [0 0]; dirSupp = 2;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp,CP);

% Upper wall is of no-slip condition
xiSupp = [0 1]; etaSupp = [1 1]; dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);
xiSupp = [0 1]; etaSupp = [1 1]; dirSupp = 2;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% Inhomogeneous Dirichlet Boundary Conditions :
inhomDOFs = [];

propIDBC.numCnd = 1;
propIDBC.xiExtension = zeros(propIDBC.numCnd,2);
propIDBC.etaExtension = zeros(propIDBC.numCnd,2);
propIDBC.prescribedDirection = zeros(propIDBC.numCnd,1);
propIDBC.isUniqueOnBoundary = zeros(propIDBC.numCnd,1);

% Left boundary is inflow boundary
propIDBC.xiExtension(1,:) = [0 0];
propIDBC.etaExtension(1,:) = [0 1];
propIDBC.prescribedDirection(1) = 1;
propIDBC.isUniqueOnBoundary(1) = 1;
propIDBC.irb(1,:) = findDofs3D(inhomDOFs,propIDBC.xiExtension(1,:), ...
    propIDBC.etaExtension(1,:), propIDBC.prescribedDirection(1), CP);

% The prescribed values by function pointers
% IDBC.prescribedValue = {@quadraticInletDistributionForVectorTransportProblems2D};
propIDBC.prescribedValue = {6e1};

% Flag on the dominance of the inhomogeneous bc's to the homogeneous
propIDBC.isDominant = 1;

% Find the DOFs where inhomogeneous Dirichlet boundary conditions are
% applied
for i = 1:propIDBC.numCnd
    inhomDOFs = mergesorted(inhomDOFs,propIDBC.irb(i,:));
end
inhomDOFs = unique(inhomDOFs);

%% Neumann boundary conditions

% Transient Neumann boundary conditions

% Initialize the boundary conditions
propNBC.noCnd = 0; 
propNBC.xiLoadExtension = {};
propNBC.etaLoadExtension = {};
propNBC.loadAmplitude = {};
propNBC.loadDirection = {};
propNBC.isFollower = [false
                      false];
propNBC.computeLoadVct = {'computeLoadVctLineIGAIncompressibleNavierStokesFlow'};
                     
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
propFldDynamics.noTimeSteps = 1e3;

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
[upHistoryStokes, resHistStokes, ~] = solve_IGATransientFlow ...
    (analysis, BSplinePatch, computeBodyForces, computeInitCnds, ...
    @computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4StokesE2D, ...
    solve_LinearSystem, @solve_IGALinearSystem, propIDBC, ...
    propFldDynamics, propNLinearAnalysis,'outputEnabled');

%% Solve the transient Navier-Stokes Problem using the Bossak scheme
[upHistoryNavierStokes, resHistNavierStokes, ~] = solve_IGATransientFlow ...
    (analysis, BSplinePatch, computeBodyForces, computeInitCnds, ...
    @computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4NSE2D, ...
    solve_LinearSystem, @solve_IGANLinearSystem, propIDBC, ...
    propFldDynamics, propNLinearAnalysis,'outputEnabled');

%% Postprocessing

% Plot the analytical solution versus time
% plot_transientTaylorGreenVortices2D ...
%     (p, Xi, q, Eta, CP, isNURBS, parameters, propFldDynamics, propGraph);

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