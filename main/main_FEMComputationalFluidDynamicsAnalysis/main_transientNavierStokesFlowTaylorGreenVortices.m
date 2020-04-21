%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Script documentation
%
% Task : Solves the transient incompressible Navier-Stokes equations
%
% Date : 06.04.2020
%
%% Preamble
clear;
clc;
close all;

%% Includes

% Add transient analysis functions
addpath('../../transientAnalysis/');

% Add functions related to equation system solvers
addpath('../../equationSystemSolvers/');

% Add general math functions
addpath('../../generalMath/');

% Add the classical finite element basis functions
addpath('../../basisFunctions/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/postProcessing/');

% Add all functions related to the Finite Element Methods for Computational
% Fluid Dynamics problems
addpath('../../FEMComputationalFluidDynamicsAnalysis/solutionMatricesAndVectors/',...
        '../../FEMComputationalFluidDynamicsAnalysis/initialConditions',...
        '../../FEMComputationalFluidDynamicsAnalysis/solvers/',...
        '../../FEMComputationalFluidDynamicsAnalysis/graphics/',...
        '../../FEMComputationalFluidDynamicsAnalysis/loads/',...
        '../../FEMComputationalFluidDynamicsAnalysis/output/',...
        '../../FEMComputationalFluidDynamicsAnalysis/ALEMotion/',...
        '../../FEMComputationalFluidDynamicsAnalysis/postProcessing/',...
        '../../FEMComputationalFluidDynamicsAnalysis/transientAnalysis/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');

%% Parse the data from the GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
%caseName = 'taylorGreenVortices_pi_domain';
%caseName = 'taylorGreenVortices_2pi_domain';
caseName = 'unitTest_taylorGreenVortices_2pi_domain';

% Parse the data
[fldMsh, homDOFs, inhomDOFs, ~, nodesALE, propNBC, ...
    propAnalysis, parameters, propNLinearAnalysis, propFldDynamics, ...
    propGaussInt] = parse_FluidModelFromGid...
    (pathToCase, caseName, 'outputEnabled');

% On the postprocessing properties
%'xVelocity','yVelocity','pressure','2normVelocity','velocityVectorPlot'
propPostproc.postProcComponent = 'pressure';
if strcmp(propPostproc.postProcComponent, 'xVelocity')
    propPostproc.computeAnalytical = ...
        @(x, y, t, parameters) -cos(x)*sin(y)*exp(-2*t*parameters.nue);
elseif strcmp(propPostproc.postProcComponent, 'yVelocity')
    propPostproc.computeAnalytical = ...
        @(x, y, t, parameters) sin(x)*cos(y)*exp(-2*t*parameters.nue);
elseif strcmp(propPostproc.postProcComponent, 'pressure')
    propPostproc.computeAnalytical = ...
        @(x, y, t, parameters) -.25*(cos(2*x) + cos(2*y))*exp(-4*t*parameters.nue);
elseif strcmp(propPostproc.postProcComponent, '2normVelocity')
    propPostproc.computeAnalytical = ...
        @(x, y, t, parameters) norm([-cos(x)*sin(y)*exp(-2*t*parameters.nue)
                                     sin(x)*cos(y)*exp(-2*t*parameters.nue)]);
end

%% UI

% On the graph
propGraph.index = 1;

% On the computation of the body forces
computeBodyForces = @computeConstantVerticalFluidBodyForceVct;

% On the writing the output function
% propVTK.isOutput = true;
% propVTK.writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;
% propVTK.VTKResultFile = 'undefined';
propVTK.isOutput = false;
propVTK.writeOutputToFile = 'undefined';
propVTK.VTKResultFile = 'undefined'; % '_contourPlots_75'

%% GUI
if strcmp(propFldDynamics.method, 'BOSSAK')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE;
end
   
%% On transient inhomogeneous Dirichlet boundary conditions
computeTaylorGreenBCs = @(fldMsh,propIDBC,t) reshape([-cos(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)).*sin(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2))*exp(-2*t*propIDBC.nue),...
                                                       sin(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)).*cos(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2))*exp(-2*t*propIDBC.nue),...
                                                      -0.25*( cos(2*fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)) + cos(2*fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2)) )*exp(-4*t*propIDBC.nue)]',1,[]);
% Assign anonymous function as a function handle
updateInhomDOFs = computeTaylorGreenBCs;
propIDBC = [];

% Assign the taylor-Green boundary conditions function
valuesInhomDOFs = computeTaylorGreenBCs(fldMsh,parameters,propFldDynamics.T0);

%% Choose the equation system solver
if strcmp(propAnalysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen')
end
   
%% Define the initial condition function
computeInitialConditions = @computeInitialConditionsForTaylorGreenVorticesFEM4NSE2D;

%% Solve the CFD problem
[upHistory, minElSize] = solve_FEMVMSStabTransientNSEBossakTI ...
    (fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, updateInhomDOFs, ...
    nodesALE, parameters, computeBodyForces, propAnalysis, ...
    computeInitialConditions, solve_LinearSystem, propFldDynamics, ...
    propNLinearAnalysis, propIDBC, propGaussInt, propVTK, caseName, ...
    'outputEnabled');

%% Visualize analytical solution
% propGraph.index = plot_transientTaylorGreenVortices2D ... 
% (fldMsh, parameters, propFldDynamics.TEnd ,propGraph, 'outputEnabled');

%% Compute the selected resultant at the chosen Cartesian location over time
x = pi/2; % pi/3, pi/2;
y = -pi/2; % pi/3, -pi/2;
[timeSpaceDiscrete, resultantNumerical, resultantAnalytical] = ...
    computeResultantAtPointOverTime...
    (x, y, fldMsh, parameters, upHistory, ...
    propFldDynamics, propPostproc, 'outputEnabled');

%% Plot the selected resultant at the chosen Cartesian location over time
figure(propGraph.index)
if ~ischar(resultantAnalytical)
    plot(timeSpaceDiscrete, resultantAnalytical, 'black',...
         timeSpaceDiscrete, resultantNumerical, 'blue');
    legend('Analytical', 'FEM', 'Orientation', 'horizontal', 'Location', 'southoutside');
else
    plot(timeSpaceDiscrete, resultantNumerical, 'blue');
end
xlabel('time (seconds)');
if strcmp(propPostproc.postProcComponent, 'xVelocity')
    yLabelString = 'x-velocity component u_x (m/s)';
elseif strcmp(propPostproc.postProcComponent, 'yVelocity')
    yLabelString = 'y-velocity component u_y ';
elseif strcmp(propPostproc.postProcComponent, 'pressure')
    yLabelString = 'pressure p (Pa)';
elseif strcmp(propPostproc.postProcComponent,'2normVelocity')
    yLabelString = '||u|| (m/s)';
end
ylabel(yLabelString);
title(sprintf('Evolution of %s at point X = (%d, %d)', ...
    propPostproc.postProcComponent, x, y));
propGraph.index = propGraph.index + 1;

%% END OF THE SCRIPT