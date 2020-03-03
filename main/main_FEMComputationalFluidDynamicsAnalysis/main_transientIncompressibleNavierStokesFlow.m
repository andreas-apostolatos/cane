%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Solves the transient incompressible Navier-Stokes equations
%
% Date : 21.04.2014
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
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/');

% Add all functions related to the Finite Element Methods for Computational
% Fluid Dynamics problems
addpath('../../FEMComputationalFluidDynamicsAnalysis/solutionMatricesAndVectors/',...
        '../../FEMComputationalFluidDynamicsAnalysis/initialConditions',...
        '../../FEMComputationalFluidDynamicsAnalysis/solvers/',...
        '../../FEMComputationalFluidDynamicsAnalysis/loads/',...
        '../../FEMComputationalFluidDynamicsAnalysis/output/',...
        '../../FEMComputationalFluidDynamicsAnalysis/ALEMotion/',...
        '../../FEMComputationalFluidDynamicsAnalysis/transientAnalysis/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the efficient computation functions
addpath('../../efficientComputation/');

%% GUI

% On the body forces
computeBodyForces = @computeConstantVerticalFluidBodyForceVct;

% On the writing the output function
writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
% caseName = 'flowAroundCylinder';
% caseName = 'LShapedChannel';
% caseName = 'channelFlow';
% caseName = 'flowAroundSquareObstacle';
% caseName = 'flowAroundCylinderAdaptive'; % best case I have
% caseName = 'flowAroundCylinderAdaptiveFine';
% caseName = 'BenchmarkHigStrRefined';
caseName = 'flowAroundCylinderAdaptiveALE';
% caseName = 'NACA2412_AoA5_CFD';
% caseName = 'flowAroundCylinder3D'; % need to find the case
% caseName = 'unitTest_semisphere';
% caseName = 'semisphereEl150000';
% caseName = 'squareObstacleInFlow';
% caseName = 'flowAroundSquareObjectBoundaryLayerPowerLaw'; % problemZero, needs then ALE module

%% Parse the data from the GiD input file
[fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,NBC,analysis,parameters,...
    propNLinearAnalysis,propFldDynamics,gaussInt] = ...
    parse_FluidModelFromGid...
    (pathToCase,caseName,'outputEnabled');

%% Apply a non-constant inlet
if strcmp(caseName,'flowAroundSquareObjectBoundaryLayerPowerLaw')
    % Define number of DOFs per node
    noDOFsPerNode = 3;

    % Define the corresponding law
    u_max = 10;
    law = @(x,y,z) [u_max*y^(1/7)
                    0
                    0
                    0];

    % Loop over all inlet DOFs
    for i = 1:length(inhomDBC)
        % Find the inlet DOF
        indexDOF = inhomDBC(1,i);

        % Find the corresponding node
        indexNode = ceil(indexDOF/noDOFsPerNode);

        % Get the nodal coordinates
        coordsNode = fldMsh.nodes(indexNode,:);

        % Compute the value according to the law
        presValue = law(coordsNode(1,1),coordsNode(1,2),coordsNode(1,3));

        % Cartesian direction
        cartDir = indexDOF - (noDOFsPerNode*ceil(indexDOF/noDOFsPerNode) - noDOFsPerNode);

        if cartDir == 1
            valuesInhomDBC(1,i) = presValue(1,1);
        end
    end
end

%% GUI
if strcmp(propFldDynamics.method,'BOSSAK')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE2D;
end

%% Choose the equation system solver
if strcmp(analysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen')
end
   
%% Define the initial condition function
computeInitialConditions = @computeNullInitialConditionsFEM4NSE;

%% Define the name of the vtk file from where to resume the simulation
% VTKResultFile = '_contourPlots_75';
VTKResultFile = 'undefined';

%% Solve the CFD problem
[upHistory,minElSize] = solve_FEMVMSStabTransientNSEBossakTI2D...
    (fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,parameters,...
    computeBodyForces,analysis,computeInitialConditions,...
    VTKResultFile,solve_LinearSystem,propFldDynamics,...
    propNLinearAnalysis,writeOutputToFile,gaussInt,caseName,...
    'outputEnabled');

%% END OF THE SCRIPT
