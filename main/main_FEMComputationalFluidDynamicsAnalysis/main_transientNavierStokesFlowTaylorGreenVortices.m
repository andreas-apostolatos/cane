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

%% Parse the data from the GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
%caseName = 'TaylorGreenVortices';
caseName = 'TaylorGreenVortices2';

% Parse the data
[fldMsh, homDOFs, inhomDOFs, ~, nodesALE, propNBC, ...
    propAnalysis, parameters, propNLinearAnalysis, propFldDynamics, ...
    propGaussInt] = parse_FluidModelFromGid...
    (pathToCase, caseName, 'outputEnabled');

%% UI

% On the computation of the body forces
computeBodyForces = @computeConstantVerticalFluidBodyForceVct;

% On the writing the output function
propVTK.isOutput = true;
propVTK.writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;
propVTK.VTKResultFile = 'undefined'; % '_contourPlots_75'

%% GUI
if strcmp(propFldDynamics.method, 'BOSSAK')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE;
end
   
%% On transient inhomogeneous Dirichlet boundary conditions

% Anonymous function to Apply Taylor-Green boundary conditions
% computeTaylorGreenBCs = @(propIDBC,t) reshape([-cos(propIDBC.coordsNode(:,1)).*sin(propIDBC.coordsNode(:,2))*exp(-2*t*propIDBC.nue),...
%                                                 sin(propIDBC.coordsNode(:,1)).*cos(propIDBC.coordsNode(:,2))*exp(-2*t*propIDBC.nue),...
%                                                -0.25*( cos(2*propIDBC.coordsNode(:,1)) + cos(2*propIDBC.coordsNode(:,2)) )*exp(-4*t*propIDBC.nue)]',1,[]);

computeTaylorGreenBCs = @(fldMsh,propIDBC,t) reshape([-cos(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)).*sin(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2))*exp(-2*t*propIDBC.nue),...
                                                       sin(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)).*cos(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2))*exp(-2*t*propIDBC.nue),...
                                                      -0.25*( cos(2*fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)) + cos(2*fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2)) )*exp(-4*t*propIDBC.nue)]',1,[]);

%% Find the mesh nodes where inhomDBCs apply

% % Number of nodes where inhomDBC apply
% noInhomNodes = length(inhomDOFs)/propAnalysis.noFields;
%                                    
% % Initialize vector of nodes and coordinates
% inhomNodes = zeros(noInhomNodes,1);
% coordsNode = zeros(noInhomNodes,2);
% 
% % Find all the coresponding nodes
% for n = 1:length(inhomDOFs)
%     % Find the corresponding node
%     indexDOF = inhomDOFs(n);
%     inhomNodes(n) = ceil(indexDOF/propAnalysis.noFields);
% end
% inhomNodes = unique(inhomNodes);
% 
% % Find correct coordinates
% count  = 1;
% for m = 1:noInhomNodes
%     nodeIndex = inhomNodes(m);
%     coordsNode(m,:) = fldMsh.nodes(nodeIndex,1:2);
%     
%     % Compute node coordinates for test solution
%     x = fldMsh.nodes(nodeIndex,1);
%     y = fldMsh.nodes(nodeIndex,2);
%     
%     % Compute test solution
%     t = propFldDynamics.T0;
%     expectedSolution(count) = -exp(-2*parameters.nue*t)*cos(x)*sin(y);
%     expectedSolution(count+1) = exp(-2*parameters.nue*t)*sin(x)*cos(y);
%     expectedSolution(count+2) = -0.25*( (cos(2*x)+cos(2*y))*exp(-4*parameters.nue*t) );
%     count = count + 3;
% end

%% On transient inhomogeneous Dirichlet boundary conditions
updateInhomDOFs = computeTaylorGreenBCs;
propIDBC = [];

%% Define the update boundary conditions function
%valuesInhomDOFs = computeTaylorGreenBCs(propIDBC,propFldDynamics.T0);
valuesInhomDOFs = computeTaylorGreenBCs(fldMsh,parameters,propFldDynamics.T0);

% For testing purpose -> if arrays are the same the value must be 1
% validate = min(abs(expectedSolution) - abs(valuesInhomDOFs) < 1e-10);

%% Choose the equation system solver
if strcmp(propAnalysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen')
end
   
%% Define the initial condition function
computeInitialConditions = @computeNullInitialConditionsFEM4NSE;
% computeInitialConditions = @computeInitialConditionsFromVTKFileFEM4NSE;

%% Solve the CFD problem
[upHistory, minElSize] = solve_FEMVMSStabTransientNSEBossakTI ...
    (fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, updateInhomDOFs, ...
    nodesALE, parameters, computeBodyForces, propAnalysis, ...
    computeInitialConditions, solve_LinearSystem, propFldDynamics, ...
    propNLinearAnalysis, propIDBC, propGaussInt, propVTK, caseName, ...
    'outputEnabled');

%% END OF THE SCRIPT
