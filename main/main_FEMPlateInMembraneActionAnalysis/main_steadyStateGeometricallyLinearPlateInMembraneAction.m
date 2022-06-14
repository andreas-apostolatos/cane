%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Plane stress analysis for a rectangular plate subject to uniform
%        pressure on its top edge
%
% Date : 19.02.2014
%
%% Preamble
clear;
clc;
close all;

%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the low order basis functions
addpath('../../basisFunctions/');

% Add all equation system solvers
addpath('../../equationSystemSolvers/');

% Add all the efficient computation functions
addpath('../../efficientComputation/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/',...
        '../../FEMPlateInMembraneActionAnalysis/errorComputation/');
    
% Include performance optimzed functions
addpath('../../efficientComputation/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMPlateInMembraneActionAnalysis/';
% caseName = 'infinitePlateWithHoleQuadrilaterals';
caseName = 'cantileverBeamPlaneStress';
% caseName = 'PlateWithAHolePlaneStress';
% caseName = 'PlateWithMultipleHolesPlaneStress';
% caseName = 'InfinitePlateWithAHolePlaneStress';
% caseName = 'unitTest_curvedPlateTipShearPlaneStress';
% caseName = 'gammaStructureMixedElementsPlaneStress';
% caseName = 'NACA2412_AoA5_CSD';

% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, ~, propGaussInt] = ...
    parse_StructuralModelFromGid(pathToCase, caseName, 'outputEnabled');

%% GUI

% On the body forces
computeBodyForces = @computeConstantVerticalStructureBodyForceVct;

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Output properties
propOutput.isOutput = true;
propOutput.writeOutputToFile = @writeOutputFEMPlateInMembraneActionToVTK;
propOutput.VTKResultFile = 'undefined';

% Choose computation of the stiffness matrix
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST;
% computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed;

% Quadrature for the stiffness matrix and the load vector of the problem
% 'default', 'user'
propIntDomain.type = 'default';
propIntDomain.noGP = 1;

% Quadrature for the L2-norm of the error
intError.type = 'user';
intError.noGP = 4;

% Linear analysis
propStrDynamics = 'undefined';

% Initialize graphics index
graph.index = 1;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Compute the load vector
t = 0;
F = computeLoadVctFEMPlateInMembraneAction ...
    (strMsh, propAnalysis, propNBC, t, propGaussInt, 'outputEnabled');

%% Visualization of the configuration
graph.index = plot_referenceConfigurationFEMPlateInMembraneAction ...
    (strMsh, propAnalysis, F, homDOFs, [], graph, 'outputEnabled');

%% Initialize solution
numNodes = length(strMsh.nodes(:,1));
numDOFs = 2*numNodes;
dHat = zeros(numDOFs,1);

%% Solve the plate in membrane action problem
[dHat, FComplete, minElSize] = solve_FEMPlateInMembraneAction ...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    propNBC, computeBodyForces, parameters, computeStiffMtxLoadVct, ...
    solve_LinearSystem, propNLinearAnalysis, propIntDomain, propOutput, ...
    caseName, pathToOutput, 'outputEnabled');

%% Postprocessing
graph.visualization.geometry = 'reference_and_current';
resultant = 'stress';
component = 'xy';
nodeIDs_active = 'undefined';
contactSegments = 'undefined';
graph.index = plot_currentConfigurationAndResultants ...
    (propAnalysis, strMsh, homDOFs, dHat, nodeIDs_active, contactSegments, ...
    parameters, resultant, component, graph);

% Compute the error in the L2-norm for the case of the plane stress
% analysis over a quarter annulus plate subject to tip shear force
if strcmp(caseName,'unitTest_curvedPlateTipShearPlaneStress')
    nodeNeumann = strMsh.nodes(propNBC.nodes(1, 1), 2:end);
    funHandle = str2func(propNBC.fctHandle(1, :));
    forceAmplitude = norm(funHandle(nodeNeumann(1, 1),nodeNeumann(1, 2),nodeNeumann(1, 3), 0));
    internalRadius = 4;
    externalRadius = 5;
    propError.resultant = 'stress';
    propError.component = 'tensor';
    errorL2 = computeRelErrorL2CurvedBeamTipShearFEMPlateInMembraneAction ...
        (strMsh, dHat, parameters, internalRadius, externalRadius, ...
        forceAmplitude, propError, intError, 'outputEnabled');
end

%% END OF THE SCRIPT
