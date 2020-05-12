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
% Date : 11.05.2020
%
%% Preamble
clear;
clc;
close all;

%% Includes

% path prefix
path_prefix = '../../';

% Add transient analysis functions
addpath([path_prefix 'transientAnalysis/']);

% Add functions related to equation system solvers
addpath([path_prefix 'equationSystemSolvers/']);

% Add general math functions
addpath([path_prefix 'generalMath/']);

% Add the classical finite element basis functions
addpath([path_prefix 'basisFunctions/']);

% Add all functions related to plate in membrane action analysis
addpath([path_prefix 'FEMPlateInMembraneActionAnalysis/solvers/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/loads/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/graphics/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/output/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/postprocessing/']);

% Add all functions related to the Finite Element Methods for Computational
% Fluid Dynamics problems
addpath([path_prefix 'FEMComputationalFluidDynamicsAnalysis/solutionMatricesAndVectors/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/initialConditions'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/solvers/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/loads/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/ALEMotion/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/transientAnalysis/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/postProcessing/']);

% Add all functions related to Fluid-Structure interaction
addpath([path_prefix 'FEMComputationalFluidStructureInteractionAnalysis/ALEMotion'], ...
        [path_prefix 'FEMComputationalFluidStructureInteractionAnalysis/solvers']);
    
% Add all functions related to parsing
addpath([path_prefix 'parsers/']);

% Add all functions related to the efficient computation functions
addpath([path_prefix 'efficientComputation/']);

%% Parse the data from the GiD input file
pathToCase = '../../inputGiD/FEMComputationalFluidStructureInteraction/';
caseName = 'testTurek';

% Parse fluid model from GiD
[fldMsh, homDOFsFld, inhomDOFsFld, valuesInhomDOFsFld, propALE, ~, ...
    propAnalysisFld, propParametersFld, propNLinearAnalysisFld, ...
    propFldDynamics, propGaussIntFld, propPostProcFld, propFldFSI] = ...
    parse_FluidModelFromGid ...
    (pathToCase, caseName, 'outputEnabled');

% Parse structural model from GiD
[strMsh, homDOFsStr, inhomDOFsStr, valuesInhomDOFsStr, propNBCStr, ...
    propAnalysisStr, propParametersStr, propNLinearAnalysisStr, ...
    propStrDynamics, propGaussIntStr, ~, propStrFSI] = ...
    parse_StructuralModelFromGid...
    (pathToCase, caseName, 'outputEnabled');


%% Combine the data from parsers

% Combine the nodes of both meshes together
allNodes = [[fldMsh.nodeIDs;strMsh.nodeIDs],[fldMsh.nodes;strMsh.nodes]];
[globalMsh.nodes,IA] = unique(allNodes,'rows');
globalMsh.nodes = globalMsh.nodes(:,2:4);

% Combine the elements of both meshes together
allElements = [[fldMsh.elementIDs;strMsh.elementIDs],[fldMsh.elements;strMsh.elements(:,1:3)]];
allElements = sort(allElements,1);
globalMsh.elements = allElements(:,2:4);

% Find the shared nodes between both meshes
allNodes = allNodes(:,1);
allNodes(IA) = [];
sharedNodeIDs = allNodes;

% Check if the coupled boundary nodes are equal for testing purposes
checkEqual = isequal(propStrFSI.coupledNodeIDs, propFldFSI.coupledNodeIDs, sharedNodeIDs);

%% Plot coupled boundary nodes
figure(1);
hold on
title('Coupled Boudary Nodes');
axis equal

% Plot all mesh nodes
plot(globalMsh.nodes(:,1),globalMsh.nodes(:,2),'gx');

% Plot all coupled structure nodes
coupledNodesStr = globalMsh.nodes(propStrFSI.coupledNodeIDs,:);
plot(coupledNodesStr(:,1),coupledNodesStr(:,2),'ro');

% Plot all coupled fluid nodes
coupledNodesFld = globalMsh.nodes(propFldFSI.coupledNodeIDs,:);
plot(coupledNodesFld(:,1),coupledNodesFld(:,2),'b+');

% Plot coupled nodes that were automaticaly calculated
sharedNodes = globalMsh.nodes(sharedNodeIDs,:);
plot(sharedNodes(:,1),sharedNodes(:,2),'bd');
hold off

%% Plot all structure and fluid nodes
figure(2);
hold on
title('Fluid and Structure Nodes');
axis equal

% Plot Structure nodes
strNodes = globalMsh.nodes(strMsh.nodeIDs,:);
plot(strNodes(:,1),strNodes(:,2),'ro');

% Plot Fluid nodes
fldNodes = globalMsh.nodes(fldMsh.nodeIDs,:);
plot(fldNodes(:,1),fldNodes(:,2),'gx');
hold off

%% END OF SCRIPT
