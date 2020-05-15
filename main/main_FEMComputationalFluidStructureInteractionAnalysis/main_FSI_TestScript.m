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


%% Create a global mesh from both parsers

% Combine the nodes of both meshes together
allNodes = [fldMsh.nodes;strMsh.nodes];
[globalMsh.nodes,uniqueNodeIDs] = unique(allNodes,'rows');

% Combine the elements of both meshes together
allElements = [fldMsh.elements;strMsh.elements(:,1:4)];
globalMsh.elements = sort(allElements,1);

% Find the shared nodes between both meshes
allNodes = allNodes(:,1);
allNodes(uniqueNodeIDs) = [];
globalMsh.coupledNodeIDs = allNodes;

% Check if the coupled boundary nodes are equal for testing purposes
checkEqual = isequal(propStrFSI.coupledNodeIDs, propFldFSI.coupledNodeIDs, globalMsh.coupledNodeIDs);

%% Convert to local numbering for structure meshe
for iEl = 1:size(strMsh.elements,1)
    strMsh.elements(iEl,1) = iEl; 
    strMsh.elements(iEl,2) = find(strMsh.nodes(:,1) == strMsh.elements(iEl,2));
    strMsh.elements(iEl,3) = find(strMsh.nodes(:,1) == strMsh.elements(iEl,3));
    strMsh.elements(iEl,4) = find(strMsh.nodes(:,1) == strMsh.elements(iEl,4));   
end

% Change the node numbering of FSI coupled boundary
for iNode = 1:propStrFSI.numCoupledNodes
    propStrFSI.coupledNodeIDs(iNode) = find(strMsh.nodes(:,1) == propStrFSI.coupledNodeIDs(iNode));
end

% Change the nodes numbering of structure mesh
strMsh.nodes(:,1) = (1:size(strMsh.nodes,1))';

%% Convert to local numbering for fluid mesh
for iEl = 1:size(fldMsh.elements,1)
    fldMsh.elements(iEl,1) = iEl; 
    fldMsh.elements(iEl,2) = find(fldMsh.nodes(:,1) == fldMsh.elements(iEl,2));
    fldMsh.elements(iEl,3) = find(fldMsh.nodes(:,1) == fldMsh.elements(iEl,3));
    fldMsh.elements(iEl,4) = find(fldMsh.nodes(:,1) == fldMsh.elements(iEl,4));   
end

% Change the node numbering of FSI coupled boundary
for iNode = 1:propStrFSI.numCoupledNodes
    propFldFSI.coupledNodeIDs(iNode) = find(fldMsh.nodes(:,1) == propFldFSI.coupledNodeIDs(iNode));
end

% Change the nodes numbering of structure mesh
fldMsh.nodes(:,1) = (1:size(fldMsh.nodes,1))';


%% Plot coupled boundary nodes
figure(1);
hold on
title('Coupled Boudary Nodes');
axis equal

% Plot all mesh nodes
plot(globalMsh.nodes(:,2),globalMsh.nodes(:,3),'gx');

% Plot all coupled structure nodes
coupledNodesStr = strMsh.nodes(propStrFSI.coupledNodeIDs,:);
plot(coupledNodesStr(:,2),coupledNodesStr(:,3),'ro');

% Plot all coupled fluid nodes
coupledNodesFld = fldMsh.nodes(propFldFSI.coupledNodeIDs,:);
plot(coupledNodesFld(:,2),coupledNodesFld(:,3),'b+');

% Plot coupled nodes that were automaticaly calculated
sharedNodes = globalMsh.nodes(globalMsh.coupledNodeIDs,:);
plot(sharedNodes(:,2),sharedNodes(:,3),'bd');
hold off

%% Plot all structure and fluid nodes
figure(2);
hold on
title('Fluid and Structure Nodes');
axis equal

% Plot Structure nodes
plot(strMsh.nodes(:,2),strMsh.nodes(:,3),'ro');

% Plot Fluid nodes
plot(fldMsh.nodes(:,2),fldMsh.nodes(:,3),'gx');
hold off

%% Plot all structure and fluid mesh
figure(3);
hold on
title('Fluid and Structure Meshes');
axis equal

% Plot Structure mesh
trimesh(strMsh.elements(:,2:end), strMsh.nodes(:,2), strMsh.nodes(:,3), strMsh.nodes(:,4), ...
    'edgecolor', 'black', 'facecolor', 'red');

% Plot Fluid mesh
trimesh(fldMsh.elements(:,2:end), fldMsh.nodes(:,2), fldMsh.nodes(:,3), fldMsh.nodes(:,4), ...
    'edgecolor', 'black', 'facecolor', 'green');

hold off

%% END OF SCRIPT
