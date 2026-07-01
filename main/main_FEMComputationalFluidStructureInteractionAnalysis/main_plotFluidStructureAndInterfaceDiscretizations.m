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
% Task : Plots the fluid and structural meshes along with the interface
%        nodes shared along the wet surface.
%
% Date : 17.05.2020
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
% caseName = 'testTurek';
caseName = 'turek_fsi';

% Parse fluid model from GiD
[fldMsh, homDOFsFld, inhomDOFsFld, valuesInhomDOFsFld, propALE, ~, ...
    propAnalysisFld, propParametersFld, propNLinearAnalysisFld, ...
    propFldDynamics, propGaussIntFld, propPostProcFld, propFSIFld] = ...
    parse_FluidModelFromGid ...
    (pathToCase, caseName, 'outputEnabled');

% Parse structural model from GiD
[strMsh, homDOFsStr, inhomDOFsStr, valuesInhomDOFsStr, propNBCStr, ...
    propAnalysisStr, propParametersStr, propNLinearAnalysisStr, ...
    propStrDynamics, propGaussIntStr, ~, propFSIStr] = ...
    parse_StructuralModelFromGid...
    (pathToCase, caseName, 'outputEnabled');


%% Create a global mesh from both parsers

% Combine the nodes of both meshes together
allNodes = [fldMsh.nodes
            strMsh.nodes];
[globalMsh.nodes, uniqueNodeIDs] = unique(allNodes, 'rows');

% Combine the elements of both meshes together
allElements = [fldMsh.elements
               strMsh.elements(:, 1:4)];
globalMsh.elements = sort(allElements, 1);

% Find the shared nodes between both meshes
allNodes = allNodes(:, 1);
allNodes(uniqueNodeIDs) = [];
globalMsh.coupledNodeIDs = allNodes;

% Check if the coupled boundary nodes are equal for testing purposes
checkEqual = isequal(propFSIStr.nodes, propFSIFld.nodes, globalMsh.coupledNodeIDs);

%% Convert to local numbering for structure meshe
for iEl = 1:size(strMsh.elements, 1)
    strMsh.elements(iEl, 1) = iEl; 
    strMsh.elements(iEl, 2) = find(strMsh.nodes(:, 1) == ...
        strMsh.elements(iEl, 2));
    strMsh.elements(iEl, 3) = find(strMsh.nodes(:, 1) == ...
        strMsh.elements(iEl, 3));
    strMsh.elements(iEl, 4) = find(strMsh.nodes(:, 1) == ...
        strMsh.elements(iEl, 4));   
end

% Change the node numbering of FSI coupled boundary
for iNode = 1:length(propFSIStr.nodes)
    propFSIStr.nodes(iNode) = ...
        find(strMsh.nodes(:,1) == propFSIStr.nodes(iNode));
end

% Change the nodes numbering of structure mesh
strMsh.nodes(:, 1) = (1:size(strMsh.nodes, 1))';

%% Convert to local numbering for fluid mesh
for iEl = 1:size(fldMsh.elements, 1)
    fldMsh.elements(iEl, 1) = iEl; 
    fldMsh.elements(iEl, 2) = find(fldMsh.nodes(:, 1) == ...
        fldMsh.elements(iEl,2));
    fldMsh.elements(iEl, 3) = find(fldMsh.nodes(:, 1) == ...
        fldMsh.elements(iEl,3));
    fldMsh.elements(iEl, 4) = find(fldMsh.nodes(:, 1) == ...
        fldMsh.elements(iEl,4));   
end

% Change the node numbering of FSI coupled boundary
for iNode = 1:length(propFSIStr.nodes)
    propFSIFld.nodes(iNode) = ...
        find(fldMsh.nodes(:, 1) == propFSIFld.nodes(iNode));
end

% Change the nodes numbering of structure mesh
fldMsh.nodes(:, 1) = (1:length(fldMsh.nodes))';

%% Plot coupled boundary nodes
figure(1);
hold on;
title('Coupled boundary nodes');
axis equal;

% Plot all mesh nodes
plot(globalMsh.nodes(:,2), globalMsh.nodes(:, 3), '.', ...
    'Color', [0.70 0.70 0.70], 'MarkerSize', 3);

% Plot all coupled structure nodes
coupledNodesStr = strMsh.nodes(propFSIStr.nodes, :);
plot(coupledNodesStr(:, 2), coupledNodesStr(:, 3), '^', ...
    'Color', [0.85 0.10 0.10], 'LineWidth', 1.5, 'MarkerSize', 5);

% Plot all coupled fluid nodes
coupledNodesFld = fldMsh.nodes(propFSIFld.nodes, :);
plot(coupledNodesFld(:, 2), coupledNodesFld(:, 3), '+', ...
    'Color', [0.00 0.25 0.80], 'LineWidth', 1.5, 'MarkerSize', 5);

% Plot coupled nodes that were automaticaly calculated
sharedNodes = globalMsh.nodes(globalMsh.coupledNodeIDs, :);
plot(sharedNodes(:, 2), sharedNodes(:, 3), 'ko', ...
    'LineWidth', 1.5, 'MarkerSize', 4);
xlim([0.10 0.70]);
ylim([0.13 0.27]);
grid on;
hold off;

%% Plot all structure and fluid nodes
figure(2);
hold on
title('Fluid and Structure Nodes');
axis equal;

% Plot Structure nodes
plot(strMsh.nodes(:, 2), strMsh.nodes(:, 3), '.', ...
    'Color', [0.85 0.10 0.10], 'MarkerSize', 4);

% Plot Fluid nodes
plot(fldMsh.nodes(:, 2), fldMsh.nodes(:, 3), '.', ...
    'Color', [0.30 0.30 0.30], 'MarkerSize', 3);
xlim([0.10 0.70]);
ylim([0.13 0.27]);
grid on;
hold off;

%% Plot all structure and fluid mesh
figure(3);
hold on;
title('FSI interface mesh near the cylinder and beam');
axis equal;

% Plot Structure mesh
trimesh(strMsh.elements(:, 2:end), strMsh.nodes(:,2), strMsh.nodes(:, 3), ...
    strMsh.nodes(:,4), 'edgecolor', [0.05 0.05 0.05], ...
    'facecolor', [0.82 0.82 0.82], 'LineWidth', 0.25);

% Plot Fluid mesh
trimesh(fldMsh.elements(:,2:end), fldMsh.nodes(:,2), fldMsh.nodes(:,3), ...
    fldMsh.nodes(:,4), 'edgecolor', [0.15 0.15 0.15], ...
    'facecolor', [0.88 0.94 0.97], 'LineWidth', 0.15);

plot(coupledNodesFld(:, 2), coupledNodesFld(:, 3), '.', ...
    'Color', [0.90 0.20 0.05], 'MarkerSize', 8);
xlim([0.13 0.36]);
ylim([0.145 0.255]);
view(2);
grid on;
hold off;

%% Plot Turek FSI benchmark setup schematic
figure(4);
clf;
hold on;
axis equal;
box on;
grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 11);
title('Turek FSI benchmark setup', 'FontName', 'Arial', 'FontSize', 15);

% Channel and immersed structure
rectangle('Position', [0.05 0.05 0.65 0.30], 'EdgeColor', [0.15 0.15 0.15], ...
    'LineWidth', 1.8);
rectangle('Position', [0.15 0.15 0.10 0.10], 'Curvature', [1 1], ...
    'EdgeColor', [0.05 0.05 0.05], 'FaceColor', [0.92 0.92 0.92], ...
    'LineWidth', 1.8);
rectangle('Position', [0.25 0.19 0.35 0.02], 'EdgeColor', [0.05 0.05 0.05], ...
    'FaceColor', [0.70 0.70 0.70], 'LineWidth', 1.5);

% Inflow arrows
for yArrow = [0.13 0.20 0.27]
    quiver(0.075, yArrow, 0.055, 0, 0, 'Color', [0.00 0.25 0.80], ...
        'LineWidth', 2, 'MaxHeadSize', 0.8);
end

% Interface highlight along the beam
plot([0.25 0.60], [0.21 0.21], 'Color', [0.10 0.10 0.10], 'LineWidth', 2);
plot([0.25 0.60], [0.19 0.19], 'Color', [0.10 0.10 0.10], 'LineWidth', 2);

% Labels
text(0.085, 0.305, 'inflow', 'Color', [0.00 0.25 0.80], ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
text(0.165, 0.105, 'cylinder', 'FontName', 'Arial', 'FontSize', 12);
text(0.405, 0.225, 'flexible beam', 'Color', [0.00 0.00 0.00], ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
text(0.61, 0.305, 'outflow', 'FontName', 'Arial', 'FontSize', 12);

xlim([0.04 0.71]);
ylim([0.04 0.36]);
xlabel('x');
ylabel('y');
hold off;

%% END OF SCRIPT
