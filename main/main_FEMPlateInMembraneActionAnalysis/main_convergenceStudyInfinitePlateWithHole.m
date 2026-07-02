%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Convergence study for the plane stress analysis over the a plate
%        with a hole subject to a constant force at -infty (X-axis)
%
% Date : 02.01.2016
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

%% GUI

% Pick up the coordinates of the node for which to plot the displacement
% field over the refinements steps
% Coordinates of the node on which to compute the displacement field
x = -4;
y = 0;
z = 0;
nodeCoord = [x y z];

% Function handle to the body force vector computation
computeBodyForces = @computeConstantVerticalStructureBodyForceVct;

% Function handle to the linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Output properties
propOutput.isOutput = false;
propOutput.writeOutputToFile = 'undefined';
propOutput.VTKResultFile = 'undefined';

% Function handle to the computation of the linear stiffness matrix
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST;
% computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed;

% Geometry characteristics for the case
internalRadius = 4;
externalRadius = 5;

% Quadrature for the integration of the stiffness matrix and the load
% vector
intLoad.type = 'default';
intDomain.type = 'default';
intLoad.noGP = 1;
intDomain.noGP = 1;

% Quadrature for the computation of the error
intError.type = 'user';
intError.noGP = 8;

% General problem parameters
radiusHole = 1;
forceAmplitude = 10;

% Define the properties for the error computation
propError.resultant = 'stress';
propError.component = '2norm';

% Initialize graphics index
graph.index = 1;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/refinementStudyInfinitePlateWithHole/';

%% Compute an overkill solution

% Name of the case
% caseName = 'infinitePlateWithHoleCoarse';
% caseName = 'infinitePlateWithHole';
% caseName = 'infinitePlateWithHoleFine';
% caseName = 'infinitePlateWithHoleVeryFine';
caseName = 'infinitePlateWithHoleQuadrilaterals';

% Parse the case
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, ~, ~] = ...
    parse_StructuralModelFromGid...
    (pathToCase, caseName, 'outputEnabled');

% Find the global numbering of the selected for postprocessing node
for iNodes = 1:length(strMsh.nodes)
    if strMsh.nodes(iNodes, 2) == nodeCoord(1, 1) && ...
        strMsh.nodes(iNodes, 3) == nodeCoord(1, 2) && ...
        strMsh.nodes(iNodes, 4) == nodeCoord(1, 3)
        nodeID = iNodes;
    end
end
if ~exist('nodeID','var')
    error('The node over which to compute the displacement field was not found');
end

% plot the reference configuration
F = computeLoadVctFEMPlateInMembraneAction...
    (strMsh, propAnalysis, propNBC, 0, intDomain, 'outputEnabled');
plot_referenceConfigurationFEMPlateInMembraneAction...
    (strMsh, propAnalysis, F, homDOFs, [], graph, 'outputEnabled');

% Initialize solution
numNodes = length(strMsh.nodes(:,1));
numDOFs = 2*numNodes;
dHat = zeros(numDOFs,1);

% Solve for the discrete displacement field of the overkill solution
[dHat, FComplete, minElEdgeSizeOverkill] = ...
    solve_FEMPlateInMembraneAction...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    propNBC, computeBodyForces, parameters, computeStiffMtxLoadVct, ...
    solve_LinearSystem, propNLinearAnalysis, intDomain, propOutput, ...
    caseName, pathToOutput, 'outputEnabled');

% Compute the displacement field of the selected for postprocessing node
displacementOverkill = sqrt(dHat(2*nodeID - 1)^2 + dHat(2*nodeID)^2);

%% Perform a convergence study

% Number of refinement steps
noRef = 6;

% Define the meshes corresponding to the refinement
caseNames = {'infinitePlateWithHole_El36' 'infinitePlateWithHole_El52' ...
    'infinitePlateWithHole_El97' 'infinitePlateWithHole_El252' ...
    'infinitePlateWithHole_El586' 'infinitePlateWithHole_El1026'};

% Number of elements for each refinement
noElemnts = [36; 52; 97; 252; 586; 1026];

% Initialize arrays related to the graphs for the convergence study
relErrorStress = zeros(noRef, 1);
relErrorDisplacement = zeros(noRef, 1);
minElEdgeSize = zeros(noRef, 1);
displacement = zeros(noRef, 1);

% Loop over all the refinement steps
for iRefStep = 1:noRef
    % Get the corresponding case name
    caseNameCurrent = strcat('refinementStudyInfinitePlateWithHole/', caseNames{iRefStep});
    
    % Parse the corresponding case
    [strMsh,homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
        parameters, propNLinearAnalysis, ~, ~, ~] = ...
        parse_StructuralModelFromGid ...
        (pathToCase, caseNameCurrent, 'outputEnabled');
    
    % Find the global node numbering of the slected for postprocessing node
    for iNodes = 1:length(strMsh.nodes)
        if strMsh.nodes(iNodes, 2) == nodeCoord(1, 1) && ...
            strMsh.nodes(iNodes, 3) == nodeCoord(1, 2) && ...
            strMsh.nodes(iNodes, 4) == nodeCoord(1, 3)
            nodeID = iNodes;
        end
    end
    if ~exist('nodeID','var')
        error('The node over which to compute the displacement field was not found');
    end
    
    % Initialize solution
    numNodes = length(strMsh.nodes(:,1));
    numDOFs = 2*numNodes;
    dHat = zeros(numDOFs,1);
    
    % Solve the plane stress problem for the current refinement step
    [dHat, FComplete, minElEdgeSize(iRefStep,1)] = ...
        solve_FEMPlateInMembraneAction ...
        (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
        propNBC, computeBodyForces, parameters, computeStiffMtxLoadVct, ...
        solve_LinearSystem, propNLinearAnalysis, intDomain, propOutput, ...
        caseName, pathToOutput, 'outputEnabled');
    
    % Compute the displacement field for the slected for postprocessing node
    displacement(iRefStep, 1) = sqrt(dHat(2*nodeID - 1)^2 + dHat(2*nodeID)^2);
    
    % Compute the relative error for the selected for postprocessing node
    % in the displacement field
    relErrorDisplacement(iRefStep, 1) = ...
        norm(displacement(iRefStep, 1) - displacementOverkill)/ ...
        norm(displacementOverkill);
    
    % Compute the error in the L2-norm over the domain for the selected
    % resultant component
    relErrorStress(iRefStep, 1) = ...
        computeRelErrorL2InfinitePlateWithHoleFEMPlateInMembraneAction...
        (strMsh, dHat, parameters, radiusHole, forceAmplitude, propError, ...
        intError, 'outputEnabled');
end

%% Postprocessing

% Plot the von Mises stress field on the deformed configuration for the
% last refinement mesh.
figure(graph.index);
clf;
hold on;
coords = strMsh.nodes(:, 2:3);
elements = strMsh.elements(:, 2:end);
dHatMatrix = reshape(dHat, 2, []).';
nodalDisplacementMagnitude = sqrt(sum(dHatMatrix.^2, 2));
maxDisplacement = max(nodalDisplacementMagnitude);
domainSize = max(max(coords) - min(coords));
if maxDisplacement > 0
    deformationScale = 0.10*domainSize/maxDisplacement;
else
    deformationScale = 1.0;
end
coordsDeformed = coords + deformationScale*dHatMatrix;

[~, sigma] = computePostprocFEMPlateInMembraneActionCSTLinear ...
    (strMsh, propAnalysis, parameters, dHat);
vonMisesElement = sqrt(sigma(1, :).^2 - sigma(1, :).*sigma(2, :) + ...
    sigma(2, :).^2 + 3*sigma(3, :).^2);
vonMisesNodal = zeros(size(coords, 1), 1);
nodalCounter = zeros(size(coords, 1), 1);
for iElement = 1:size(elements, 1)
    elementNodes = elements(iElement, :);
    elementNodes(isnan(elementNodes)) = [];
    vonMisesNodal(elementNodes) = vonMisesNodal(elementNodes) + ...
        vonMisesElement(iElement);
    nodalCounter(elementNodes) = nodalCounter(elementNodes) + 1;
end
vonMisesNodal = vonMisesNodal./max(nodalCounter, 1);

patch('Faces', elements, 'Vertices', coordsDeformed, ...
    'FaceVertexCData', vonMisesNodal, 'FaceColor', 'interp', ...
    'EdgeColor', [0.05 0.05 0.05], 'LineWidth', 0.20);
axis equal;
axis tight;
grid on;
box on;
colormap(parula);
colorbar;
title(sprintf('von Mises stress on deformed configuration, scale %.1f', ...
    deformationScale));
xlabel('x');
ylabel('y');
hold off;
graph.index = graph.index + 1;

%% Convergence graphs

% Plot the relative error of the stresses in the L2-norm against the
% minimum element edge size
figure(graph.index)
loglog(minElEdgeSize, relErrorStress, '-ob', 'LineWidth', 2, ...
    'MarkerSize', 7);
grid on;
xlabel('Minimum element edge size');
ylabel('$\| \sigma - \sigma_h \|_{L^2}/\|\sigma\|_{L^2}$', ...
    'Interpreter', 'latex');
graph.index = graph.index + 1;

% Plot the relative error of the displacement of the selected node against
% the minimum element edge size
figure(graph.index)
loglog(minElEdgeSize, relErrorDisplacement, '-ob', 'LineWidth', 2, ...
    'MarkerSize', 7);
grid on;
xlabel('Minimum element edge size');
ylabel('$|u-u_h|/|u|$', 'Interpreter', 'latex');
graph.index = graph.index + 1;

% Plot the displacememt of the selected node against the number of elements
figure(graph.index)
plot(noElemnts, displacement, '-ob', 'LineWidth', 2, 'MarkerSize', 7);
grid on;
xlabel('No. elements');
ylabel('$\|u_h\|_2$', 'Interpreter', 'latex');
graph.index = graph.index + 1;

%% END OF THE SCRIPT
