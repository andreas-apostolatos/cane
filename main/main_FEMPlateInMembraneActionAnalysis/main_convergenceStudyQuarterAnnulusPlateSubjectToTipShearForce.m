%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Convergence study for the plane stress analysis over the quarter
%        of annulus plate subject to tip shear force
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
x = 5;
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

% Define the properties for the error computation
propError.resultant = 'stress';
propError.component = 'tensor';

% Initialize graphics index
graph.index = 1;

%% Compute an overkill solution

% Name of the case
caseName = strcat('refinementStudyCurvedBeamTipShear/', 'curvedBeamTipShear_overkill');

% Parse the case
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, ~, ~, ~] = ...
    parse_StructuralModelFromGid...
    (pathToCase, caseName, 'outputEnabled');

% Find the global numbering of the selected for postprocessing node
for counterNodes = 1:length(strMsh.nodes)
    if strMsh.nodes(counterNodes,1) == nodeCoord(1,1) && ...
        strMsh.nodes(counterNodes,2) == nodeCoord(1,2) && ...
        strMsh.nodes(counterNodes,3) == nodeCoord(1,3)
        nodeID = counterNodes;
    end
end
if ~exist('nodeID','var')
    error('The node over which to compute the displacement field was not found');
end

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
    caseName, 'undefined', 'outputEnabled');

% Compute the displacement field of the selected for postprocessing node
displacementOverkill = sqrt(dHat(2*nodeID - 1)^2 + dHat(2*nodeID)^2);

%% Perform a convergence study

% Number of refinement steps
noRef = 10;

% Define the meshes corresponding to the refinement
caseNames = {'curvedBeamTipShear_El2' 'curvedBeamTipShear_El4' ...
    'curvedBeamTipShear_El8' 'curvedBeamTipShear_El16' ...
    'curvedBeamTipShear_El62' 'curvedBeamTipShear_El374' ...
    'curvedBeamTipShear_El1656' 'curvedBeamTipShear_El6730' ...
    'curvedBeamTipShear_El25768' 'curvedBeamTipShear_El162164'};

% Number of elements for each refinement
noElemnts = [2; 4; 8; 16; 62; 374; 1656; 16730; 25768; 162164];

% Initialize arrays related to the graphs for the convergence study
relErrorStress = zeros(noRef, 1);
relErrorDisplacement = zeros(noRef, 1);
minElEdgeSize = zeros(noRef, 1);
displacement = zeros(noRef, 1);

% Loop over all the refinement steps
for counterRefStep = 1:noRef
    % Get the corresponding case name
    caseName = strcat('refinementStudyCurvedBeamTipShear/', caseNames{counterRefStep});
    
    % Parse the corresponding case
    [strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
        parameters, propNLinearAnalysis, propStrDynamics] = ...
        parse_StructuralModelFromGid...
        (pathToCase, caseName, 'outputEnabled');
    
    % Find the global node numbering of the slected for postprocessing node
    for counterNodes = 1:length(strMsh.nodes)
        if strMsh.nodes(counterNodes, 1) == nodeCoord(1, 1) && ...
            strMsh.nodes(counterNodes, 2) == nodeCoord(1, 2) && ...
            strMsh.nodes(counterNodes, 3) == nodeCoord(1, 3)
            nodeID = counterNodes;
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
    [dHat, FComplete, minElEdgeSize(counterRefStep,1)] = ...
        solve_FEMPlateInMembraneAction...
        (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
        propNBC, computeBodyForces, parameters, computeStiffMtxLoadVct, ...
        solve_LinearSystem, propNLinearAnalysis, intDomain, propOutput, ...
        caseName, 'undefined', 'outputEnabled');
    
    % Compute the displacement field for the slected for postprocessing node
    displacement(counterRefStep, 1) = sqrt(dHat(2*nodeID - 1)^2 + dHat(2*nodeID)^2);
    
    % Compute the relative error for the selected for postprocessing node
    % in the displacement field
    relErrorDisplacement(counterRefStep, 1) = ...
        norm(displacement(counterRefStep, 1) - displacementOverkill)/...
        norm(displacementOverkill);
    
    % Get a node on the Neumann boundary
    nodeNeumann = strMsh.nodes(propNBC.nodes(1, 1), :);
    
    % Get the corresponding function handle for the computation of the load
    funHandle = str2func(propNBC.fctHandle(1, :));
    
    % Compute the force amplitude for the selected node on the Neumann
    % boundary
    forceAmplitude = norm(funHandle(nodeNeumann(1, 1), nodeNeumann(1, 2), ...
        nodeNeumann(1, 3), 0, propNBC));
    
    % Compute the error in the L2-norm over the domain for the selected
    % resultant component
    relErrorStress(counterRefStep, 1) = ...
        computeRelErrorL2CurvedBeamTipShearFEMPlateInMembraneAction ...
        (strMsh, dHat, parameters, internalRadius, externalRadius, ...
        forceAmplitude, propError, intError, 'outputEnabled');
end

%% Plot the corresponding convergence graphs

% Plot the relative error of the stresses in the L2-norm against the
% minimum element edge size
figure(graph.index)
loglog(minElEdgeSize, relErrorStress);
grid on;
graph.index = graph.index + 1;

% Plot the relative error of the displacement of the selected node against
% the minimum element edge size
warning('There is a bug in the computation of the error in terms of stresses');
figure(graph.index)
loglog(minElEdgeSize, relErrorDisplacement);
grid on;
graph.index = graph.index + 1;

% Plot the displacememt of the selected node against the number of elements
figure(graph.index)
semilogx(noElemnts, displacement);
grid on;
graph.index = graph.index + 1;

%% END OF THE SCRIPT
