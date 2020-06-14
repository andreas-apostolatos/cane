function writeOutputFEMIncompressibleFlowToVTK ...
    (propAnalysis, propNLinearAnalysis, propFldDynamics, fldMsh, ...
    propParameters, up, upDot, upDDot, DOF4Output, caseName, pathToOutput, ...
    title, noTimeStep)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation 
%
% Writes out the results of the solution to an incompressible flow problem
% using finite elements into a VTK file for the given time step.
%
%               Input :
%        propAnalysis : Structure containing general information on the 
%                       analysis,
%                           .type : Analysis type
% propNLinearAnalysis : On the nonlinear analysis scheme
%     propFldDynamics : Information on transient analysis,
%                           .timeDependence : 'STEADY_STATE' or 'TRANSIENT'
%              fldMsh : Nodes and elements in the mesh
%      propParameters : Structure containing information on the material 
%                       parameters,
%                           .nue : The kinematic viscosity
%                           .rho : The flow density
%                  up : The discrete nodal solution vector
%               upDot : The rate of the nodal solution vector
%              upDDot : The second order rate of the nodal solution vector
%                       (dummy variable for this function)
%          DOF4Output : Array containing the arrangment of the DOFs for
%                       printing them out
%            caseName : The name of the case in the inputGiD case folder
%        pathToOutput : Path to the output file
%      outputFilename : The name of the output file
%               title : The title of the VTK file
%          noTimeStep : The number of the time step
%
%              Output :
%                       Write results into file
%
% Function layout :
%
% 0. Read input
%
% 1. Write out the data for the color plots
%
% 2. Close file
%
% 3. Write out the data for the rates of the primary fields
%
% 4. Close the file for writting out the results of the rates to the primary fields
% 
%% Function main body

%% 0. Read input

% Check the dimensionality of the analysis
if strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
    isAnalysis3D = true;
    noDOFsNode = 4;
elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
    isAnalysis3D = false;
    noDOFsNode = 3;
else
    error('Error in the dimensionality of the analysis');
end

% Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput, caseName), 'dir');
if ~isExistent
    mkdir(strcat(pathToOutput, caseName));
end

%  Number of nodes in the mesh
[numNodes, ~] = size(fldMsh.nodes(:,2:end));

% Number of DOFs
numDOFs = noDOFsNode*numNodes;

% Number of elements in the mesh
[numElements, elementOrder] = size(fldMsh.elements(:,2:end));

% Arrange the velocity and pressure field arranged into a 2D array as per
% [[x-comp y-comp z-comp pressure],noNodes]
if isAnalysis3D
    solutionVct = [up(DOF4Output(1,:)), up(DOF4Output(2,:)), up(DOF4Output(3,:)), up(DOF4Output(4,:))];
else
	solutionVct = [up(DOF4Output(1,:)), up(DOF4Output(2,:)), zeros(numNodes, 1), up(DOF4Output(3,:))];
end

% Open the file to write out the results
outputUnitColorPlots = fopen(strcat(pathToOutput, caseName, '/', ...
    caseName, '_contourPlots_', num2str(noTimeStep), '.vtk'), 'w');

% Transpose the nodal coordinates array
XYZ = fldMsh.nodes(1:numNodes,2:end)';

% Decide according to the element order
if isAnalysis3D
    elOrder = 10;
else
    elOrder = 6;
end
if ( elementOrder == elOrder )
	fprintf(1, '\n' );
    fprintf(1, 'TWO_TO_VTK:\n' );
    fprintf(1, '  The input data uses quadratic elements.\n' );
    fprintf(1, '  The output data will use linear elements.\n' );
end

% Get the indices of the nodes in the element list with respect to their
% ordering in the nodes array (necessary step when the node IDs are not 
% sequentially numbered starting from 1)
[~, idxElements] = ...
    ismember(fldMsh.elements(:, 2:elementOrder + 1), fldMsh.nodes(:, 1));

% Re-arrange the element numbering to start from zero
elements = zeros(elementOrder, numElements);
elements(1:elementOrder, 1:numElements) = ...
    idxElements(1:numElements, 1:elementOrder)' - 1;

%% 1. Write out the data for the color plots

% Write out the preamble
fprintf(outputUnitColorPlots, '# vtk DataFile Version 2.0\n');
fprintf(outputUnitColorPlots, '%s\n',title);
fprintf(outputUnitColorPlots, 'ASCII\n');
fprintf(outputUnitColorPlots, '\n');
fprintf(outputUnitColorPlots, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(outputUnitColorPlots, 'POINTS %d double\n', numNodes);

% Write out the nodal coordinates
str = [num2str(XYZ'), repmat(' \n', numNodes, 1)]'; 
fprintf(outputUnitColorPlots, reshape(str, 1, size(str, 1)*size(str, 2)));

% Note that CELL_SIZE uses ELEMENT_ORDER+1 because the order of each 
% element is included as a data item.
cellSize = numElements*(elementOrder + 1);

% Output the element connectivities to the nodes
fprintf(outputUnitColorPlots,'\n' );
fprintf(outputUnitColorPlots,'CELLS  %d  %d\n', numElements, cellSize);

% Write out the elements in the mesh and their order
str = [num2str([ones(numElements, 1)*elementOrder, elements']), repmat(' \n', numElements, 1)]'; 
fprintf(outputUnitColorPlots, reshape(str, 1, size(str, 1)*size(str, 2)));

% Write out the elements, the nodal coordinates and the element 
% connectivities according to the element order
fprintf (outputUnitColorPlots, '\n');
fprintf (outputUnitColorPlots, 'CELL_TYPES %d\n', numElements);
if elementOrder == noDOFsNode
    if isAnalysis3D
        fprintf(outputUnitColorPlots, repmat('10\n', 1, numElements));
    else
        fprintf(outputUnitColorPlots, repmat('5\n', 1, numElements));
    end
elseif elementOrder == elOrder
    if isAnalysis3D
        fprintf(outputUnitColorPlots, repmat('24\n', 1, numElements));
    else
        fprintf(outputUnitColorPlots, repmat('22\n', 1, numElements));
    end
end

% Write out the pressure field at each npde
fprintf(outputUnitColorPlots, '\n' );
fprintf(outputUnitColorPlots, 'POINT_DATA %d\n',numNodes);
fprintf(outputUnitColorPlots, 'SCALARS pressure double\n');
fprintf (outputUnitColorPlots, 'LOOKUP_TABLE default\n' );
str = [num2str(solutionVct(:, 4)/propParameters.rho), repmat(' \n', numNodes, 1)]'; 
fprintf(outputUnitColorPlots, reshape(str, 1, size(str, 1) * size(str, 2) ) );

% Write out the displacement vector field at each node (the actual pressure 
% field, that is divided by the density of the fluid)
fprintf(outputUnitColorPlots, 'VECTORS velocity double\n');
str = [num2str(solutionVct(:, 1:3)), repmat(' \n', numNodes, 1)]'; 
fprintf(outputUnitColorPlots, reshape(str, 1, size(str, 1)*size(str, 2) ) );

%% 2. Close the file for writting out the results of the primary fields
fclose(outputUnitColorPlots);

%% 3. Return if the simulation is steady-state
if strcmp(propFldDynamics.timeDependence,'STEADY_STATE')
    return;
end

%% 3. Write out the data for the rates of the primary fields

% Check whether there is an output unit, if not create one
directoryName = strcat(pathToOutput, caseName, '/upRates');
isExistent = exist(directoryName, 'dir');
if ~isExistent
    mkdir(directoryName);
end

% Open the file to write the restart data
outputUnitUPRates = fopen(strcat(directoryName, '/', caseName, '_upRates_', ...
    num2str(noTimeStep), '.txt'),'w');
fprintf(outputUnitUPRates, '# Restart file\n');
fprintf(outputUnitUPRates, '# The rates of the velocity and pressure fields\n');
fprintf(outputUnitUPRates, '\n');
fprintf(outputUnitUPRates, 'upRates \n');
str = [num2str(upDot), repmat(' \n', numDOFs, 1)]'; 
fprintf(outputUnitUPRates, reshape(str, 1, size(str, 1) * size(str, 2)));

%% 4. Close the file for writting out the results of the rates to the primary fields
fclose(outputUnitUPRates);

end
