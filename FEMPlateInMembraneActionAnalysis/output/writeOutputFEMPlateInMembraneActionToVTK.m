function writeOutputFEMPlateInMembraneActionToVTK ...
    (propAnalysis, propNLinearAnalysis, propStrDynamics, strMsh, ...
    parameters, d, dDot, dDDot, DOF4Output, caseName, pathToOutput, ...
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
% Writes out the results of a classical FE plate in membrane action
% analysis into a VTK file to be read by paraview.
%
%                 Input :
%          propAnalysis : Structure containing general information on the 
%                         analysis,
%                          .type: Analysis type
%   propNLinearAnalysis : Properties of the nonlinear method
%                        .method : 'NEWTON_RAPHSON', 'UNDEFINED', etc
%       propStrDynamics : Properties of the transient analysis (dummy 
%                         variable for this function)
%                strMsh : Nodes and elements in the mesh
%            parameters : The technical parameters of the problem
%                     d : The displacement field arranged into a 2D array 
%                         [[x-comp y-comp],noNodes]
%            dDot,dDDot : Dummy variables
%            DOF4Output : Arrangment of the DOFs for writing out the 
%                         results
%              caseName : The name of the case in the inputGiD case folder
%          pathToOutput : Path to the output file
%                 title : The name of the output file
%            noTimeStep : Number of the current time step
%
%                Output :
%                         Write results into file
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the strain [3,[epsilonXX epsilonYY epsilonXY]] and stress vectors [3,[sigmaXX sigmaYY sigmaXY]]
%
% 2. Re-arrange the solution vector into a 2D array as [[x-comp y-comp z-comp],noNodes]
%
% 3. Write out the data for the color plots
%
% 4. Close file
%
% 5. Return if the simulation is steady-state
%
% 6. Write out the data for the rates of the primary fields
%
% 7. Close the file for writting out the results of the rates to the primary fields
%
%% Function main body

%% 0. Read input

% Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput, caseName), 'dir');
if ~isExistent
    mkdir(strcat(pathToOutput, caseName));
end

%  Number of nodes in the mesh
numNodes = length(strMsh.nodes(:, 1));

% Number of DOFs in the computational domain
if strcmp(propAnalysis.type, 'planeStress') || ...
        strcmp(propAnalysis.type, 'planeStrain')
    numDOFs = 2*numNodes;
else
   error('Analysis type %s not supported in this function', ...
       propAnalysis.type);
end

% Number of elements in the mesh
[numElements, elementOrder] = size(strMsh.elements(:, 2:end));

% Open file to write out the results
outputUnitColorPlots = fopen(strcat(pathToOutput, caseName, '/', ...
    caseName, '_contourPlots_', num2str(noTimeStep),'.vtk'), 'w');

% Transpose the nodal coordinates array
XYZ = strMsh.nodes(1:numNodes, 2:end)';

% Get the indices of the nodes in the element list with respect to their
% ordering in the nodes array (necessary step when the node IDs are not 
% sequentially numbered starting from 1)
[~, idxElements] = ...
    ismember(strMsh.elements(:, 2:elementOrder + 1), strMsh.nodes(:, 1));
idxElements(idxElements == 0) = NaN;

% Re-arrange the element numbering to start from zero
elements = idxElements(1:numElements,1:elementOrder)' - 1;

%% 1. Compute the strain [3,[epsilonXX epsilonYY epsilonXY]] and stress vectors [3,[sigmaXX sigmaYY sigmaXY]]
if strcmp(propNLinearAnalysis.method, 'UNDEFINED')
    [epsilon, sigma] = computePostprocFEMPlateInMembraneActionCSTLinear ...
        (strMsh, propAnalysis, parameters, d);
elseif strcmp(propNLinearAnalysis.method, 'NEWTON_RAPHSON')
    [epsilon, sigma] = computePostprocFEMPlateInMembraneActionCSTNLinear ...
        (strMsh, propAnalysis, parameters, d);
else
    error('Selected nonlinear method unknown');
end

%% 2. Re-arrange the solution vector into a 2D array as [[x-comp y-comp z-comp],noNodes]
nodalDisplacement = [d(DOF4Output(1, :))'
                     d(DOF4Output(2, :))'
                     zeros(1, numNodes)];

%% 3. Write out the data for the color plots

% Write out the preamble
fprintf(outputUnitColorPlots, '# vtk DataFile Version 2.0\n');
fprintf(outputUnitColorPlots, '%s\n',title);
fprintf(outputUnitColorPlots, 'ASCII\n');
fprintf(outputUnitColorPlots, '\n');
fprintf(outputUnitColorPlots, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(outputUnitColorPlots, 'POINTS %d double\n', numNodes);

% Write out the nodal coordinates
for nodeID = 1:numNodes
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n', XYZ(:, nodeID) + nodalDisplacement(:, nodeID));
end

% Note that CELL_SIZE uses ELEMENT_ORDER + 1 because the order of each 
% element is included as a data item.
cellSize = 0;
for iElmtId = 1:numElements
    check = isnan(elements(:,iElmtId));
    if check(4) == 1
    cellSize = cellSize + 4;
    else
    cellSize = cellSize + 5;
    end
end

% Output the element connectivities to the nodes
fprintf(outputUnitColorPlots, '\n' );
fprintf(outputUnitColorPlots, 'CELLS  %d  %d\n', numElements, cellSize);

% Loop over all the elements in the mesh
for iElmtId = 1:numElements
    check = isnan(elements(:, iElmtId));
    if check(4) == 1
        elementOrder = 3;
        % Write out the element order
        fprintf (outputUnitColorPlots, '  %d', elementOrder);
        
        % Loop over all the polynomial orders
        for order = 1:elementOrder
            % Write out the polynomial order
            fprintf(outputUnitColorPlots, '  %d', elements(order, iElmtId));
        end
        
        % Change line
        fprintf (outputUnitColorPlots, '\n');
    else
        % Write out the element order
        fprintf (outputUnitColorPlots, '  %d', elementOrder);
        
        % Loop over all the polynomial orders
        for order = 1:elementOrder
            % Write out the polynomial order
            fprintf(outputUnitColorPlots, '  %d', elements(order, iElmtId));
        end
        
        % Change line
        fprintf(outputUnitColorPlots, '\n' );
    end
end

% VTK has a cell type 22 for quadratic triangles.  However, we
% are going to strip the data down to linear triangles for now,
% which is cell type 5.

fprintf (outputUnitColorPlots, '\n' );
fprintf (outputUnitColorPlots, 'CELL_TYPES %d\n', numElements);

% Loop over all the elements and write out the nodal coordinates and the
% element connectivities according to the element order
for iElmtId = 1:numElements
    check = isnan(elements(:, iElmtId));
    if check(4) == 1
        fprintf (outputUnitColorPlots, '5\n');
    else
        fprintf (outputUnitColorPlots, '9\n');
    end
end

% Write out the strain tensor
fprintf(outputUnitColorPlots, '\n' );
fprintf(outputUnitColorPlots, 'CELL_DATA %d\n', numElements);
fprintf(outputUnitColorPlots, 'TENSORS strain double\n');
for iElmtId = 1 : numElements
    fprintf(outputUnitColorPlots, '  %f  %f  %f\n', epsilon(1, iElmtId), epsilon(3, iElmtId), 0.0);
    fprintf(outputUnitColorPlots, '  %f  %f  %f\n', epsilon(3, iElmtId), epsilon(2, iElmtId), 0.0);
    fprintf(outputUnitColorPlots, '  %f  %f  %f\n', 0.0, 0.0, 0.0);
    if iElmtId ~= numElements
        fprintf(outputUnitColorPlots, '\n');
    end
end

% Write out the displacement vector field at each node
fprintf(outputUnitColorPlots,'POINT_DATA %d\n',numNodes);
fprintf(outputUnitColorPlots,'VECTORS displacements double\n');
for nodeID = 1:numNodes
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n',nodalDisplacement(:,nodeID));
end

% Write out the stress tensor
fprintf(outputUnitColorPlots, '\n' );
fprintf(outputUnitColorPlots, 'CELL_DATA %d\n', numElements);
fprintf(outputUnitColorPlots, 'TENSORS stress double\n');
for iElmtId = 1 : numElements
    fprintf(outputUnitColorPlots, '  %f  %f  %f\n', sigma(1, iElmtId), sigma(3, iElmtId), 0.0);
    fprintf(outputUnitColorPlots, '  %f  %f  %f\n', sigma(3, iElmtId), sigma(2, iElmtId), 0.0);
    fprintf(outputUnitColorPlots, '  %f  %f  %f\n\n', 0.0, 0.0, 0.0);
    if iElmtId ~= numElements
        fprintf(outputUnitColorPlots, '\n');
    end
end

%% 4. Close file
fclose(outputUnitColorPlots);

%% 5. Return if the simulation is steady-state
if strcmp(propStrDynamics.timeDependence, 'STEADY_STATE')
    return;
end

%% 6. Write out the data for the rates of the primary fields

% Check whether there is an output unit, if not create one
directoryName = strcat(pathToOutput, caseName, '/dRates');
isExistent = exist(directoryName, 'dir');
if ~isExistent
    mkdir(directoryName);
end

% Open the file to write the restart data
outputUnitDRates = fopen(strcat(directoryName, '/', caseName, '_dRates_', ...
    num2str(noTimeStep), '.txt'),'w');
fprintf(outputUnitDRates, '# Restart file\n');
fprintf(outputUnitDRates, '# The rates of the velocity and pressure fields\n');
fprintf(outputUnitDRates, '\n');
fprintf(outputUnitDRates, 'dDot \n');
str = [num2str(dDot), repmat(' \n', numDOFs, 1)]'; 
fprintf(outputUnitDRates, reshape(str, 1, size(str, 1) * size(str, 2)));
fprintf(outputUnitDRates, '\n');
fprintf(outputUnitDRates, 'dDDot \n');
str = [num2str(dDDot), repmat(' \n', numDOFs, 1)]'; 
fprintf(outputUnitDRates, reshape(str, 1, size(str, 1) * size(str, 2)));

%% 7. Close the file for writting out the results of the rates to the primary fields
fclose(outputUnitDRates);

end