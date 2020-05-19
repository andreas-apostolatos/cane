function writeOutputFEMPlateInMembraneActionToVTK ...
    (analysis, propNLinearAnalysis, propTransientAnalysis, strMsh, ...
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
%                  Input :
%               analysis : Information on the analysis
%                           .type: Analysis type
%    propNLinearAnalysis : Properties of the nonlinear method
%                         .method : 'NEWTON_RAPHSON', 'UNDEFINED', etc
% propTransientAnalysis : Properties of the transient analysis (dummy 
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
% 4. Close files
%
%% Function main body

%% 0. Read input

% Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput,caseName),'dir');
if ~isExistent
    mkdir(strcat(pathToOutput,caseName));
end

%  Number of nodes in the mesh
[noNodes,~] = size(strMsh.nodes);

% Number of elements in the mesh
[noElements,elementOrder] = size(strMsh.elements(:,2:end));

output = fopen(strcat(pathToOutput,caseName,'/',caseName,'_',...
    num2str(noTimeStep),'.vtk'),'w');

% Transpose the nodal coordinates array
XYZ = strMsh.nodes(1:noNodes,2:end)';


% Get the indices of the nodes in the element list with respect to their
% ordering in the nodes array (necessary step when the node IDs are not 
% sequentially numbered starting from 1)
[~, idxElements] = ...
    ismember(strMsh.elements(:, 2:elementOrder + 1), strMsh.nodes(:, 1));

% Re-arrange the element numbering to start from zero
elements = zeros(elementOrder,noElements);
% elements(1:elementOrder,1:noElements) = ...
%     idxElements(1:noElements,1:elementOrder)' - 1;
elements(1:elementOrder,1:noElements) = strMsh.elements(1:noElements,2:elementOrder+1)' - 1;

%% 1. Compute the strain [3,[epsilonXX epsilonYY epsilonXY]] and stress vectors [3,[sigmaXX sigmaYY sigmaXY]]
if strcmp(propNLinearAnalysis.method, 'UNDEFINED')
    [epsilon,sigma] = computePostprocFEMPlateInMembraneActionCSTLinear...
        (strMsh,analysis,parameters,d);
elseif strcmp(propNLinearAnalysis.method, 'NEWTON_RAPHSON')
    [epsilon,sigma] = computePostprocFEMPlateInMembraneActionCSTNLinear...
        (strMsh,analysis,parameters,d);
    warning('Postpsocessing tensors are computed with the linear strain measures');
else
    error('Selected nonlinear method unknown');
end

%% 2. Re-arrange the solution vector into a 2D array as [[x-comp y-comp z-comp],noNodes]
nodalDisplacement = [d(DOF4Output(1,:))'; d(DOF4Output(2,:))'; zeros(1,noNodes)];

%% 3. Write out the data for the color plots

% % Print out analysis information
% if strcmp(propNLinearAnalysis.method,'UNDEFINED')
%     if strcmp(analysis.type,'PLANE_STRESS')
%         fprintf(output, 'Geometrically linear plane stress analysis\n');
%     elseif strcmp(analysis.type,'PLANE_STRAIN')
%         fprintf(output, 'Geometrically linear plane strain analysis\n');
%     else
%         fprintf(output, 'Error in the analysis.type variable\n');
%         error('Either PLANE_STRESS or PLANE_STRAIN must be selected as analysis type');
%     end
% elseif strcmp(propNLinearAnalysis.method,'NEWTON_RAPHSON')
%     if strcmp(analysis.type,'PLANE_STRESS')
%         fprintf(output, 'Geometrically linear plane stress analysis\n');
%     elseif strcmp(analysis.type,'PLANE_STRAIN')
%         fprintf(output, 'Geometrically linear plane strain analysis\n');
%     else
%         fprintf(output, 'Error in the analysis.type variable\n');
%         error('Either PLANE_STRESS or PLANE_STRAIN must be selected as analysis type');
%     end
% end
% fprintf(output, '\n');

% Write out the preamble
fprintf(output, '# vtk DataFile Version 2.0\n');
fprintf(output, '%s\n',title);
fprintf(output, 'ASCII\n');
fprintf(output, '\n');
fprintf(output, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(output, 'POINTS %d double\n', noNodes);

% Write out the nodal coordinates
for nodeID = 1:noNodes
    fprintf(output,'  %f  %f  %f\n', XYZ(:,nodeID) + nodalDisplacement(:,nodeID));
end

% Note that CELL_SIZE uses ELEMENT_ORDER + 1 because the order of each 
% element is included as a data item.
cellSize = 0;
for iElmtId = 1:noElements
    check = isnan(elements(:,iElmtId));
    if check(4) == 1
    cellSize = cellSize + 4;
    else
    cellSize = cellSize + 5;
    end
end

% Output the element connectivities to the nodes
fprintf(output,'\n' );
fprintf(output,'CELLS  %d  %d\n',noElements,cellSize);

% Loop over all the elements in the mesh
for iElmtId = 1:noElements
    check = isnan(elements(:,iElmtId));
    if check(4) == 1
        elementOrder = 3;
        % Write out the element order
        fprintf (output,'  %d',elementOrder);
        
        % Loop over all the polynomial orders
        for order = 1:elementOrder
            % Write out the polynomial order
            fprintf(output,'  %d',elements(order,iElmtId));
        end
        
        % Change line
        fprintf (output,'\n' );
    else
        % Write out the element order
        fprintf (output,'  %d',elementOrder);
        
        % Loop over all the polynomial orders
        for order = 1:elementOrder
            % Write out the polynomial order
            fprintf(output,'  %d',elements(order,iElmtId));
        end
        
        % Change line
        fprintf (output,'\n' );
    end
end

% VTK has a cell type 22 for quadratic triangles.  However, we
% are going to strip the data down to linear triangles for now,
% which is cell type 5.

fprintf (output,'\n' );
fprintf (output,'CELL_TYPES %d\n',noElements);

% Loop over all the elements and write out the nodal coordinates and the
% element connectivities according to the element order
for iElmtId = 1:noElements
    check = isnan(elements(:,iElmtId));
    if check(4) == 1
        fprintf (output,'5\n');
    else
        fprintf (output,'9\n');
    end
end

% Write out the strain tensor
fprintf(output,'\n' );
fprintf(output,'CELL_DATA %d\n',noElements);
fprintf(output,'TENSORS strain double\n');
for iElmtId = 1 : noElements
    fprintf(output,'  %f  %f  %f\n',epsilon(1,iElmtId),epsilon(3,iElmtId),0.0);
    fprintf(output,'  %f  %f  %f\n',epsilon(3,iElmtId),epsilon(2,iElmtId),0.0);
    fprintf(output,'  %f  %f  %f\n',0.0,0.0,0.0);
    if iElmtId ~= noElements
        fprintf(output,'\n');
    end
end

% Write out the displacement vector field at each node
fprintf(output,'POINT_DATA %d\n',noNodes);
fprintf(output,'VECTORS displacements double\n');
for nodeID = 1:noNodes
    fprintf(output,'  %f  %f  %f\n',nodalDisplacement(:,nodeID));
end

% Write out the stress tensor
fprintf(output,'\n' );
fprintf(output,'CELL_DATA %d\n',noElements);
fprintf(output,'TENSORS stress double\n');
for iElmtId = 1 : noElements
    fprintf(output,'  %f  %f  %f\n',sigma(1,iElmtId),sigma(3,iElmtId),0.0);
    fprintf(output,'  %f  %f  %f\n',sigma(3,iElmtId),sigma(2,iElmtId),0.0);
    fprintf(output,'  %f  %f  %f\n\n',0.0,0.0,0.0);
    if iElmtId~=noElements
        fprintf(output,'\n');
    end
end

%% 4. Close files
fclose(output);

return;

end
