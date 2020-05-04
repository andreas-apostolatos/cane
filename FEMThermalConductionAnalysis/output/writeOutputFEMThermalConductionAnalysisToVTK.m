function writeOutputFEMThermalConductionAnalysisToVTK ...
    (propAnalysis, propNLinearAnalysis, propTransientAnalysis, strMsh, ...
    parameters, d, dDot, dDDot, DOF4Output, caseName, pathToOutput, ...
    title, noTimeStep)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documentation 
%
% Writes out the results of a 2D thermal conduction transfer analysis into 
% a VTK file to be read by paraview.
%
%                  Input :
%           propAnalysis : Structure containing general information on the 
%                          analysis,
%                               .type: Analysis type
%    propNLinearAnalysis : Properties of the nonlinear method
%                         .method : 'NEWTON_RAPHSON', 'UNDEFINED', etc
% propTransientAnalysis : Properties of the transient analysis (dummy 
%                         variable for this function)
%                strMsh : Nodes and elements in the mesh
%            parameters : The technical parameters of the problem
%                     d : The solution (temperature) field
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
% 1. Re-arrange the solution vector into a 2D array as [[x-comp y-comp z-comp],noNodes]
%
% 2. Write out the data for the color plots
%
% 3. Close files
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
[noElements,elementOrder] = size(strMsh.elements);

output = fopen(strcat(pathToOutput,caseName,'/',caseName,'_',...
    num2str(noTimeStep),'.vtk'),'w');

% Transpose the nodal coordinates array
XYZ = strMsh.nodes(1:noNodes,:)';

% Re-arrange the element numbering to start from zero
elements = zeros(elementOrder,noElements);
elements(1:elementOrder,1:noElements) = strMsh.elements(1:noElements,1:elementOrder)' - 1;

%% 1. Re-arrange the solution vector into a 2D array as [[x-comp y-comp z-comp],noNodes]
nodalTemperature = [d(DOF4Output(1,:))'; zeros(1, noNodes); zeros(1,noNodes)];

%% 2. Write out the data for the color plots

% Write out the preamble
fprintf(output, '# vtk DataFile Version 2.0\n');
fprintf(output, '%s\n',title);
fprintf(output, 'ASCII\n');
fprintf(output, '\n');
fprintf(output, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(output, 'POINTS %d double\n', noNodes);

% Write out the nodal coordinates
for nodeID = 1:noNodes
    fprintf(output,'  %f  %f  %f\n', XYZ(:,nodeID));
end

% Note that CELL_SIZE uses ELEMENT_ORDER + 1 because the order of each 
% element is included as a data item.
cellSize = 0;
for elementID = 1:noElements
    check = isnan(elements(:,elementID));
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
for elementID = 1:noElements
    check = isnan(elements(:,elementID));
    if check(4) == 1
        elementOrder = 3;
        % Write out the element order
        fprintf (output,'  %d',elementOrder);
        
        % Loop over all the polynomial orders
        for order = 1:elementOrder
            % Write out the polynomial order
            fprintf(output,'  %d',elements(order,elementID));
        end
        
        % Change line
        fprintf (output,'\n' );
    else
        % Write out the element order
        fprintf (output,'  %d',elementOrder);
        
        % Loop over all the polynomial orders
        for order = 1:elementOrder
            % Write out the polynomial order
            fprintf(output,'  %d',elements(order,elementID));
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
for elementID = 1:noElements
    check = isnan(elements(:,elementID));
    if check(4) == 1
        fprintf (output,'5\n');
    else
        fprintf (output,'9\n');
    end
end

% Write out the temperature field at each node
fprintf(output,'POINT_DATA %d\n',noNodes);
fprintf(output,'SCALARS Temperature double\n');
fprintf(output, 'LOOKUP_TABLE default\n');
for nodeID = 1:noNodes
    fprintf(output,'  %f\n', nodalTemperature(1,nodeID));
end

%% 3. Close files
fclose(output);

end
