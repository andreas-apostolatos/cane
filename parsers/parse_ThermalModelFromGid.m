function [strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, ...
    propAnalysis, propParameters, propNLinearAnalysis, propThermalDynamics, ...
    propGaussInt] = parse_ThermalModelFromGid(pathToCase, caseName, outMsg)
%% Licensingprop
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documentation
%
% Parses data from an input file created using GiD for a structural 
% boundary value problem.
%
%               Input :
%          pathToCase : The absolute path to the inputGiD case folder
%            caseName : The name of the case in the inputGiD case folder
%              outMsg : On the output information on the command window
%
%              Output :
%              strMsh : On the structural mesh   
%                           .nodes : The nodes in the FE mesh
%                        .elements : The elements in the FE mesh
%             homDOFs : The global numbering of the nodes where homogeneous
%                       Dirichlet boundary conditions are applied
%           inhomDOFs : The global numbering of the nodes where 
%                       inhomogeneous Dirichlet boundary conditions are 
%                       applied
%     valuesInhomDOFs : The prescribed values for the inhomogeneous 
%                       Dirichlet boundary conditions
%             propNBC : Structure containing information on the Neumann
%                       boundary conditions (fluxes)
%                           .nodes : The nodes where Neumann boundary 
%                                    conditions are applied
%                        .loadType : The type of the load for each Neumann 
%                                    node
%                       .fctHandle : The function handle for each Neumann 
%                                    node for the computation of the load 
%                                    vector (these functions are under the 
%                                    folder load)
%        propAnalysis : Structure defining general properties of the
%                       analysis,
%                             .type : The analysis type
%      propParameters : Problem specific technical (physical) parameters
%                              .rho : Material density
%                                .k : Thermal conductivity
%                               .cp : Specific heat capacity
% propNLinearAnalysis :     .method : The employed nonlinear method
%                        .tolerance : The residual tolerance
%                          .maxIter : The maximum number of the nonlinear 
%                                     iterations
% propThermalDynamics : .timeDependence : On the transient analysis
%                               .method : The time integration method
%                                   .T0 : The start time of the simulation
%                                 .TEnd : The end time of the simulation
%                          .noTimeSteps : The number of the time steps
%        propGaussInt : On the Gauss Point integration
%                               .type : 'default', 'user'
%                         .domainNoGP : Number of Gauss Points for the 
%                                       domain integration
%                       .boundaryNoGP : Number of Gauss Points for the 
%                                       boundary integration
%
%% Function layout :
%
% 1. Load the input file from GiD
%
% 2. Load the analysis type
%
% 3. Load the material properties
%
% 4. Load the nonlinear method (heat transfer without internal heat generation is always linear)
%
% 5. Load the time integration method
%
% 6. Load the Gauss Point integration method
%
% 7. Load the heat nodes
%
% 8. Load the heat elements by connectivity arrays
%
% 9. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
%
% 10. Load the nodes on the Neumann boundary together with the load application information
%
% 11. Get edge connectivity arrays for the Neumann edges
%
% 12. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________\n');
    fprintf('#############################################################\n');
    fprintf('Parsing data from GiD input file for a heat transfer boundary\n');
    fprintf('value problem has been initiated\n');
    fprintf('_____________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Initialize output arrays
homDOFs = [];
inhomDOFs = [];
valuesInhomDOFs = [];

%% 1. Load the input file from GiD
fstring = fileread([pathToCase caseName '.dat']); 

%% 2. Load the analysis type
block = regexp(fstring, 'THERMAL_CONDUCTION_ANALYSIS', 'split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
propAnalysis.type = out{1}{2};
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Analysis type: %s \n', propAnalysis.type);
end

%% 3. Load the material properties
block = regexp(fstring, 'THERMAL_MATERIAL_PROPERTIES','split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ',', 'MultipleDelimsAsOne', 1);

% Read material density
propParameters.rho = str2double(out{1}{2});

% Read thermal conductivity
propParameters.k = str2double(out{1}{4});

% Read specific heat
propParameters.cp = str2double(out{1}{6});

%% 4. Load the nonlinear method (heat transfer without internal heat generation is always linear)
propNLinearAnalysis.method = 'UNDEFINED';
propNLinearAnalysis.noLoadSteps = [];
propNLinearAnalysis.eps = [];
propNLinearAnalysis.maxIter = 1;

%% 5. Load the time integration method
block = regexp(fstring, 'THERMAL_TRANSIENT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
propThermalDynamics.timeDependence = out{1}{2};
if ~strcmp(propThermalDynamics.timeDependence, 'STEADY_STATE')
    propThermalDynamics.method = out{1}{4};
    propThermalDynamics.T0 = str2double(out{1}{6});
    propThermalDynamics.TEnd = str2double(out{1}{8});
    propThermalDynamics.noTimeSteps = str2double(out{1}{10});
    propThermalDynamics.isAdaptive = out{1}{12};
    propThermalDynamics.dt = (propThermalDynamics.TEnd - propThermalDynamics.T0)/propThermalDynamics.noTimeSteps;
end
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Thermal dynamics: %s \n', propThermalDynamics.timeDependence);
    if ~strcmp(propThermalDynamics.timeDependence, 'STEADY_STATE')
        fprintf('\t>> Time integration method: %s \n', propThermalDynamics.method);
        fprintf('\t>> Start time of the simulation: %f \n', propThermalDynamics.T0);
        fprintf('\t>> End time of the simulation: %f \n', propThermalDynamics.TEnd);
        fprintf('\t>> Number of time steps: %f \n', propThermalDynamics.noTimeSteps);
        fprintf('\t>> Time step size: %f \n', propThermalDynamics.dt);
    end
end

%% 6. Load the Gauss Point integration method
block = regexp(fstring, 'THERMAL_INTEGRATION', 'split');
block(1) = [];
out = textscan(block{1}, '%s','delimiter', ' ', 'MultipleDelimsAsOne', 1);
propGaussInt.type = out{1}{2};
if strcmp(propGaussInt.type, 'user')
    propGaussInt.domainNoGP = str2double(out{1}{4});
    propGaussInt.boundaryNoGP = str2double(out{1}{6});
end
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Gauss integration type: %s \n', propGaussInt.type);
    if strcmp(propGaussInt.type, 'user')
        fprintf('\t>> No. Gauss Points for the domain integration: %d \n', propGaussInt.domainNoGP);
        fprintf('\t>> No. Gauss Points for the boundary integration: %d \n', propGaussInt.boundaryNoGP);
    end
end

%% 7. Load the structural nodes
block = regexp(fstring, 'THERMAL_NODES', 'split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %f %f %f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
strMsh.nodes = out(:, 2:4);
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Number of nodes in the mesh: %d \n', length(strMsh.nodes(:, 1)));
end

%% 8. Load the structural elements by connectivity arrays
block = regexp(fstring, 'THERMAL_ELEMENTS', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %f %f %f %f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
noNodes = length(out(1, :));
strMsh.elements = out(:, 2:noNodes);
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Number of elements in the mesh: %d \n', length(strMsh.elements));
end

%% 9. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
block = regexp(fstring, 'THERMAL_DIRICHLET_NODES', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);

% Filter out the actual DOFs
numDOFsNodeGiD = 1;
if strcmp(propAnalysis.type, 'THERMAL_CONDUCTION_2D')
    numDOFsNode = 1;
end

% Initialize counter
counterHomDBC = 1;
counterInhomDBC = 1;

% Find the number of nodes where Dirichlet boundary conditions are applied
numDBCNodes = length(out)/(numDOFsNodeGiD + 1);

% Loop overll all nodes where Dirchlet boundary conditions are applied
for i = 1:numDBCNodes
    % Get the Dirichlet node ID
    nodeID = out((numDOFsNodeGiD + 1)*i - numDOFsNodeGiD);
    
    % Get the x-component of the prescribed value
    for j = 1:numDOFsNode
        presValue = out((numDOFsNodeGiD + 1)*i - numDOFsNodeGiD + j);
        if ~isnan(presValue)
            if presValue == 0
                homDOFs(counterHomDBC) = numDOFsNode*nodeID - numDOFsNode + j;
                counterHomDBC = counterHomDBC + 1;
            else
                inhomDOFs(counterInhomDBC) = numDOFsNode*nodeID - numDOFsNode + j;
                valuesInhomDOFs(counterInhomDBC) = presValue;
                counterInhomDBC = counterInhomDBC + 1;
            end
        end
    end
end

% Sort out the vectors
homDOFs = sort(homDOFs);
[inhomDOFs, indexSorting] = sort(inhomDOFs);
valuesInhomDOFs = valuesInhomDOFs(indexSorting);

%% 10. Load the nodes on the Neumann boundary together with the load application information
block = regexp(fstring,'THERMAL_FLUX_NODES', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %s %s', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
end
if ~isempty(out)
    out = out{1};
    propNBC.nodes = cell2mat(out(:, 1));
    outLoadType = out(:, 2);
    propNBC.loadType = cell2mat(outLoadType{1});
    outFctHandle = out(:, 3);
    propNBC.fctHandle = cell2mat(outFctHandle{1});
end

%% 11. Get edge connectivity arrays for the Neumann edges
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Neumann boundary edges: %d \n', length(propNBC.nodes) - 1);
end

% Initialize the Neumann boundary lines
propNBC.lines = zeros(length(unique(propNBC.nodes)) - 1, 3);

% Initialize line counter
counterLines = 1;

% Loop over each node pair
for i = 1:length(propNBC.nodes)
    for j = i + 1:length(propNBC.nodes)
        % Get the node index in the element array
        nodeI = propNBC.nodes(i);
        nodeJ = propNBC.nodes(j);

        % Find the element indices to which the nodes belong
        [indexI, ~] = find(nodeI == strMsh.elements);
        [indexJ, ~] = find(nodeJ == strMsh.elements);

        % Find the common elements to which the nodes belong to
        [idComElmnt, ~] = intersect(indexI, indexJ);

        % Store the Neumann information on the common lines only if the
        % nodes on the Neumann boundary share one common elements
        if length(idComElmnt) == 1 && ...
            strcmp(propNBC.fctHandle(i, :), propNBC.fctHandle(j, :))
            propNBC.lines(counterLines,:) = ...
                [propNBC.nodes(i) propNBC.nodes(j) idComElmnt];
            fctHandle(counterLines,:) = propNBC.fctHandle(i,:);
            counterLines = counterLines + 1;
        end
    end
end

%% 12. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('\nParsing took %.2d seconds \n\n', computationalTime);
    fprintf('_________________________Parsing Ended__________________________\n');
    fprintf('################################################################\n\n\n');
end

end