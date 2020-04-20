function [strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, ...
    propAnalysis, parameters, propNLinearAnalysis, propHeatDynamics, ...
    propGaussInt] = parse_HeatModelFromGid(pathToCase, caseName, outMsg)
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
%             propNBC :     .nodes : The nodes where Neumann boundary 
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
%          parameters : Problem specific technical (physical) parameters
%                       (.rho,.E,.nue,.t)
% propNLinearAnalysis :     .method : The employed nonlinear method
%                        .tolerance : The residual tolerance
%                          .maxIter : The maximum number of the nonlinear 
%                                     iterations
%     propStrDynamics : .timeDependence : On the transient analysis
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
% 4. Load the nonlinear method
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
    fprintf('________________________________________________________________\n');
    fprintf('################################################################\n');
    fprintf('Parsing data from GiD input file for a structural boundary value\n');
    fprintf('problem has been initiated\n');
    fprintf('________________________________________________________________\n\n');

    % start measuring computational time
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
block = regexp(fstring,'HEAT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
propAnalysis.type = out{1}{2};
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Analysis type: %s \n',propAnalysis.type);
end

%% 3. Load the material properties
block = regexp(fstring,'HEAT_MATERIAL_PROPERTIES','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
parameters.rho = str2double(out{1}{2}); % density
parameters.k = str2double(out{1}{4});   % thermal conductivity
parameters.cp = str2double(out{1}{6});  % specific heat
% compute thermal diffusivity
parameters.alpha = parameters.k / (parameters.rho*parameters.cp);

%% 4. Load the nonlinear method
% heat transfer without internal heat generation is always linear
propNLinearAnalysis.method = 'UNDEFINED';
propNLinearAnalysis.noLoadSteps = [];
propNLinearAnalysis.eps = [];
propNLinearAnalysis.maxIter = 1;

%% 5. Load the time integration method
block = regexp(fstring,'HEAT_TRANSIENT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',' ','MultipleDelimsAsOne', 1);
propHeatDynamics.timeDependence = out{1}{2};
if ~strcmp(propHeatDynamics.timeDependence,'STEADY_STATE')
    propHeatDynamics.method = out{1}{4};
    propHeatDynamics.T0 = str2double(out{1}{6});
    propHeatDynamics.TEnd = str2double(out{1}{8});
    propHeatDynamics.noTimeSteps = str2double(out{1}{10});
    propHeatDynamics.isAdaptive = out{1}{12};
    propHeatDynamics.dt = (propHeatDynamics.TEnd - propHeatDynamics.T0)/propHeatDynamics.noTimeSteps;
end

if strcmp(outMsg,'outputEnabled')
    fprintf('>> Structural dynamics: %s \n',propHeatDynamics.timeDependence);
    if ~strcmp(propHeatDynamics.timeDependence,'STEADY_STATE')
        fprintf('\t>> Time integration method: %s \n',propHeatDynamics.method);
        fprintf('\t>> Start time of the simulation: %f \n',propHeatDynamics.T0);
        fprintf('\t>> End time of the simulation: %f \n',propHeatDynamics.TEnd);
        fprintf('\t>> Number of time steps: %f \n',propHeatDynamics.noTimeSteps);
        fprintf('\t>> Time step size: %f \n',propHeatDynamics.dt);
    end
end

%% 6. Load the Gauss Point integration method
block = regexp(fstring,'HEAT_INTEGRATION','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',' ','MultipleDelimsAsOne', 1);
propGaussInt.type = out{1}{2};
if strcmp(propGaussInt.type,'user')
    propGaussInt.domainNoGP = str2double(out{1}{4});
    propGaussInt.boundaryNoGP = str2double(out{1}{6});
end
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Gauss integration type: %s \n',propGaussInt.type);
    if strcmp(propGaussInt.type,'user')
        fprintf('\t>> No. Gauss Points for the domain integration: %d \n',propGaussInt.domainNoGP);
        fprintf('\t>> No. Gauss Points for the boundary integration: %d \n',propGaussInt.boundaryNoGP);
    end
end

%% 7. Load the structural nodes
block = regexp(fstring,'HEAT_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
strMsh.nodes = out(:,2:4);
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Number of nodes in the mesh: %d \n',length(strMsh.nodes(:,1)));
end

%% 8. Load the structural elements by connectivity arrays
block = regexp(fstring,'HEAT_ELEMENTS','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
noNodes = length(out(1,:));
strMsh.elements = out(:,2:noNodes);
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Number of elements in the mesh: %d \n',length(strMsh.elements));
end

%% 9. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
block = regexp(fstring,'HEAT_DIRICHLET_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);

% Filter out the actual DOFs
nDOFsPerNodeFromGiD = 1;
if strcmp(propAnalysis.type,'HEAT_TRANSFER_2D')
    nDOFsPerNode = 1;
end

% Initialize counter
counterHomDBC = 1;
counterInhomDBC = 1;

% Find the number of nodes where Dirichlet boundary conditions are applied
nDBCNodes = length(out)/(nDOFsPerNodeFromGiD+1);

for i = 1:nDBCNodes
    % Get the Dirichlet node ID
    nodeID = out((nDOFsPerNodeFromGiD+1)*i-nDOFsPerNodeFromGiD);
    
    % Get the x-component of the prescribed value
    for j = 1:nDOFsPerNode
        presValue = out((nDOFsPerNodeFromGiD+1)*i-nDOFsPerNodeFromGiD+j);
        if ~isnan(presValue)
            if presValue == 0
                homDOFs(counterHomDBC) = nDOFsPerNode*nodeID-nDOFsPerNode+j;
                counterHomDBC = counterHomDBC + 1;
            else
                inhomDOFs(counterInhomDBC) = nDOFsPerNode*nodeID-nDOFsPerNode+j;
                valuesInhomDOFs(counterInhomDBC) = presValue;
                counterInhomDBC = counterInhomDBC + 1;
            end
        end
    end
end

% Sort out the vectors
homDOFs = sort(homDOFs);
[inhomDOFs,indexSorting] = sort(inhomDOFs);
valuesInhomDOFs = valuesInhomDOFs(indexSorting);

%% 10. Load the nodes on the Neumann boundary together with the load application information
block = regexp(fstring,'HEAT_FLUX_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %s %s','delimiter',' ','MultipleDelimsAsOne', 1);
end

if ~isempty(out)
    out = out{1};
    propNBC.nodes = cell2mat(out(:,1));
    outLoadType = out(:,2);
    propNBC.loadType = cell2mat(outLoadType{1});
    outFctHandle = out(:,3);
    propNBC.fctHandle = cell2mat(outFctHandle{1});
end

%% 11. Get edge connectivity arrays for the Neumann edges

if strcmp(outMsg,'outputEnabled')
    fprintf('>> Neumann boundary edges: %d \n',length(propNBC.nodes) - 1);
end

% Initialize the Neumann boundary lines
propNBC.lines = zeros(length(unique(propNBC.nodes)) - 1,3);

% Initialize line counter
counterLines = 1;

% Loop over each node pair
for i = 1:length(propNBC.nodes)
    for j = i:length(propNBC.nodes)
        % If we are not in the same node
        if i ~= j
            % Get the node index in the element array
            nodeI = propNBC.nodes(i);
            nodeJ = propNBC.nodes(j);


            % Find the element indices to which the nodes belong
            [indexI,~] = find(nodeI == strMsh.elements);
            [indexJ,~] = find(nodeJ == strMsh.elements);

            % For all the element indices to which indexJ belongs to
            for k = 1:length(indexJ)
                % Find the common elements to which both nodes belong to
                commonElmnts = find(indexJ(k) == indexI);

                % If there are common elements to which the nodes belong
                if norm(commonElmnts) ~= 0 && strcmp(propNBC.fctHandle(i,:),propNBC.fctHandle(j,:))
                    % Get the common element index
                    elementIndex = indexI(commonElmnts);

                    % Store the line into the NBC.line array with the same
                    % ordering as the are stored in the element array
                    propNBC.lines(counterLines,:) = ...
                        [propNBC.nodes(i) propNBC.nodes(j) elementIndex];

                    % Get the appropriate function handle for the line
                    fctHandle(counterLines,:) = propNBC.fctHandle(i,:);

                    % Update counter
                    counterLines = counterLines + 1;
                end
            end
        end
    end
end
%propNBC.fctHandle = fctHandle;
%% 12. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nParsing took %.2d seconds \n\n',computationalTime);
    fprintf('_________________________Parsing Ended__________________________\n');
    fprintf('################################################################\n\n\n');
end

end
