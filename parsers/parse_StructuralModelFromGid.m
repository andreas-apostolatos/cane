function [strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,...
    propNLinearAnalysis,propStrDynamics,gaussInt] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
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
%              homDBC : The global numbering of the nodes where homogeneous
%                       Dirichlet boundary conditions are applied
%            inhomDBC : The global numbering of the nodes where 
%                       inhomogeneous Dirichlet boundary conditions are 
%                       applied
%      valuesInhomDBC : The prescribed values for the inhomogeneous 
%                       Dirichlet boundary conditions
%                 NBC :     .nodes : The nodes where Neumann boundary 
%                                    conditions are applied
%                        .loadType : The type of the load for each Neumann 
%                                    node
%                       .fctHandle : The function handle for each Neumann 
%                                    node for the computation of the load 
%                                    vector (these functions are under the 
%                                    folder load)
%            analysis : On the analysis type
%                             .type : The analysis type
%          parameters : Problem specific technical parameters
% propNLinearAnalysis :     .method : The employed nonlinear method
%                        .tolerance : The residual tolerance
%                          .maxIter : The maximum number of the nonlinear 
%                                     iterations
%     propStrDynamics : .timeDependence : On the transient analysis
%                               .method : The time integration method
%                                   .T0 : The start time of the simulation
%                                 .TEnd : The end time of the simulation
%                          .noTimeSteps : The number of the time steps
%            gaussInt : On the Gauss Point integration
%                               .type : 'default', 'user'
%                         .domainNoGP : Number of Gauss Points for the 
%                                       domain integration
%                       .boundaryNoGP : Number of Gauss Points for the 
%                                       boundary integration
%
% Function layout :
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
% 7. Load the structural nodes
%
% 8. Load the structural elements by connectivity arrays
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
homDBC = [];
inhomDBC = [];
valuesInhomDBC = [];

%% 1. Load the input file from GiD
fstring = fileread([pathToCase caseName '.dat']); 

%% 2. Load the analysis type
block = regexp(fstring,'STRUCTURE_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
analysis.type = out{1}{2};
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Analysis type: %s \n',analysis.type);
end

%% 3. Load the material properties
block = regexp(fstring,'STRUCTURE_MATERIAL_PROPERTIES','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
parameters.rho = str2double(out{1}{2});
parameters.E = str2double(out{1}{4});
parameters.nue = str2double(out{1}{6});

%% 4. Load the nonlinear method
block = regexp(fstring,'STRUCTURE_NLINEAR_SCHEME','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
propNLinearAnalysis.method = out{1}{2};
propNLinearAnalysis.noLoadSteps = str2double(out{1}{4});
propNLinearAnalysis.eps = str2double(out{1}{6});
propNLinearAnalysis.maxIter = str2double(out{1}{8});
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Nonlinear method: %s \n',propNLinearAnalysis.method);
    fprintf('\t>> No. load steps = %d \n',propNLinearAnalysis.noLoadSteps);
    fprintf('\t>> Convergence tolerance = %d \n',propNLinearAnalysis.eps);
    fprintf('\t>> Maximum number of iterations = %d \n',propNLinearAnalysis.maxIter);
end

%% 5. Load the time integration method
block = regexp(fstring,'STRUCTURE_TRANSIENT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',' ','MultipleDelimsAsOne', 1);
propStrDynamics.timeDependence = out{1}{2};
propStrDynamics.method = out{1}{4};
if strcmp(propStrDynamics.method,'bossak')
    propStrDynamics.alphaBeta =  str2double(out{1}{6});
    propStrDynamics.gamma =  str2double(out{1}{8});
end
propStrDynamics.T0 = str2double(out{1}{10});
propStrDynamics.TEnd = str2double(out{1}{12});
propStrDynamics.noTimeSteps = str2double(out{1}{14});
propStrDynamics.isAdaptive = out{1}{16};
propStrDynamics.dt = (propStrDynamics.TEnd - propStrDynamics.T0)/propStrDynamics.noTimeSteps;
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Structural dynamics: %s \n',propStrDynamics.timeDependence);
    if ~strcmp(propStrDynamics.timeDependence,'STEADY_STATE')
        fprintf('\t>> Time integration method: %s \n',propStrDynamics.method);
        if strcmp(propStrDynamics.method,'bossak')
            fprintf('\t \t>> alphaBeta =  %s \n',propStrDynamics.alphaBeta);
            fprintf('\t \t>> gamma =  %s \n',propStrDynamics.gamma);
        end
        fprintf('\t>> Start time of the simulation: %f \n',propStrDynamics.T0);
        fprintf('\t>> End time of the simulation: %f \n',propStrDynamics.TEnd);
        fprintf('\t>> Number of time steps: %f \n',propStrDynamics.noTimeSteps);
        fprintf('\t>> Time step size: %f \n',propStrDynamics.dt);
    end
end

%% 6. Load the Gauss Point integration method
block = regexp(fstring,'STRUCTURE_INTEGRATION','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',' ','MultipleDelimsAsOne', 1);
gaussInt.type = out{1}{2};
if strcmp(gaussInt.type,'user')
    gaussInt.domainNoGP = str2double(out{1}{4});
    gaussInt.boundaryNoGP = str2double(out{1}{6});
end
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Gauss integration type: %s \n',gaussInt.type);
    if strcmp(gaussInt.type,'user')
        fprintf('\t>> No. Gauss Points for the domain integration: %d \n',gaussInt.domainNoGP);
        fprintf('\t>> No. Gauss Points for the boundary integration: %d \n',gaussInt.boundaryNoGP);
    end
end

%% 7. Load the structural nodes
block = regexp(fstring,'STRUCTURE_NODES','split'); 
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
block = regexp(fstring,'STRUCTURE_ELEMENTS','split'); 
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
block = regexp(fstring,'STRUCTURE_DIRICHLET_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);

% Filter out the actual DOFs
nDOFsPerNodeFromGiD = 3;
if strcmp(analysis.type,'planeStress')||strcmp(analysis.type,'planeStrain')
    nDOFsPerNode = 2;
end

% Initialize counter
counterHomDBC = 1;
counterInhomDBC = 1;

% Find the number of nodes where Dirichlet boundary conditions are applied
nDBCNodes = length(out)/(nDOFsPerNodeFromGiD+1);

for i=1:nDBCNodes
    % Get the Dirichlet node ID
    nodeID = out((nDOFsPerNodeFromGiD+1)*i-nDOFsPerNodeFromGiD);
    
    % Get the x-component of the prescribed value
    for j = 1:nDOFsPerNode
        presValue = out((nDOFsPerNodeFromGiD+1)*i-nDOFsPerNodeFromGiD+j);
        if ~isnan(presValue)
            if presValue == 0
                homDBC(counterHomDBC) = nDOFsPerNode*nodeID-nDOFsPerNode+j;
                counterHomDBC = counterHomDBC + 1;
            else
                inhomDBC(counterInhomDBC) = nDOFsPerNode*nodeID-nDOFsPerNode+j;
                valuesInhomDBC(counterInhomDBC) = presValue;
                counterInhomDBC = counterInhomDBC + 1;
            end
        end
    end
end

% Sort out the vectors
homDBC = sort(homDBC);
[inhomDBC,indexSorting] = sort(inhomDBC);
valuesInhomDBC = valuesInhomDBC(indexSorting);

%% 10. Load the nodes on the Neumann boundary together with the load application information
block = regexp(fstring,'STRUCTURE_FORCE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %s %s','delimiter',' ','MultipleDelimsAsOne', 1);
end
out = out{1};
NBC.nodes = cell2mat(out(:,1));
outLoadType = out(:,2);
NBC.loadType = cell2mat(outLoadType{1});
outFctHandle = out(:,3);
NBC.fctHandle = cell2mat(outFctHandle{1});

%% 11. Get edge connectivity arrays for the Neumann edges
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Neumann boundary edges: %d \n',length(NBC.nodes) - 1);
end

% Initialize the Neumann boundary lines
NBC.lines = zeros(length(unique(NBC.nodes)) - 1,3);

% Initialize line counter
counterLines = 1;

% Loop over each node pair
for i = 1:length(NBC.nodes)
    for j = i:length(NBC.nodes)
        % If we are not in the same node
        if i ~= j
            % Get the node index in the element array
            nodeI = NBC.nodes(i);
            nodeJ = NBC.nodes(j);
            
            
            % Find the element indices to which the nodes belong
            [indexI,~] = find(nodeI == strMsh.elements);
            [indexJ,~] = find(nodeJ == strMsh.elements);
            
            % For all the element indices to which indexJ belongs to
            for k = 1:length(indexJ)
                % Find the common elements to which both nodes belong to
                commonElmnts = find(indexJ(k) == indexI);
                
                % If there are commont elements to which the nodes belong
                if norm(commonElmnts) ~= 0 && strcmp(NBC.fctHandle(i,:),NBC.fctHandle(j,:))
                    % Get the common element index
                    elementIndex = indexI(commonElmnts);
                    
                    % Store the line into the NBC.line array with the same
                    % ordering as the are stored in the element array
                    NBC.lines(counterLines,:) = ...
                        [NBC.nodes(i) NBC.nodes(j) elementIndex];
                    
                    % Get the appropriate function handle for the line
                    fctHandle(counterLines,:) = NBC.fctHandle(i,:);
                    
                    % Update counter
                    counterLines = counterLines + 1;
                end
            end
        end
    end
end
NBC.fctHandle = fctHandle;

%% 12. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nParsing took %.2d seconds \n\n',computationalTime);
    fprintf('_________________________Parsing Ended__________________________\n');
    fprintf('################################################################\n\n\n');
end

end
