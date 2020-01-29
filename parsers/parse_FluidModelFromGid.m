function [fldMsh,homDBC,inhomDBC,valuesInhomDBC,propALE,NBC,analysis,...
    parameters,propNLinearAnalysis,propFldDynamics,gaussInt,postProc] = ...
    parse_FluidModelFromGid(pathToCase,caseName,outMsg)
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
% Parses data from an input file created using GiD for a fluid boundary 
% value problem.
%
%                Input :
%           pathToCase : The absolute path to the inputGiD case folder
%             caseName : The name of the case in the inputGiD case folder
%               outMsg : On the output information on the command window
%
%              Output :
%              fldMsh :     .nodes : The nodes in the FE mesh
%                        .elements : The elements in the FE mesh
%              homDBC : The global numbering of the nodes where homogeneous
%                       Dirichlet boundary conditions are applied
%            inhomDBC : The global numbering of the nodes where 
%                       inhomogeneous Dirichlet boundary conditions are 
%                       applied
%      valuesInhomDBC : The prescribed values for the inhomogeneous 
%                       Dirichlet boundary conditions
%             propALE : The nodes on the ALE boundary
%                           .nodes : The coordinates of the nodes
%                       .fctHandle : The function handle to the computation
%                                    of the prescribed motion on the ALE
%                                    boundary nodes
%                 NBC :     .nodes : The nodes where Neumann boundary 
%                                    conditions are applied
%                        .loadType : The type of the load for each Neumann 
%                                    node
%                       .fctHandle : The function handle for each Neumann 
%                                    node for the computation of the load 
%                                    vector (these functions are unde the 
%                                    folder load)
%            analysis : .type      : The analysis type
%             .noSpatialDimensions : Number of spatial dimensions
%                       .noFields  : Number of DOFs per node
%          parameters : Problem specific technical parameters
% propNLinearAnalysis :     .method : The employed nonlinear method
%                      .noLoadSteps : Number of load steps (typically used 
%                                     only in structural steady-state 
%                                     dynamics)
%                              .eps : The residual tolerance
%                          .maxIter : The maximum number of the nonlinear 
%                                     iterations
%     propFldDynamics : .timeDependence : Steady-state or transient analysis
%                           .method : The time integration method
%                               .T0 : The start time of the simulation
%                             .TEnd : The end time of the simulation
%                      .noTimeSteps : The number of the time steps
%            gaussInt : On the Gauss Point integration
%                                 .type : 'default', 'user'
%                           .domainNoGP : Number of Gauss Points for the 
%                                         domain integration
%                         .boundaryNoGP : Number of Gauss Points for the
%                                         boundary integration
%            postProc : Post-processing properties
%                           .nameDomain : names of all the domains for post
%                          .nodesDomain : The global numbering of nodes
%                                         that are part of the domains above
%                      .computePostProc : function handles for calculation
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
% 10. Load the nodes on which ALE conditions are applied
%
% 11. Load the nodes, body names and function handles for post processing
%
% 12. Load the nodes on the Neumann boundary together with the load application information
% 
% 13. Get edge connectivity arrays for the Neumann edges
%
% 14. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________\n');
    fprintf('###########################################################\n');
    fprintf('Parsing data from GiD input file for a fluid boundary value\n');
    fprintf('problem has been initiated\n');
    fprintf('___________________________________________________________\n\n');

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
block = regexp(fstring,'FLUID_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
analysis.type = out{1}{2};

% save the number of dimensions of the problem
if strcmp(analysis.type,'NAVIER_STOKES_2D')
        analysis.noSpatialDimensions = 2;
        analysis.noFields = 3;
elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
        analysis.noSpatialDimensions = 3;
        analysis.noFields = 4;
else
        error('Wrong analysis type selected');
end

if strcmp(outMsg,'outputEnabled')
    fprintf('>> Analysis type: %s \n',analysis.type);
end

%% 3. Load the material properties
block = regexp(fstring,'FLUID_MATERIAL_PROPERTIES','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
parameters.rho = str2double(out{1}{2});
parameters.nue = str2double(out{1}{4});

%% 4. Load the nonlinear method
block = regexp(fstring,'FLUID_NLINEAR_SCHEME','split');
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
block = regexp(fstring,'FLUID_TRANSIENT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',' ','MultipleDelimsAsOne', 1);
propFldDynamics.timeDependence = out{1}{2};
propFldDynamics.method = out{1}{4};
if strcmp(propFldDynamics.method,'bossak')
    propFldDynamics.alphaBeta =  str2double(out{1}{6});
    propFldDynamics.gamma =  str2double(out{1}{8});
end
propFldDynamics.T0 = str2double(out{1}{10});
propFldDynamics.TEnd = str2double(out{1}{12});
propFldDynamics.noTimeSteps = str2double(out{1}{14});
propFldDynamics.isAdaptive = out{1}{16};
propFldDynamics.dt = (propFldDynamics.TEnd - propFldDynamics.T0)/propFldDynamics.noTimeSteps;
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Fluid dynamics: %s \n',propFldDynamics.timeDependence);
    if ~strcmp(propFldDynamics.timeDependence,'STEADY_STATE')
        fprintf('\t>> Time integration method: %s \n',propFldDynamics.method);
        if strcmp(propFldDynamics.method,'bossak')
            fprintf('\t \t>> alphaBeta =  %s \n',propFldDynamics.alphaBeta);
            fprintf('\t \t>> gamma =  %s \n',propFldDynamics.gamma);
        end
        fprintf('\t>> Start time of the simulation: %f \n',propFldDynamics.T0);
        fprintf('\t>> End time of the simulation: %f \n',propFldDynamics.TEnd);
        fprintf('\t>> Number of time steps: %f \n',propFldDynamics.noTimeSteps);
        fprintf('\t>> Time step size: %f \n',propFldDynamics.dt);
    end
end

%% 6. Load the Gauss Point integration method
block = regexp(fstring,'FLUID_INTEGRATION','split');
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

%% 7. Load the fluid nodes
block = regexp(fstring,'FLUID_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
fldMsh.nodes = out(:,2:4);
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Number of nodes in the mesh: %d \n',length(fldMsh.nodes(:,1)));
end

%% 8. Load the fluid elements by connectivity arrays
block = regexp(fstring,'FLUID_ELEMENTS','split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    if strcmp(analysis.type,'NAVIER_STOKES_2D')
        out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
        out{k} = textscan(block{k},'%f %f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    else
        error('Wrong analysis type selected');
    end
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
fldMsh.elements = out(:,2:end);
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Number of elements in the mesh: %d \n',length(fldMsh.elements));
end

%% 9. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
block = regexp(fstring,'FLUID_DIRICHLET_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);

% Filter out the actual DOFs
noDOFsPerNodeFromGiD = 4;
if strcmp(analysis.type,'NAVIER_STOKES_2D')
    noDOFsPerNode = 3;
    noDOFsPerNodeArray = [1 2 4];
elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
    noDOFsPerNode = 4;
    noDOFsPerNodeArray = [1 2 3 4];
else
	error('Wrong analysis type selected');
end

% Initialize counter
counterHomDBC = 1;
counterInhomDBC = 1;

% Find the number of nodes where Dirichlet boundary conditions are applied
noDBCNodes = length(out)/(noDOFsPerNodeFromGiD+1);

for i = 1:noDBCNodes
    % Get the Dirichlet node ID
    nodeID = out((noDOFsPerNodeFromGiD+1)*i-noDOFsPerNodeFromGiD);
    
    % Get the x-component of the prescribed value
    for jCounter = 1:length(noDOFsPerNodeArray)
        % Get the actual j-counter
        j = noDOFsPerNodeArray(jCounter);
        
        % Get the prescribed value at the current DOF
        presValue = out((noDOFsPerNodeFromGiD+1)*i-noDOFsPerNodeFromGiD+j);
        if ~isnan(presValue)
            if presValue == 0
                homDBC(counterHomDBC) = noDOFsPerNode*nodeID-noDOFsPerNode+j;
                counterHomDBC = counterHomDBC + 1;
            else
                inhomDBC(counterInhomDBC) = noDOFsPerNode*nodeID-noDOFsPerNode+j;
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

%% 10. Load the nodes on which ALE conditions are applied
block = regexp(fstring,'FLUID_DIRICHLET_ALE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %s %s','delimiter',' ','MultipleDelimsAsOne', 1);
end
if ~isempty(out)
    out = out{1};
    propALE.nodes = cell2mat(out(:,1));
    outFctHandle = out(:,2);
    propALE.fctHandle = cell2mat(outFctHandle{1});
    fldMsh.initialNodes = fldMsh.nodes;
else
    propALE = [];
end

if strcmp(outMsg,'outputEnabled') && ~isempty(out)
    fprintf('>> Arbitrary Lagrangian-Eulerian method selected');
end

%% 11. Load the nodes, body names and function handles for post processing
block = regexp(fstring,'FLUID_POST_PROC_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %s %s','delimiter',' ','MultipleDelimsAsOne', 1);
end

if ~isempty(out)
    out = out{1};
    outNodes = out(:,1);
    nodes = outNodes{1};
    outDomains = out(:,2);
    domains = string(outDomains{1}); 
    outFctHandle = out(:,3);
    fctHandles = string(outFctHandle{1});
    
    % get only the unique body names
    postProc.nameDomain = unique(domains)';
    
    % loop over the number of unique body names
    for k = 1:length(postProc.nameDomain)
        % find inxed of nodes
        indexArray = (domains == postProc.nameDomain(k));
        postProc.nodesDomain{k} = nodes(indexArray);
        % find corresponding function handle
        firstIndex = find(indexArray,1);
        postProc.computePostProc(k) = fctHandles(firstIndex);
    end
    
else
    postProc = [];
end

%% 12. Load the nodes on the Neumann boundary together with the load application information
block = regexp(fstring,'FLUID_FORCE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %s %s','delimiter',' ','MultipleDelimsAsOne', 1);
end

if ~isempty(out)
    out = out{1};
    NBC.nodes = cell2mat(out(:,1));
    outLoadType = out(:,2);
    NBC.loadType = cell2mat(outLoadType{1});
    outFctHandle = out(:,3);
    NBC.fctHandle = cell2mat(outFctHandle{1});
end

%% 13. Get edge connectivity arrays for the Neumann edges
if ~isempty(out)
    if strcmp(outMsg,'outputEnabled')
        fprintf('>> Neumann boundary edges: %d \n',length(NBC.nodes)-1);
    end
    
    % Initialize the Neumann boundary lines
    NBC.lines = zeros(length(NBC.nodes)-1,3);

    % Initialize line counter
    counterLines = 1;

    % Loop over each node pair
    for i=1:length(NBC.nodes)
        for j=i:length(NBC.nodes)
            % If we are not in the same node
            if i~=j
                % Get the node index in the element array
                nodeI = NBC.nodes(i);
                nodeJ = NBC.nodes(j);


                % Find the element indices to which the nodes belong
                [indexI,~] = find(nodeI == fldMsh.elements);
                [indexJ,~] = find(nodeJ == fldMsh.elements);

                % For all the element indices to which indexJ belongs to
                for k=1:length(indexJ)
                    % Find the common elements to which both nodes belong to
                    commonElmnts = find(indexJ(k) == indexI);

                    % If there are commont elements to which the nodes belong
                    if norm(commonElmnts)~=0
                        % Get the common element index
                        elementIndex = indexI(commonElmnts);

                        % Store the line into the NBC.line array with the same
                        % ordering as the are stored in the element array
                        NBC.lines(counterLines,:) = ...
                            [NBC.nodes(i) NBC.nodes(j) elementIndex];

                        % Update counter
                        counterLines = counterLines + 1;
                    end
                end
            end
        end
    end
else
    NBC = 'undefined';
end

%% 14. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nParsing took %.2d seconds \n\n',computationalTime);
    fprintf('_______________________Parsing Ended_______________________\n');
    fprintf('###########################################################\n\n\n');
end

end
