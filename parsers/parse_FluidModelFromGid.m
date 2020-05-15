function [fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propALE, propNBC, ...
    propAnalysis, propParameters, propNLinearAnalysis, propFldDynamics, ...
    propGaussInt, postProc, propFSI] = ...
    parse_FluidModelFromGid(pathToCase, caseName, outMsg)
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
%              fldMsh : Structure containing information on the mesh,
%                           .nodes : The nodes in the FE mesh
%                        .elements : The elements in the FE mesh
%             homDOFs : The global numbering of the nodes where homogeneous
%                       Dirichlet boundary conditions are applied
%           inhomDOFs : The global numbering of the nodes where 
%                       inhomogeneous Dirichlet boundary conditions are 
%                       applied
%     valuesInhomDOFs : The prescribed values for the inhomogeneous 
%                       Dirichlet boundary conditions
%             propALE : The nodes on the ALE boundary
%                           .nodes : The coordinates of the nodes
%                       .fctHandle : The function handle to the computation
%                                    of the prescribed motion on the ALE
%                                    boundary nodes
%                 NBC : Structure on the Neumann boundary conditions,
%                           .nodes : The nodes where Neumann boundary 
%                                    conditions are applied
%                        .loadType : The type of the load for each Neumann 
%                                    node
%                       .fctHandle : The function handle for each Neumann 
%                                    node for the computation of the load 
%                                    vector (these functions are unde the 
%                                    folder load)
%        propAnalysis : Structure containing general information on the
%                       analysis,
%                                 .type      : The analysis type
%                       .noSpatialDimensions : Number of spatial dimensions
%                                 .noFields  : Number of DOFs per node
%      propParameters : Problem specific technical parameters,
%                              .nue : Dynamic viscosity
% propNLinearAnalysis : Structure containing information on the
%                       geometrically nonlinear analysis,
%                           .method : The employed nonlinear method
%                      .noLoadSteps : Number of load steps (typically used 
%                                     only in structural steady-state 
%                                     dynamics)
%                              .eps : The residual tolerance
%                          .maxIter : The maximum number of the nonlinear 
%                                     iterations
%     propFldDynamics : Structure containing information on the time
%                       integration regarding the fluid dynamics,
%                   .timeDependence : Steady-state or transient analysis
%                           .method : The time integration method
%                               .T0 : The start time of the simulation
%                             .TEnd : The end time of the simulation
%                      .noTimeSteps : The number of the time steps
%        propGaussInt : Structure containing information on the Gaussian
%                       quadrature,
%                                 .type : 'default', 'user'
%                           .domainNoGP : Number of Gauss Points for the 
%                                         domain integration
%                         .boundaryNoGP : Number of Gauss Points for the
%                                         boundary integration
%            postProc : Structure containing information on the
%                       postprocessing,
%                           .nameDomain : names of all the domains for post
%                          .nodesDomain : The global numbering of nodes
%                                         that are part of the domains above
%                      .computePostProc : function handles for calculation
%             propFSI : Structure containing information on Fluid-Structure
%                       interaction
%                       .coupledNodeIDs : Global numbering of the FSI nodes
%                      .numCoupledNodes : Number of FSI coupled nodes
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
% 14. Load the coupled fluid nodes for FSI
%
% 15. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________\n');
    fprintf('###########################################################\n');
    fprintf('Parsing data from GiD input file for a fluid boundary value\n');
    fprintf('problem has been initiated\n');
    fprintf('___________________________________________________________\n\n');
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
block = regexp(fstring, 'FLUID_ANALYSIS','split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
propAnalysis.type = out{1}{2};

% save the number of dimensions of the problem
if strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
        propAnalysis.noSpatialDimensions = 2;
        propAnalysis.noFields = 3;
elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
        propAnalysis.noSpatialDimensions = 3;
        propAnalysis.noFields = 4;
else
        error('Wrong analysis type selected');
end

if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Analysis type: %s \n', propAnalysis.type);
end

%% 3. Load the material properties
block = regexp(fstring, 'FLUID_MATERIAL_PROPERTIES', 'split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
propParameters.rho = str2double(out{1}{2});
propParameters.nue = str2double(out{1}{4});

%% 4. Load the nonlinear method
block = regexp(fstring, 'FLUID_NLINEAR_SCHEME', 'split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
propNLinearAnalysis.method = out{1}{2};
propNLinearAnalysis.noLoadSteps = str2double(out{1}{4});
propNLinearAnalysis.eps = str2double(out{1}{6});
propNLinearAnalysis.maxIter = str2double(out{1}{8});
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Nonlinear method: %s \n', propNLinearAnalysis.method);
    fprintf('\t>> No. load steps = %d \n', propNLinearAnalysis.noLoadSteps);
    fprintf('\t>> Convergence tolerance = %d \n', propNLinearAnalysis.eps);
    fprintf('\t>> Maximum number of iterations = %d \n', propNLinearAnalysis.maxIter);
end

%% 5. Load the time integration method
block = regexp(fstring, 'FLUID_TRANSIENT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
propFldDynamics.timeDependence = out{1}{2};
propFldDynamics.method = out{1}{4};
if strcmp(propFldDynamics.method, 'BOSSAK')
    propFldDynamics.alphaBeta =  str2double(out{1}{6});
    propFldDynamics.gamma =  str2double(out{1}{8});
end
propFldDynamics.T0 = str2double(out{1}{10});
propFldDynamics.TEnd = str2double(out{1}{12});
propFldDynamics.noTimeSteps = str2double(out{1}{14});
propFldDynamics.isAdaptive = out{1}{16};
propFldDynamics.dt = (propFldDynamics.TEnd - propFldDynamics.T0)/propFldDynamics.noTimeSteps;
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Fluid dynamics: %s \n', propFldDynamics.timeDependence);
    if ~strcmp(propFldDynamics.timeDependence, 'STEADY_STATE')
        fprintf('\t>> Time integration method: %s \n', propFldDynamics.method);
        if strcmp(propFldDynamics.method, 'BOSSAK')
            fprintf('\t \t>> alphaBeta =  %s \n', propFldDynamics.alphaBeta);
            fprintf('\t \t>> gamma =  %s \n', propFldDynamics.gamma);
        end
        fprintf('\t>> Start time of the simulation: %f \n', propFldDynamics.T0);
        fprintf('\t>> End time of the simulation: %f \n', propFldDynamics.TEnd);
        fprintf('\t>> Number of time steps: %f \n', propFldDynamics.noTimeSteps);
        fprintf('\t>> Time step size: %f \n', propFldDynamics.dt);
    end
end

%% 6. Load the Gauss Point integration method
block = regexp(fstring, 'FLUID_INTEGRATION', 'split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
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

%% 7. Load the fluid nodes
block = regexp(fstring,'FLUID_NODES', 'split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %f %f %f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
fldMsh.nodes = out(:,1:4);
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Number of nodes in the mesh: %d \n', length(fldMsh.nodes(:,1)));
end

%% 8. Load the fluid elements by connectivity arrays
block = regexp(fstring,'FLUID_ELEMENTS', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    if strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
        out{k} = textscan(block{k}, '%f %f %f %f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
        out{k} = textscan(block{k}, '%f %f %f %f %f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    else
        error('Wrong analysis type selected');
    end
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
fldMsh.elements = out(:, 1:end);
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Number of elements in the mesh: %d \n', length(fldMsh.elements));
end

%% 9. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
block = regexp(fstring, 'FLUID_DIRICHLET_NODES', 'split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);

% Filter out the actual DOFs
numDOFsNodeGiD = 4;
if strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
    numDOFsNode = 3;
    numDOFsNodeArray = [1 2 4];
elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
    numDOFsNode = 4;
    numDOFsNodeArray = [1 2 3 4];
else
	error('Wrong analysis type selected');
end

% Initialize counter
counterHomDBC = 1;
counterInhomDBC = 1;

% Find the number of nodes where Dirichlet boundary conditions are applied
noDBCNodes = length(out)/(numDOFsNodeGiD + 1);

for i = 1:noDBCNodes
    % Get the Dirichlet node ID
    nodeID = out((numDOFsNodeGiD + 1)*i - numDOFsNodeGiD);
    
    % Get the x-component of the prescribed value
    for jCounter = 1:length(numDOFsNodeArray)
        % Get the actual j-counter
        j = numDOFsNodeArray(jCounter);
        
        % Get the prescribed value at the current DOF
        presValue = out((numDOFsNodeGiD + 1)*i - numDOFsNodeGiD + j);
        if ~isnan(presValue)
            if presValue == 0
                homDOFs(counterHomDBC) = numDOFsNode*nodeID-numDOFsNode + jCounter;
                counterHomDBC = counterHomDBC + 1;
            else
                inhomDOFs(counterInhomDBC) = numDOFsNode*nodeID - numDOFsNode + jCounter;
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

%% 10. Load the nodes on which ALE conditions are applied
block = regexp(fstring, 'FLUID_DIRICHLET_ALE_NODES', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %s %d', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
end
if ~isempty(out)
    if strcmp(outMsg, 'outputEnabled')
        fprintf('>> Arbitrary Lagrangian-Eulerian method selected');
    end
    out = out{1};
    propALE.nodes = cell2mat(out(:, 1));
    outFctHandle = out(:, 2);
    propALE.fctHandle = outFctHandle{1};
    propALE.isFree = cell2mat(out(:, 3));
    fldMsh.initialNodes = fldMsh.nodes;
else
    propALE = [];
end

%% 11. Load the nodes, body names and function handles for post processing
block = regexp(fstring, 'FLUID_POST_PROC_NODES', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %s %s', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
end
if ~isempty(out)
    out = out{1};
    outNodes = out(:, 1);
    nodes = outNodes{1};
    outDomains = out(:, 2);
    domains = string(outDomains{1}); 
    outFctHandle = out(:, 3);
    fctHandles = string(outFctHandle{1});
    
    % get only the unique body names
    postProc.nameDomain = unique(domains)';
    
    % loop over the number of unique body names
    for k = 1:length(postProc.nameDomain)
        % find inxed of nodes
        indexArray = (domains == postProc.nameDomain(k));
        postProc.nodesDomain{k} = nodes(indexArray);
        % find corresponding function handle
        firstIndex = find(indexArray, 1);
        postProc.computePostProc(k) = fctHandles(firstIndex);
    end
else
    postProc = [];
end

%% 12. Load the nodes on the Neumann boundary together with the load application information
block = regexp(fstring, 'FLUID_FORCE_NODES', 'split');
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

%% 13. Get edge connectivity arrays for the Neumann edges
if ~isempty(out)
    if strcmp(outMsg, 'outputEnabled')
        fprintf('>> Neumann boundary edges: %d \n', length(propNBC.nodes) - 1);
    end
    
    % Initialize the Neumann boundary lines
    propNBC.lines = zeros(length(propNBC.nodes) - 1, 3);

    % Initialize line counter
    counterLines = 1;

    % Loop over each node pair
    for i = 1:length(propNBC.nodes)
        for j = i + 1:length(propNBC.nodes)
            % Get the node index in the element array
            nodeI = propNBC.nodes(i);
            nodeJ = propNBC.nodes(j);

            % Find the element indices to which the nodes belong
            [indexI, ~] = find(nodeI == fldMsh.elements(:,2:end));
            [indexJ, ~] = find(nodeJ == fldMsh.elements(:,2:end));

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
else
    propNBC = 'undefined';
end

%% 14. Load coupled fluid nodes for FSI
block = regexp(fstring, 'FLUID_FSI_NODES', 'split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f');
end
if ~isempty(out)
    out = out{1};
    propFSI.coupledNodeIDs = cell2mat(out(:, 1));
    propFSI.numCoupledNodes = length(propFSI.coupledNodeIDs);
else
    propFSI.coupledNodeIDs = [];
    propFSI.numCoupledNodes = 0;
end

%% 15. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('\nParsing took %.2d seconds \n\n', computationalTime);
    fprintf('_______________________Parsing Ended_______________________\n');
    fprintf('###########################################################\n\n\n');
end

end
