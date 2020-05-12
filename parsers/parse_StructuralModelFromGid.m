function [strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, ...
    propAnalysis, propParameters, propNLinearAnalysis, propStrDynamics, ...
    propGaussInt, propContact, propFSI] = ...
    parse_StructuralModelFromGid(pathToCase, caseName, outMsg)
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
%              strMsh : Structure containing information on the mesh,
%                           .nodes : The nodes in the FE mesh
%                         .nodeIDs : Global numbering of the nodes
%                        .elements : The elements in the FE mesh
%                      .elementIDs : Global numbering of the elements
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
%                              .rho : Material density,
%                                .E : Young's modulus
%                              .nue : Poisson ration 
%                                .t : thickness
% propNLinearAnalysis : Structure containing information on the
%                       geometrically nonlinear analysis,
%                           .method : The employed nonlinear method
%                        .tolerance : The residual tolerance
%                          .maxIter : The maximum number of the nonlinear 
%                                     iterations
%     propStrDynamics : Structure containing information on the time
%                       integration regarding the structural dynamics,
%                       .timeDependence : On the transient analysis
%                               .method : The time integration method
%                                   .T0 : The start time of the simulation
%                                 .TEnd : The end time of the simulation
%                          .noTimeSteps : The number of the time steps
%        propGaussInt : Structure containing information on the Gaussian
%                       quadrature,
%                               .type : 'default', 'user'
%                         .domainNoGP : Number of Gauss Points for the 
%                                       domain integration
%                       .boundaryNoGP : Number of Gauss Points for the 
%                                       boundary integration
%         propContact : Structure containing information on the nodes which
%                       are going to be in potential contact
%                            .nodeIDs : Global numbering of contact nodes
%                       numberOfNodes : number of contact nodes
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
% 10. Load the nodes on the Neumann boundary together with the load application information
%
% 11. Load the nodes that are candidates for contact
%
% 12. Load the coupled structure nodes for FSI
%
% 13. Get edge connectivity arrays for the Neumann edges
%
% 14. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('________________________________________________________________\n');
    fprintf('################################################################\n');
    fprintf('Parsing data from GiD input file for a structural boundary value\n');
    fprintf('problem has been initiated\n');
    fprintf('________________________________________________________________\n\n');
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
block = regexp(fstring, 'STRUCTURE_ANALYSIS','split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
propAnalysis.type = out{1}{2};
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Analysis type: %s \n', propAnalysis.type);
end

%% 3. Load the material properties
block = regexp(fstring, 'STRUCTURE_MATERIAL_PROPERTIES', 'split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ',', 'MultipleDelimsAsOne', 1);
propParameters.rho = str2double(out{1}{2});
propParameters.E = str2double(out{1}{4});
propParameters.nue = str2double(out{1}{6});
propParameters.t = str2double(out{1}{8});

%% 4. Load the nonlinear method
block = regexp(fstring, 'STRUCTURE_NLINEAR_SCHEME', 'split');
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
block = regexp(fstring, 'STRUCTURE_TRANSIENT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1}, '%s', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
propStrDynamics.timeDependence = out{1}{2};
propStrDynamics.method = out{1}{4};
if strcmp(propStrDynamics.method, 'bossak')
    propStrDynamics.alphaBeta =  str2double(out{1}{6});
    propStrDynamics.gamma =  str2double(out{1}{8});
end
propStrDynamics.T0 = str2double(out{1}{10});
propStrDynamics.TEnd = str2double(out{1}{12});
propStrDynamics.noTimeSteps = str2double(out{1}{14});
propStrDynamics.isAdaptive = out{1}{16};
propStrDynamics.dt = (propStrDynamics.TEnd - propStrDynamics.T0)/propStrDynamics.noTimeSteps;
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Structural dynamics: %s \n', propStrDynamics.timeDependence);
    if ~strcmp(propStrDynamics.timeDependence, 'STEADY_STATE')
        fprintf('\t>> Time integration method: %s \n', propStrDynamics.method);
        if strcmp(propStrDynamics.method, 'bossak')
            fprintf('\t \t>> alphaBeta =  %s \n', propStrDynamics.alphaBeta);
            fprintf('\t \t>> gamma =  %s \n', propStrDynamics.gamma);
        end
        fprintf('\t>> Start time of the simulation: %f \n', propStrDynamics.T0);
        fprintf('\t>> End time of the simulation: %f \n', propStrDynamics.TEnd);
        fprintf('\t>> Number of time steps: %f \n', propStrDynamics.noTimeSteps);
        fprintf('\t>> Time step size: %f \n', propStrDynamics.dt);
    end
end

%% 6. Load the Gauss Point integration method
block = regexp(fstring, 'STRUCTURE_INTEGRATION', 'split');
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

%% 7. Load the structural nodes
block = regexp(fstring, 'STRUCTURE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %f %f %f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
strMsh.nodes = out(:, 2:4);
strMsh.nodeIDs = out(:,1);
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Number of nodes in the mesh: %d \n', length(strMsh.nodes(:, 1)));
end

%% 8. Load the structural elements by connectivity arrays
block = regexp(fstring, 'STRUCTURE_ELEMENTS', 'split');
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %f %f %f %f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
noNodes = length(out(1, :));
strMsh.elements = out(:, 2:noNodes);
strMsh.elementIDs = out(:, 1);
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Number of elements in the mesh: %d \n', length(strMsh.elements));
end

%% 9. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
block = regexp(fstring, 'STRUCTURE_DIRICHLET_NODES', 'split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f', 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);

% Filter out the actual DOFs
numDOFsNodeGiD = 3;
if strcmp(propAnalysis.type, 'planeStress') || ...
        strcmp(propAnalysis.type, 'planeStrain')
    numDOFsPerNode = 2;
end

% Initialize counter
counterHomDBC = 1;
counterInhomDBC = 1;

% Find the number of nodes where Dirichlet boundary conditions are applied
numDBCNodes = length(out)/(numDOFsNodeGiD + 1);

for i = 1:numDBCNodes
    % Get the Dirichlet node ID
    nodeID = out((numDOFsNodeGiD + 1)*i-numDOFsNodeGiD);
    
    % Get the x-component of the prescribed value
    for j = 1:numDOFsPerNode
        presValue = out((numDOFsNodeGiD + 1)*i-numDOFsNodeGiD + j);
        if ~isnan(presValue)
            if presValue == 0
                homDOFs(counterHomDBC) = numDOFsPerNode*nodeID - numDOFsPerNode + j;
                counterHomDBC = counterHomDBC + 1;
            else
                inhomDOFs(counterInhomDBC) = numDOFsPerNode*nodeID-numDOFsPerNode + j;
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
block = regexp(fstring, 'STRUCTURE_FORCE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f %s %s','delimiter', ' ', 'MultipleDelimsAsOne', 1);
end
out = out{1};
propNBC.nodes = cell2mat(out(:, 1));
outLoadType = out(:, 2);
propNBC.loadType = cell2mat(outLoadType{1});
outFctHandle = out(:, 3);
propNBC.fctHandle = cell2mat(outFctHandle{1});

%% 11. Load the nodes that are candidates for contact
block = regexp(fstring, 'STRUCTURE_CONTACT_NODES', 'split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k}, '%f');
end
if ~isempty(out)
    out = out{1};
    propContact.nodeIDs = cell2mat(out(:, 1));
    propContact.numNodes = length(propContact.nodeIDs);
else
    propContact.nodeIDs = [];
    propContact.numNodes = 0;
end

%% 12. Load coupled structure nodes for FSI
block = regexp(fstring, 'STRUCTURE_FSI_NODES', 'split'); 
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

%% 13. Get edge connectivity arrays for the Neumann edges
if strcmp(outMsg, 'outputEnabled')
    fprintf('>> Neumann boundary edges: %d \n', length(propNBC.nodes) - 1);
end
propNBC.lines = zeros(length(unique(propNBC.nodes)) - 1, 3);
counterLines = 1;
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
propNBC.fctHandle = fctHandle;

%% 14. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('\nParsing took %.2d seconds \n\n', computationalTime);
    fprintf('_________________________Parsing Ended__________________________\n');
    fprintf('################################################################\n\n\n');
end

end
