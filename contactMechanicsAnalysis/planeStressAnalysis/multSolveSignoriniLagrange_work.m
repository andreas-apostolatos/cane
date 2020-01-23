%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %
%   Authors                                                               %
%   _______                                                               %
%                                                                         %
%   Fabien Pean, Andreas Hauso, Georgios Koroniotis                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [displacement,lagrange] = multSolveSignoriniLagrange_work...
    (mesh,homDBC,contactNodes,F,segments,materialProperties,analysis,maxIteration)
%% Function documentation
%
% Returns the displacement field and the Lagrange multipliers corresponding 
% to a plain stress/strain analysis for the given mesh of the geometry 
% together with its Dirichlet and Neumann boundary conditions and the 
% contact constraints for MULTIPLE rigids walls by applying the Lagrange
% multiplier method.
% In the structure array candidateNodes(j) can be specified which canditates 
% for contact nodes are related to a certain wall segment(j)
% 
%              Input :
%               mesh : Elements and nodes of the mesh
%             homDBC : Vector of the prescribed DoFs (by global numbering)
%     candidateNodes : STRUCTURE ARRAY '(j=1..n).indices' 
%                      containing the global numbering of the canditate-nodes 
%                      for contact to segment j [segmentPoints(:,:,j)] 
%                      in the field 'indices'
%                  F : Global load vector
%      segmentPoints : Matrix with the coordinates of two wall determining
%                      points, for every segment j=1..n
% materialProperties : The material properties of the structure
%
%             Output :
%       displacement : The resulting displacement field
%           lagrange : The resulting values of the Lagrange multipliers (*.multipliers)
%                      and the node numbers of the active nodes (*.active_nodes)
%
% Function layout :
%
% 0. Remove fully constrained nodes
%
% 1. Compute the gap function
%
% 2. Compute the master stiffness matrix of the structure
%
% 3. Create the expanded system of equations
%
% 4. Solve system iteratively
% |->
%    4.1 Determine the inactive_nodes nodes
%  
%    4.2 Reduce the system of equations according to the constraints
%  
%    4.3 Solve the reduced system of equations
%  
%    4.3 Assemble to the expanded displacement/Lagrange  multiplier vector
% <-|
% 5. Get the values for the displacement and the Lagrange multipliers
%
% 6.  Print info
%
%% Function main body
if strcmp(analysis.type,'planeStress')
    fprintf('Plane stress analysis has been initiated \n');
elseif strcmp(analysis.type,'planeStrain')
    fprintf('Plane strain analysis has been initiated \n');
end
fprintf('\n');

%% 0. Remove fully constrained nodes

for j=1:size(contactNodes,2)
    % Remove fully constrained nodes from the tests
    % from number of indices until 1
    for i=size(contactNodes(j).indices,1):-1:1
        % Determine how many Dirichlet conditions correspond to the node:  
        nodeHasDirichlet=ismember(floor((homDBC+1)/2),contactNodes(j).indices(i));
        numberOfDirichlet=length(nodeHasDirichlet(nodeHasDirichlet==1));
        % If the 2D node has at least two Dirichlet conditions exclude it from the contact canditates :  
        if (numberOfDirichlet>=2)
           contactNodes(j).indices(i)=[];
        end
    end
contactNodes(j).positions=mesh.nodes(contactNodes(j).indices,:);
end

%% 1. Compute the gap function

% Compute normal, parallel and position vectors of the segments
segments = buildSegmentsData(segments);

% Compute for all nodes the specific gap and save it in the field .gap
contactNodes = multComputeGapFunc(contactNodes, segments);

%% 2. Compute the master stiffness matrix of the structure
fprintf('\t Computing master stiffness matrix... \n');

% Get number of DOFs in original system
nDOF = length(F);
% Master global stiffness matrix
K = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,materialProperties,analysis);

%% 3. Create the expanded system of equations
fprintf('\t Creating the expanded system of equations... \n');

% Assemble the values of the normal vector of segments to the constraint matrix
C = multBuildConstraintMatrix(nDOF,contactNodes,[],segments);

% Create a zero matrix for the bellow right side of the equation system
zero_matrix = zeros(size(C,2),size(C,2));

% Build an expanded system of equations
K_exp=[K,C;C',zero_matrix]; 

% Expand the F vector with the gap funciton for every node:
F_exp = F;
for j=1:size(contactNodes,2)
    F_exp = [F_exp;-contactNodes(j).gap(:,2)];
end

clear C;
clear zero_matrix;

% Initial values for the itaration:
it=0;
displacement_exp=zeros(length(F_exp));
inactive_nodes=[];
inactive_old=[];

% Counts the number of equations which are solved during iteration
equations_counter=0;

%% 4. Solve the system iteratively

% Repeat loop until all Lagrange multipliers are valid and system converges 
% to stable solution:
% Check whether solution is found
% If the Lagrange multiplier of the current node are valid and we had 
% no change in the set of inactive_nodes during the last iteration 
% the solution has been found and the loop is terminated
while (it==0 || (it<maxIteration && ~(isequal(inactive_old,inactive_nodes) && max(displacement_exp(nDOF+1:length(F_exp)))<=0)))
 
inactive_old=inactive_nodes;
% update iteration counter
it=it+1;
    
%% 4.1 Determine inactive nodes

% Detect non-penetrating nodes and nodes with non-compressive Lagr. multipliers
inactive_nodes = multDetectInactiveNodes(nDOF,contactNodes,displacement_exp,segments);

%% 4.2 Reduce the system of equations according to the constraints

% Merge the homDBC and inactive_nodes into a vector of nodes that are no
% longer necessary 
% Equations with this numbers will be deleted due to a dirichlet 
% boundary condition or a contact constraint
unnecessaryNodes = [homDBC,inactive_nodes];
unnecessaryNodes = unique(unnecessaryNodes);

% initialize the reduced system
K_red = K_exp;
F_red = F_exp;

% Remove constrained DOFs and inactive_nodes equations
for i = length(unnecessaryNodes):-1:1
    K_red(:,unnecessaryNodes(i)) = [];
    K_red(unnecessaryNodes(i),:) = [];
    F_red(unnecessaryNodes(i)) = [];
end
 
%% 4.3 Solve the reduced system of equations
equations_counter = equations_counter + length(K_red);
fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(K_red), cond(K_red));

% solve using the backslash operator
displacement_red = K_red\F_red;

%% 4.4 Assemble to the expanded displacement/Lagrange multiplier vector

% Build dexp out of dred:
nDOFsFull = length(F_exp);
displacement_exp = buildFullDisplacement(nDOFsFull,unnecessaryNodes,displacement_red);

%% 4.5 Visualize the structure for each step

% Select and save node numbers of active nodes
allnodes=[];
for j=1:size(contactNodes,2)
    allnodes=[allnodes,contactNodes(j).indices];
end
lagrange.active_nodes=setdiff(allnodes,allnodes(inactive_nodes-nDOF));
% The first entries correspond to the displacement
displacement=displacement_exp(1:nDOF);
% The last entries correspond to the Lagrange multipliers
lagrange.multipliers=displacement_exp(nDOF+1:length(F_exp));
% Keep only lagrange multipliers of the active nodes
lagrange.multipliers=lagrange.multipliers(lagrange.multipliers < 0);


% On the graph
graph.index = it+1;
% On the geometry visualization
graph.visualization.geometry = 'current';

graph.index = plot_currentConfigurationFEMPlateInMembraneAction(mesh,homDBC,displacement,graph);
plot_segments(segments);
plot_lagrangeMultipliers(mesh,displacement,lagrange,''); 



end %end while loop


%% 5. Get the values for the displacement and the Lagrange multipliers

% Select and save node numbers of active nodes
allnodes=[];
for j=1:size(contactNodes,2)
    allnodes=[allnodes,contactNodes(j).indices];
end
lagrange.active_nodes=setdiff(allnodes,allnodes(inactive_nodes-nDOF));

% The first entries correspond to the displacement
displacement=displacement_exp(1:nDOF);

% The last entries correspond to the Lagrange multipliers
lagrange.multipliers=displacement_exp(nDOF+1:length(F_exp));

% Keep only lagrange multipliers of the active nodes
lagrange.multipliers=lagrange.multipliers(lagrange.multipliers < 0);


%% 6. Print info

fprintf('\n');
fprintf('Output informations...\n');
fprintf('\t Constraints solved in %d iterations.A total of %d equations were solved. \n',it,equations_counter);
fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
energy=displacement'*K*displacement;
fprintf('\t #DOF: %d .Energy norm of the structure : %4.2f\n',nDOF,energy);

end