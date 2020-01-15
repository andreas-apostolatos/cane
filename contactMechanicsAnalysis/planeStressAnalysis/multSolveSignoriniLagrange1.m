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
function [displacement,lagrange] = multSolveSignoriniLagrange1...
    (mesh,homDBC,contactNodes,F,segmentPoints,materialProperties,analysis,maxIteration)
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

% Compute normal, parallel and position vectors of the segments:
segments = buildSegmentsData(segmentPoints);

% Compute for all nodes the specific gap and save it in the field 'cn.gap':
contactNodes = multComputeGapFunc(contactNodes, segments, segmentPoints);

%% 2. Compute the master stiffness matrix of the structure
fprintf('\t Computing master stiffness matrix... \n');

% Master stiffness matrix
K = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,materialProperties,analysis);

%% 3. Create the expanded system of equations
fprintf('\t Creating the expanded system of equations... \n');

% Assemble the values of the normal vector of segments to the constraint matrix:
C = multBuildConstraintMatrix(length(F),contactNodes,[],segments);

% Create a zero matrix for the bellow right side of the equation system:
N_active_node=0;
for j=1:size(contactNodes,2)
    N_active_node=N_active_node+size(contactNodes(j).indices,1);       
end
zero_matrix=zeros(N_active_node,N_active_node);

% Expand the C matrix with the zero matrix:
Ctmp=[C;zero_matrix];

% Expand the K matrix with the C matrix below:
K_expanded=[K;C'];

% Expand the K matrix with the Ctmp matrix on the right side:
K_expanded=[K_expanded,Ctmp];

% Expand the F vector with the gap constants for every node:
F_expanded = F;
for j=1:size(contactNodes,2)
    F_expanded = [F_expanded;squeeze(-contactNodes(j).gap(:,2))];
end

clear C;
clear Ctmp;
clear zero_matrix;

% Initial values for the itaration:
it=0;
d_expanded=zeros(length(F_expanded));
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
while (it==0 || (it<maxIteration && ~(isequal(inactive_old,inactive_nodes) && max(d_expanded(length(F)+1:length(F_expanded)))<=0)))
 
inactive_old=inactive_nodes;
% update iteration counter
it=it+1;
    
%% 4.1 Determine the inactive_nodes nodes

% Detect non-penetrating nodes and nodes with non-compressive Lagr. multipliers
inactive_nodes = multDetectInactiveNodes( length(F),contactNodes, d_expanded, segments );

%% 4.2 Reduce the system of equations according to the constraints

% Merge the vectors with equation numbers which will be deleted due to a
% Dirichlet boundary condition or a contact constraint:
homDBC_cb=[homDBC,inactive_nodes];
homDBC_cb=unique(homDBC_cb);
homDBC_cb=sort(homDBC_cb);


K_reduced = K_expanded;
F_reduced = F_expanded;

% Remove constrained DOFs and inactive_nodes equations
for i = length(homDBC_cb):-1:1
    K_reduced(:,homDBC_cb(i)) = [];
    K_reduced(homDBC_cb(i),:) = [];
    F_reduced(homDBC_cb(i)) = [];
end
 
%% 4.3 Solve the reduced system of equations
equations_counter=equations_counter+length(K_reduced);
fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(K_reduced), cond(K_reduced));

d_reduced = K_reduced\F_reduced;


%% 4.4 Assemble to the expanded displacement/Lagrange multiplier vector

% Build dexp out of dred:
d_expanded = buildFullDisplacement(length(F_expanded),homDBC_cb,d_reduced);

end %end while loop

%% 5. Get the values for the displacement and the Lagrange multipliers

% Select and save node numbers of active nodes
allnodes=[];
for j=1:size(contactNodes,2)
    allnodes=[allnodes,contactNodes(j).indices];
end
lagrange.active_nodes=setdiff(allnodes,allnodes(inactive_nodes-length(F)));

% The first entries correspond to the displacement
displacement=d_expanded(1:length(F));

% The last entries correspond to the Lagrange multipliers
lagrange.multipliers=d_expanded(length(F)+1:length(F_expanded));

% Keep only lagrange multipliers of the active nodes
lagrange.multipliers=lagrange.multipliers(lagrange.multipliers<0);


%% 5. Print info

fprintf('\n');
fprintf('Output informations...\n');
fprintf('\t Constraints solved in %d iterations.A total of %d equations were solved. \n',it,equations_counter);
fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
energy=displacement'*K*displacement;
fprintf('\t #DOF: %d .Energy norm of the structure : %4.2f\n',length(F),energy);

end