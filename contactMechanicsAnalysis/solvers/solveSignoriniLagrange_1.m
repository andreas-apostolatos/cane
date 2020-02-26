function [displacement,lagrange] = solveSignoriniLagrange_1...
    (mesh,homDBC,contactNodes,F,segments,materialProperties,analysis,maxIteration)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%                  -------------------
%                  Fabien Pean
%                  Andreas Hauso
%                  Georgios Koroniotis
%
%% Function documentation
%
% Returns the displacement field and the Lagrange multipliers corresponding 
% to a plain stress/strain analysis for the given mesh and geometry 
% together with its Dirichlet and Neumann boundary conditions and the 
% contact constraints for multiple rigids walls by applying the Lagrange
% multiplier method.
% 
%              Input :
%               mesh : Elements and nodes of the mesh
%             homDBC : Vector of the prescribed DoFs (by global numbering)
%       contactNodes : structure containing the global numbering of the
%                      canditate contact nodes
%                  F : Global load vector
%           segments : Matrix with the coordinates of two wall determining
%                      points, for every segment j=1..n
% materialProperties : The material properties of the structure
%           analysis : Structure about the analysis type
%       maxIteration : Maximum number of iterations
%
%             Output :
%       displacement : The resulting displacement field
%           lagrange : .multipliers  : values of the Lagrange multipliers
%                    : .active_nodes : node numbers of the active nodes
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

% Remove fully constrained nodes from the tests
% from number of indices until 1
for i=size(contactNodes.indices,1):-1:1
    % Determine how many Dirichlet conditions correspond to the node:  
    nodeHasDirichlet=ismember(floor((homDBC+1)/2),contactNodes.indices(i));
    numberOfDirichlet=length(nodeHasDirichlet(nodeHasDirichlet==1));
    % If the 2D node has at least two Dirichlet conditions exclude it from the contact canditates  
    if (numberOfDirichlet>=2)
       contactNodes.indices(i)=[];
    end
end
contactNodes.positions = mesh.nodes(contactNodes.indices,:);

%% 1. Compute the gap function

% Compute normal, parallel and position vectors of the segments
segments = buildSegmentsData(segments);

% Compute for all nodes the specific gap and save it in the field .gap
contactNodes = computeGapFunction(contactNodes,segments);

%% 2. Compute the master stiffness matrix of the structure
fprintf('\t Computing master stiffness matrix... \n');

% Get number of DOFs in original system
nDOFs = length(F);
% Master global stiffness matrix
K = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,materialProperties,analysis);

%% 3. Create the expanded system of equations
fprintf('\t Creating the expanded system of equations... \n');

% Assemble the values of the normal vector of segments to the constraint matrix
C = buildConstraintMatrix(nDOFs,contactNodes,[],segments);

% Build an expanded system of equations
K_exp = [K,C;C',zeros(size(C,2))]; 

% Expand the F vector with the gap funciton for every node:
F_exp = F;
for j=1:segments.number
    F_exp = [F_exp;-contactNodes.gap(:,2,j)];
end

clear C;

% Initial values for the itaration:
it = 0;
displacement_exp = zeros(length(F_exp),1);
inactive_DOFs = [];
inactive_old_DOFs = [];

% Counts the number of equations which are solved during iteration
equations_counter = 0;

%% 4. Solve the system iteratively

% Repeat loop until all Lagrange multipliers are valid and system converges 
% to stable solution. If the Lagrange multiplier of the current node are 
% valid and we had no change in the set of inactive_nodes during the last
% iteration the solution has been found and the loop is terminated
while (it==0 || (it<maxIteration && ~(isequal(inactive_old_DOFs,inactive_DOFs) && max(displacement_exp(nDOFs+1:length(F_exp)))<=0)))
 
inactive_old_DOFs = inactive_DOFs;
% update iteration counter
it = it + 1;
    
%% 4.1 Determine inactive nodes

% Detect non-penetrating nodes and nodes with non-compressive Lagrange multipliers
inactive_DOFs = detectInactiveNodes(nDOFs,contactNodes,displacement_exp,segments);

%% 4.2 Reduce the system of equations according to the constraints

% Merge the homDBC and inactiveDOFs into a vector of DOFs that are no
% longer necessary. Equations with this numbers will be deleted due to a 
% dirichlet boundary condition or a contact constraint
unnecessaryDOFs = [homDBC,inactive_DOFs];
unnecessaryDOFs = unique(unnecessaryDOFs);

% initialize the reduced system
K_red = K_exp;
F_red = F_exp;

% Remove constrained DOFs and inactiveDOFs equations
K_red(:,unnecessaryDOFs) = [];
K_red(unnecessaryDOFs,:) = [];
F_red(unnecessaryDOFs) = [];

 
%% 4.3 Solve the reduced system of equations
equations_counter = equations_counter + length(K_red);
fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(K_red),cond(K_red));

% solve using the backslash operator
displacement_red = K_red\F_red;

%% 4.4 Assemble to the expanded displacement/Lagrange multiplier vector

% Build dexp out of dred:
nDOFsFull = length(F_exp);
displacement_exp = buildFullDisplacement(nDOFsFull,unnecessaryDOFs,displacement_red);

end % end while loop

%% 5. Get the values for the displacement and the Lagrange multipliers

% Select and save node numbers of active nodes
allContactNodes=[];
for j=1:segments.number
    allContactNodes = [allContactNodes;contactNodes.indices];
end

% The first entries in displacement_exp correspond to the displacement
displacement = displacement_exp(1:nDOFs);

% The last entries in displacement_exp correspond to the Lagrange multipliers
lagrangeMultipliers = displacement_exp(nDOFs+1 : nDOFsFull);

% Find the active nodes where lagrange multipliers apply
activeNodes = allContactNodes(lagrangeMultipliers < 0);
lagrange.active_nodes = activeNodes;

% Keep only lagrange multipliers of the active nodes
lagrange.multipliers = lagrangeMultipliers(lagrangeMultipliers < 0);


%% 6. Print info

% energy of the structure
energy = displacement'*K*displacement;

% output
fprintf('\n');
fprintf('Output informations...\n');
fprintf('\t Constraints solved in %d iterations. A total of %d equations were solved. \n',it,equations_counter);
fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
fprintf('\t #DOF: %d \n\t Energy norm of the structure: %4.2f\n',nDOFs,energy);
if it >= maxIteration
    fprintf('\t Max number of iterations of has been reached !! Not Converged !!\n');
end

end