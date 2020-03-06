function [displacement,lagrange] = solveSignoriniLagrange_1...
    (mesh,homDBC,propContact,F,segments,materialProperties,analysis,maxIteration,outMsg)
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
%             outMsg : write 'outputEnabled' to show output during runtime
%
%             Output :
%       displacement : The resulting displacement field
%           lagrange : .multipliers  : values of the Lagrange multipliers
%                    : .active_nodes : node numbers of the active nodes
%
%% Function layout
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
%    4.4 Assemble to the expanded displacement/Lagrange  multiplier vector
%
%    4.5 Evaluate convergence conditions
% <-|
% 5. Get the values for the displacement and the Lagrange multipliers
%
% 6. Print info
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    if strcmp(analysis.type,'planeStress')
        fprintf('Plane stress analysis has been initiated \n');
    elseif strcmp(analysis.type,'planeStrain')
        fprintf('Plane strain analysis has been initiated \n');
    end
    fprintf('\n');
end    

%% 0. Remove fully constrained nodes

% Remove fully constrained nodes from the set of contact nodes
for n=size(propContact.nodeIDs,1):-1:1
    % Determine how many Dirichlet conditions correspond to the node:  
    nodeHasDirichlet = ismember(floor((homDBC+1)/2),propContact.nodeIDs(n));
    numberOfDirichlet = length(nodeHasDirichlet(nodeHasDirichlet==true));
    % If the 2D node has at least two Dirichlet conditions exclude it from the contact canditates  
    if (numberOfDirichlet>=2)
       propContact.nodeIDs(n)=[];
    end
end
propContact.numberOfNodes = length(propContact.nodeIDs);

%% 1. Compute the gap function

% Compute normal vectors of the segments
segments = buildSegmentsData(segments);

% Compute for all nodes the specific gap and save it in the field .gap
propContact = computeGapFunction(mesh,propContact,segments);

%% 2. Compute the master stiffness matrix of the structure
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Computing master stiffness matrix... \n');
end

% Get number of DOFs in original system
nDOFs = length(F);

% Master global stiffness matrix
K = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,materialProperties,analysis);

%% 3. Create the expanded system of equations
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Creating the expanded system of equations... \n');
end

% Assemble the values of the normal vector of segments to the constraint matrix
C = buildConstraintMatrix(nDOFs,propContact,[],segments);

% Build an expanded system of equations
K_exp = [K,C;C',zeros(size(C,2))]; 

% Build a gap function
gap = reshape(propContact.gap,[],1);

% Expand the F vector with the gap funciton
F_exp = [F;-gap];

% Get number of DOFs in expanded system
nDOFsFull = length(F_exp);

clear C;
clear gap;

% Initial values for the iteration:
displacement_exp = zeros(length(F_exp),1);
it = 1;
inactive_DOFs = [];

% Inialize convergence conditions
isCnd_DOFs = false;
isCnd_lagrange = false;

% Counts the number of equations which are solved during iteration
equations_counter = 0;

%% 4. Solve the system iteratively

% Repeat loop until all Lagrange multipliers are valid and system converges 
% to stable solution. If the Lagrange multiplier of the current node are 
% valid and we had no change in the set of inactive_nodes during the last
% iteration the solution has been found and the loop is terminated
while (it<maxIteration && ~(isCnd_DOFs && isCnd_lagrange))
 
    % Assign inactive DOFs to the ones from previous iteration
    inactive_old_DOFs = inactive_DOFs;
   
    %% 4.1 Determine inactive nodes

    % Detect non-penetrating nodes and nodes with non-compressive Lagrange multipliers
    inactive_DOFs = detectInactiveDOFs(nDOFs,mesh,propContact,displacement_exp,segments);

    %% 4.2 Reduce the system of equations according to the constraints

    % Merge the homDBC and inactive_DOFs into a vector of DOFs that are no
    % longer needed. Equations with this numbers will be deleted due to a
    % dirichlet boundary condition or a contact constraint
    constrained_DOFs = [homDBC,inactive_DOFs];
    constrained_DOFs = unique(constrained_DOFs);

    % initialize the reduced system
    K_red = K_exp;
    F_red = F_exp;

    % Remove constrained DOFs and inactive_DOFs equations
    K_red(:,constrained_DOFs) = [];
    K_red(constrained_DOFs,:) = [];
    F_red(constrained_DOFs) = [];

    %% 4.3 Solve the reduced system of equations
    
    % count the number of total solved equations
    equations_counter = equations_counter + length(K_red);
    
    if strcmp(outMsg,'outputEnabled')
        fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(K_red),cond(K_red));
    end
    
    % solve the reduced system using the backslash operator
    displacement_red = K_red\F_red;

    %% 4.4 Assemble to the expanded displacement/Lagrange multiplier vector

    % Build expanded displacement vector out of reduced displacement vector
    displacement_exp = buildFullDisplacement...
                       (nDOFsFull,constrained_DOFs,displacement_red);
           
    %% 4.5 Evaluate convergence conditions
    
    % check if the non-active degrees of freedom did not change
    isCnd_DOFs = isequal(inactive_old_DOFs,inactive_DOFs);
    
    % check if lagrange multipliers are negative
    isCnd_lagrange = max(displacement_exp(nDOFs+1:length(F_exp))) <= 0;
    
    % update iteration counter
    it = it + 1;

end % end while loop

%% 5. Get the values for the displacement and the Lagrange multipliers

% The first entries in displacement_exp correspond to the displacement
displacement = displacement_exp(1:nDOFs);

% The last entries in displacement_exp correspond to the Lagrange multipliers
lagrangeMultipliers = displacement_exp(nDOFs+1 : nDOFsFull);

% Create an 1D vector of all contact nodes
allContactNodes = repmat(propContact.nodeIDs,segments.number,1);

% Find the active nodes where lagrange multipliers apply
activeNodes = allContactNodes(lagrangeMultipliers < 0);
lagrange.active_nodes = activeNodes;

% Keep only lagrange multipliers of the active nodes
lagrange.multipliers = lagrangeMultipliers(lagrangeMultipliers < 0);


%% 6. Print info
if strcmp(outMsg,'outputEnabled')
    % energy of the structure
    energy = displacement'*K*displacement;
    % output
    fprintf('\n');
    fprintf('Output informations...\n');
    fprintf('\t Constraints solved in %d iterations. A total of %d equations were solved. \n',it,equations_counter);
    fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
    fprintf('\t #DOF: %d \n\t Energy norm of the structure: %4.2f\n',nDOFs,energy);
    if (it >= maxIteration)
        fprintf('\t Max number of iterations of has been reached!\n');
    end
end

end