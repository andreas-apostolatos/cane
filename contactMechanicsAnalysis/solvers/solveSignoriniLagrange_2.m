function [displacement,lagrange] = solveSignoriniLagrange_2...
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
% Function layout :
%
% 1. Remove fully constrained nodes and Compute the gap function
%
% 2. Compute the master stiffness matrix of the structure
%
% 3. Reduce the system according to the given constraints
%|->
% 4. Solve the system according to Lagrange multipliers
%   4.1 Assemble to the complete displacement vector
%   4.2 Detect active nodes
%   4.3 Rebuild system if new active nodes found
%   4.4 Relax the system until ONLY valid Lagrange multipliers are computed
%   |->
%   4.4.1 Compute the displacement and Lagrange multipliers
%   4.4.2 Detect and delete non-valid Lagrange multipliers and related rows
%   <-|
%<-|
% 5. compute the complete load vector and verify the results
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

%% 1. Remove fully constrained nodes and Compute the gap function

% Remove fully contrained nodes
propContact = removeFullyConstrainedNodes(homDBC, propContact);
% x = mesh.nodes(propContact.nodeIDs(:),1);
% y = mesh.nodes(propContact.nodeIDs(:),2);
% plot(x,y,'ro');

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

%% 3. Reduce the system according to the given constraints
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Reducing the system according to the constraints... \n');
end

K_red = K;
F_red = F;

% Remove constrained DOFs (homDBCs)
K_red(:,homDBC) = [];
K_red(homDBC,:) = [];
F_red(homDBC) = [];

%% 4. Solve the system according to Lagrange multipliers

% Initialize the displacement vector
displacement_red = zeros(length(F_red),1);

% Initial values for the iteration
isCndMain = true;
active_nodes = [];
numberOfActiveNodes = 0;
it = 1;

% Counts the number of equations which are solved during iteration   
equations_counter = 0;
    
% Iterate until no more invalid Lagrange multipliers AND no new active nodes
% are added in the pool AND max number of iterations in not reached
while(isCndMain && it <= maxIteration)    

    %% 4.1 Build the expanded (complete) displacement vector
    displacement_exp = buildFullDisplacement(nDOFs,homDBC,displacement_red);

    %% 4.2 Detect active nodes
    active_nodes_temp = detectActiveNodes(mesh,propContact,displacement_exp,segments);

    % Evaluate convergence condition
    if(isequaln(active_nodes_temp,active_nodes) && it~=1)
        isCndMain = false;
    end

    %% 4.3 Rebuild system if new active nodes found
    if(~isempty(active_nodes_temp) && isCndMain)

        % Update the active nodes
        active_nodes = active_nodes_temp;
        numberOfActiveNodes = length(active_nodes);

        % Build constraint matrix C and force vector F
        C  = buildConstraintMatrix(nDOFs,propContact,active_nodes,segments);
        F_red = buildRHS(F,propContact,active_nodes,segments);

        % Build master system matrix
        K_red = [K,C;C',zeros(size(C,2))];

        % Reduce the system according to the BCs
        K_red(:,homDBC) = [];
        K_red(homDBC,:) = [];
        F_red(homDBC) = [];
        
    end
    
    %% 4.4 Relax the system until ONLY valid Lagrange multipliers are computed
    
    % Initialize Lagrange condition
    isCndLagrange = true;
    
    while(isCndLagrange && it <= maxIteration)

        isCndLagrange = false;

        %% 4.4.1 Compute the displacement and Lagrange multipliers
        
        % Count the number of total solved equations
        equations_counter = equations_counter + length(K_red);
        
        if strcmp(outMsg,'outputEnabled')
            fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(K_red),cond(K_red));
        end
        
        % Solve using the backslash operator
        displacement_red = K_red\F_red;
        
        %% 4.4.2 Detect and delete non-valid Lagrange multipliers and related rows
        
        % Get the number of DOFs in reduced system
        nDOFs_red = length(displacement_red);
        
        % Lagrange multipliers indices
        lmDOFsIndices = nDOFs_red-numberOfActiveNodes+1 : nDOFs_red;
        
        % Check if all Lagrange multipliers are valid
        if max(displacement_red(lmDOFsIndices)) > 0
            isCndLagrange = true;
            isCndMain = true;
        end
        
        % Find the indices of only non-compressive Lagrange Multipliers
        lmDOFsIndices = lmDOFsIndices(displacement_red(lmDOFsIndices)>=0);
       
        % Delete non-valid Lagrange multipliers and related rows
        K_red(:,lmDOFsIndices) = [];
        K_red(lmDOFsIndices,:) = [];
        F_red(lmDOFsIndices) = [];
        
        % Delete active nodes related to non-valid Lagrange multipliers
        active_nodes(lmDOFsIndices - nDOFs_red+numberOfActiveNodes) = [];
        
        % Update the number of active nodes
        numberOfActiveNodes = length(active_nodes);
        
        % update the iteration counter
        it = it + 1;

    end % end of inner while loop

end % end of main while loop

%% 5. Get the values for the displacement and the Lagrange multipliers

% Build full displacement vector
displacement = buildFullDisplacement(nDOFs,homDBC,displacement_red);

% Create an 1D vector of all contact nodes
allContactNodes = repmat(propContact.nodeIDs,segments.number,1);

% Find the active nodes where lagrange multipliers apply
lagrange.active_nodes = allContactNodes(active_nodes);

% Keep only lagrange multipliers of the active nodes
lagrange.multipliers = ...
    displacement_red(nDOFs_red-numberOfActiveNodes+1 : nDOFs_red);

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
    if it >= maxIteration
        fprintf('\t Max number of iterations of has been reached !! Not Converged !!\n');
    end
end

end