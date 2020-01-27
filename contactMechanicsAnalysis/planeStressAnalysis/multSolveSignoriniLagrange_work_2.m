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
function [displacement, lagrange] = multSolveSignoriniLagrange_work_2...
    (mesh,homDBC,contactNodes,F, segments,materialProperties,analysis,maxIteration)
%% Function documentation
%
% Returns the displacement field and the Lagrange multipliers corresponding 
% to a plain stress/strain analysis for the given mesh of the geometry 
% together with its Dirichlet and Neumann boundary conditions and the 
% contact constraints for MULTIPLE rigid walls by applying the Lagrange
% multiplier method.
% In the structure array cn(j) can be specified which canditates for 
% contact nodes are related to a certain wall segment (j)
% 
%              Input :
%               mesh : Elements and nodes of the mesh
%                 rb : Vector of the prescribed DoFs (by global numbering)
%                 cn : STRUCTURE ARRAY 'cn(j=1..n).indices' 
%                      containing the global numbering of the canditate-nodes 
%                      for contact to segment j [segmentPoints(j,:,:)] 
%                      in the field 'indices'
%                  F : Global load vector
%     segmentsPoints : Matrix with the coordinates of two wall determining
%                      points, for every segment j=1..n
% materialProperties : The material properties of the structure
%              graph : Structure containing information on the graphics
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
% 3. Reduce the system according to the given constraints
%|->
% 4. Solve the system according to Lagrange multipliers
%   4.1 Assemble to the complete displacement vector
%   4.2 Detect active nodes
%   4.3 Rebuild system if new active nodes found
%       4.3.1 Update number of active nodes
%       4.3.2 Build constraint matrix C and rhs F
%       4.3.3 Build master system matrix
%       4.3.4 Reduce the system according to the BCs
%   4.4 Relax the system until ONLY valid Lagrange multipliers are computed
%   |->
%   4.4.1 Compute the displacement & Lagrange multipliers
%   4.4.2 Detect and delete non-valid Lagrange multipliers and related rows
%   <-|
%<-|
% 5. compute the complete load vector and verify the results
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
    % If the 2D node has at least two Dirichlet conditions exclude it from the contact canditates :  
    if (numberOfDirichlet>=2)
       contactNodes.indices(i)=[];
    end
end
contactNodes.positions=mesh.nodes(contactNodes.indices,:);

% test the above function by drawing the nodes 
% x = contactNodes.positions(:,1);
% y = contactNodes.positions(:,2);
% plot(x,y,'ro');

%% 1. Compute the gap function

% Compute normal, parallel and position vectors of the segments
segments = buildSegmentsData(segments);

% Compute for all nodes the specific gap and save it in the field .gap
contactNodes = multComputeGapFunc(contactNodes,segments);

%% 2. Compute the master stiffness matrix of the structure
fprintf('\t Computing master stiffness matrix... \n');

% Get number of DOFs in original system
nDOFs = length(F);
% Master global stiffness matrix
K = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,materialProperties,analysis);

%% 3. Reduce the system according to the given constraints
fprintf('\t Reducing the system according to the constraints... \n');

K_red = K;
F_red = F;

for i = length(homDBC):-1:1
    K_red(:,homDBC(i)) = [];
    K_red(homDBC(i),:) = [];
    F_red(homDBC(i)) = [];
end

%% 4. Solve the system according to Lagrange multipliers

% Initialize the displacement vector
displacement_red = zeros(length(F_red),1);

% Initial values for the iteration
isCndMain = true;
activeNodes = [];
nActiveNodes = 0;
it = 0;

% Counts the number of equations which are solved during iteration   
equations_counter = 0;
    
% Iterate until no more invalid Lagrange multipliers AND no new active nodes
% are added in the pool
while(isCndMain && it<maxIteration)    

    %% 4.1 Assemble to the complete displacement vector
    displacement = buildFullDisplacement(nDOFs,homDBC,displacement_red);

    %% 4.2 Detect active nodes
    activeNodes_tmp = multDetectActiveNodes(contactNodes,displacement,segments);

    if(isequaln(activeNodes_tmp,activeNodes) && it~=0)
        isCndMain = false;
    end

    %% 4.3 Rebuild system if new active nodes found
    if(~isempty(activeNodes_tmp) && isCndMain)

        % Update number of active nodes
        activeNodes = activeNodes_tmp;
        nActiveNodes = length(activeNodes);

        % Build constraint matrix C and rhs F
        C  = multBuildConstraintMatrix(nDOFs,contactNodes,activeNodes,segments);
        F_red = multBuildRHS(F,contactNodes,activeNodes,segments);

        % Build master system matrix
        K_red = [K  C
                 C' zeros(size(C,2))];

        % Reduce the system according to the BCs
        K_red(:,homDBC) = [];
        K_red(homDBC,:) = [];
        F_red(homDBC) = [];
    end
    %% 4.4 Relax the system until ONLY valid Lagrange multipliers are computed
    isCndLagrange = true;
    
    while(isCndLagrange && it<maxIteration)

        it = it + 1;
        isCndLagrange = false;

        % Compute the displacement & Lagrange multipliers
        equations_counter=equations_counter+length(K_red);
        fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(K_red),cond(K_red));

        % Solve using the backslash operator
        displacement_red = K_red\F_red;

        
%         lmDOFsIndices = length(displacement_red)-nActiveNodes+1:length(displacement_red);
%         
%         
%         if norm(displacement_red(lmDOFsIndices) >= 0) ~= 0
%             isCndLagrange = true;
%             isCndMain = true;
%         end
%         
%         K_red(:,displacement_red(lmDOFsIndices) >= 0 )
        
        
        % Loop through the lagrange part of displacement vector
        for i=length(displacement_red) :-1: (length(displacement_red)-nActiveNodes+1)

            % Check if lagrange multiplier is non-compressive
            if(displacement_red(i)>=0)

                % Delete non-valid Lagrange multipliers and related rows
                K_red(:,i) = [];
                K_red(i,:) = [];
                F_red(i) = [];
                activeNodes(i-length(displacement_red)+nActiveNodes)=[];

                % Set conditions to true
                isCndLagrange = true;
                isCndMain = true;
            end
        end
        nActiveNodes=length(activeNodes);



        %% 4.5 Visualize the structure for each step

        % Select and save node numbers of active nodes
        allContactNodes=[];
        for j=1:segments.number
            allContactNodes=[allContactNodes;contactNodes.indices];
        end  

        % Build full displacement vector
        displacement = buildFullDisplacement(nDOFs,homDBC,displacement_red);

        % Keep only lagrange multipliers of the active nodes
        lagrange.multipliers = displacement_red(length(displacement_red)-nActiveNodes+1:length(displacement_red));
        lagrange.active_nodes = allContactNodes(activeNodes);

        % On the graph
        graph.index = it + 1;
        % On the geometry visualization
        graph.visualization.geometry = 'current';

        graph.index = plot_currentConfigurationFEMPlateInMembraneAction(mesh,homDBC,displacement,graph);
        plot_segments(segments);
        plot_lagrangeMultipliers(mesh,displacement,lagrange); 

    end % end of innner while loop

end % end of main while loop

%% 5. Compute the complete load vector and verify the results

% Select and save node numbers of active nodes
allContactNodes=[];
for j=1:segments.number
    allContactNodes=[allContactNodes;contactNodes.indices];
end  

% Build full displacement vector
displacement = buildFullDisplacement(nDOFs,homDBC,displacement_red);

% Keep only lagrange multipliers of the active nodes
lagrange.multipliers = displacement_red(length(displacement_red)-nActiveNodes+1:length(displacement_red));
lagrange.active_nodes = allContactNodes(activeNodes);

%% 6. Print info

% energy of the structure
energy = displacement'*K*displacement;

% output
fprintf('\n');
fprintf('Output informations...\n');
fprintf('\t Constraints solved in %d iterations. A total of %d equations were solved. \n',it,equations_counter);
fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
fprintf('\t #DOF: %d \n\t Energy norm of the structure: %4.2f\n',nDOFs,energy);

end