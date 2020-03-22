function [displacement,lagrange] = solveSignoriniLagrange_1...
    (analysis,strMsh,homDOFs,inhomDOFs,valuesInhomDOFs,NBC,bodyForces,parameters,...
    contactSegments,computeStiffMtxLoadVct,solve_LinearSystem,...
    propNLinearAnalysis,propContact,gaussInt,caseName,pathToOutput,...
    isUnitTest,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%                  Fabien Pean
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
%        propContact : structure containing the global numbering of the
%                      canditate contact nodes
%                  F : Global load vector
%    contactSegments : Containts the coordinates of the vertices of the 
%                   straight segments which form the rigid body's boundary:
%                    .numSegments : Number of straight segments
%                         .points : Array containing the coordinates of the 
%                                   vertices of each line segment X0 - X1 
%                                   in the form [x0, y0 ; x1 , y1]
%                        .normals : Array of the coordinates of each 
%                                   outward unit normal vector
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
% 1. Remove fully constrained nodes and Compute the gap function
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

%% 0. Read input

% Number of nodes in the mesh
noNodes = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
noDOFs = 2*noNodes;

% GLobal DOF numbering
DOFNumbering = 1:noDOFs;

% Assign dummy variables
uSaved = 'undefined';
uDot = 'undefined';
uDotSaved = 'undefined';
propStrDynamics = 'undefined';
loadFactor = 'undefined';

% Prescribed DOFs
prescribedDOFs = sort(horzcat(homDOFs, inhomDOFs));

% Tabulation for the output in the command window
tab = '';

% Initialize output array
dHat = zeros(noDOFs,1);

% Title for the output file
title = 'geometrically linear steady-state plane stress analysis';

% Steady-state analysis
t = 0;

%% 1. Remove fully constrained nodes and Compute the gap function
fullyConstrainedNodes = false(propContact.numberOfNodes,1);
for i = 1:length(propContact.nodeIDs)
    DOFs = 2*propContact.nodeIDs(i,1)-1 : 2*propContact.nodeIDs(i,1);
    if ismember(DOFs(1),prescribedDOFs) && ismember(DOFs(2),prescribedDOFs)
        fullyConstrainedNodes(i,1) = true; 
    end
end
propContact.nodeIDs(fullyConstrainedNodes) = [];
propContact.numberOfNodes = length(propContact.nodeIDs);
% x = mesh.nodes(propContact.nodeIDs(:),1);
% y = mesh.nodes(propContact.nodeIDs(:),2);
% plot(x,y,'ro');

%% Compute initial gap function
propContact = computeGapFunction(strMsh,propContact,contactSegments);

%% Compute force vector
F = computeLoadVctFEMPlateInMembraneAction...
    (strMsh,NBC,t,gaussInt,outMsg);

%% 2. Compute the master stiffness matrix of the structure
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Computing master stiffness matrix... \n');
end
K = computeStiffMtxLoadVct(analysis,dHat,uSaved,uDot,uDotSaved, ...
    DOFNumbering,strMsh,F,loadFactor,bodyForces,propStrDynamics,...
    parameters,gaussInt);

%% 3. Create the expanded system of equations
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Creating the expanded system of equations... \n');
end

% Assemble the values of the normal vector of segments to the constraint matrix
C = buildConstraintMatrix(noDOFs,propContact,[],contactSegments);

% Build an expanded system of equations
K_exp = [K  C
         C' zeros(size(C,2))]; 

% Build a gap function
gapFunction = reshape(propContact.gap,[],1);

% Expand the F vector with the gap funciton
F_exp = [F
         - gapFunction];

% Get number of DOFs in expanded system
noDOFsFull = length(F_exp);

freeDOFs = 1:noDOFsFull;

clear C;
clear gapFunction;

% Initial values for the iteration:
displacement_exp = zeros(noDOFsFull,1);
it = 1;
inactive_nodes = [];

% Inialize convergence conditions
isCnd_DOFs = false;
isCnd_lagrange = false;

%% 4. Solve the system iteratively

% Repeat loop until all Lagrange multipliers are valid and system converges 
% to stable solution. If the Lagrange multiplier of the current node are 
% valid and we had no change in the set of inactive_nodes during the last
% iteration the solution has been found and the loop is terminated
while (it <= propContact.maxIter && ~(isCnd_DOFs && isCnd_lagrange))
 
    % Assign inactive DOFs to the ones from previous iteration
    inactive_old_nodes = inactive_nodes;
   
    %% 4.1 Determine inactive nodes

    % Detect non-penetrating nodes and nodes with non-compressive Lagrange multipliers
    inactive_nodes = detectInactiveNodes(noDOFs,strMsh,propContact,displacement_exp,contactSegments);

    %% 4.2 Reduce the system of equations according to the constraints

    % Merge the homDBC and inactive_DOFs into a vector of DOFs that are no
    % longer needed. Equations with this numbers will be deleted due to a
    % dirichlet boundary condition or a contact constraint
    constrained_DOFs = unique([homDOFs inhomDOFs inactive_nodes]);
    
    freeDOFs_iterate = freeDOFs;
    freeDOFs_iterate(ismember(freeDOFs_iterate,constrained_DOFs)) = [];

    %% 4.3 Solve the reduced system of equations
    
    if norm(valuesInhomDOFs) ~= 0
        F_exp = F_exp - K_exp(:,inhomDOFs)*valuesInhomDOFs';
    end
    
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(tab,'>> Solving the linear system of %d equations\n'),length(freeDOFs_iterate));
    end
    [displacement_red,hasLinearSystemConverged] = ...
        solve_LinearSystem(K_exp(freeDOFs_iterate,freeDOFs_iterate),F_exp(freeDOFs_iterate),zeros(length(freeDOFs_iterate),1));
    if ~hasLinearSystemConverged
        error('Linear equation solver has not converged');
    end

    %% 4.4 Assemble the expanded displacement/Lagrange multiplier vector
    displacement_exp(freeDOFs_iterate) = displacement_red;
    displacement_exp(constrained_DOFs) = 0;
    displacement_exp(inhomDOFs) = valuesInhomDOFs;
           
    %% 4.5 Evaluate convergence conditions
    
    % check if the non-active degrees of freedom did not change
    isCnd_DOFs = isequal(inactive_old_nodes,inactive_nodes);
    
    % check if lagrange multipliers are negative
    isCnd_lagrange = max(displacement_exp(noDOFs+1:length(F_exp))) <= 0;
    
    % update the iteration counter
    it = it + 1;

end % end while loop

%% 5. Get the values for the displacement and the Lagrange multipliers

% The first entries in displacement_exp correspond to the displacements
displacement = displacement_exp(1:noDOFs);

% The last entries in displacement_exp correspond to Lagrange multipliers
lagrangeMultipliers = displacement_exp(noDOFs+1 : noDOFsFull);

% Create an 1D vector of all contact nodes
allContactNodes = repmat(propContact.nodeIDs,contactSegments.number,1);

% Find the active nodes where lagrange multipliers apply
lagrange.active_nodes = allContactNodes(lagrangeMultipliers < 0);

% Keep only lagrange multipliers of the active nodes
lagrange.multipliers = lagrangeMultipliers(lagrangeMultipliers < 0);


%% 6. Print info
if strcmp(outMsg,'outputEnabled')
    % energy of the structure
    energy = displacement'*K*displacement;
    % output
    fprintf('\n');
    fprintf('Output information...\n');
    fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
    fprintf('\t #DOF: %d \n\t Energy norm of the structure: %4.2f\n',noDOFs,energy);
    if (it >= propContact.maxIter)
        fprintf('\t Max number of iterations of has been reached!\n');
    end
end

end