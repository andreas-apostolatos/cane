function [dHat,lambdaHat,nodeIDs_active] = solveSignoriniLagrange_1...
    (analysis,strMsh,homDOFs,inhomDOFs,valuesInhomDOFs,NBC,bodyForces,parameters,...
    contactSegments,computeStiffMtxLoadVct,solve_LinearSystem,...
    propNLinearAnalysis,propContact,propGaussInt,caseName,pathToOutput,...
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

% Number of DOFs in the finite element mesh
noDOFs = 2*noNodes;

% GLobal DOF numbering
DOFNumbering = 1:noDOFs;

% Assign dummy variables
uSaved = 'undefined';
uDot = 'undefined';
uDDot = 'undefined';
uDotSaved = 'undefined';
uDDotSaved = 'undefined';
uMeshALE = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
forceVct = 'undefined';
loadFactor = 'undefined';
propStrDynamics = 'undefined';
computeProblemMatricesSteadyState = 'undefined';

% Steady-state analysis
t = 0;

% Prescribed DOFs
prescribedDOFs = sort(horzcat(homDOFs, inhomDOFs));

% Tabulation for the output in the command window
tab = '\t';

% Initialize output array
dHat_stiffMtx = zeros(noDOFs,1);

% Initialize counter of contact iterations
counterContact = 1;

% Initialize array of active nodes
active_nodes = [];

% Inialize convergence conditions
isCnd_DOFs = false;
isCnd_LM = false;

% Title for the output file
title = 'geometrically linear steady-state plane stress analysis';

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

%% 2. Get number of Lagrange Multipliers DOFs, total number of DOFs, total DOF numbering and initialize solution vector
noDOFsLM = propContact.numberOfNodes*contactSegments.number;
noDOFsTotal = noDOFs + noDOFsLM;
freeDOFs = 1:noDOFsTotal;
dHat_stiffMtxLM = zeros(noDOFsTotal,1);

%% 3. Compute initial gap function
propContact = computeGapFunction(strMsh,propContact,contactSegments);
gapFunction = reshape(propContact.gap,[],1);

%% 4. Compute external force vector
F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,propGaussInt,outMsg);

%% 5. Compute the master stiffness matrix of the structure
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'>> Computing the stiffness matrix of the system\n'));
end
K = computeStiffMtxLoadVct(analysis,dHat_stiffMtx,uSaved,uDot,uDotSaved, ...
    DOFNumbering,strMsh,F,loadFactor,bodyForces,propStrDynamics,...
    parameters,propGaussInt);

%% 6. Create the expanded system of equations corresponding to the Lagrange Multipliers method
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Creating the expanded system of equations... \n');
end
C = buildConstraintMatrix(noDOFs, propContact, contactSegments);
stiffMtxLM = [K  C
              C' zeros(size(C,2))];
clear C;
resVecLM = [F
            - gapFunction];
clear gapFunction;

%% 7. Loop over all contact iterations
while counterContact <= propContact.maxIter && ~(isCnd_DOFs && isCnd_LM)
    %% 7i Assign inactive DOFs to the ones from previous contact iteration
    active_old_nodes = active_nodes;
   
    %% 7ii. Determine active contact nodes
    active_nodes = findActiveContactNodes2D(noDOFs,strMsh,propContact,dHat_stiffMtxLM,contactSegments);

    %% 7iii. Collect the DOFs of the homogeneous Dirichlet boundary conditions and the active contact nodes of the current contact iteration in one array
    homDOFs_iterate = horzcat(homDOFs, active_nodes);
    
    %% 7iv. Collect all constrained DOFs into one array
    constrained_DOFs = unique(horzcat(homDOFs_iterate, inhomDOFs));
    
    %% 7v. Find the free DOFs of the current iteration
    freeDOFs_iterate = freeDOFs;
    freeDOFs_iterate(ismember(freeDOFs_iterate,constrained_DOFs)) = [];
    
    %% 7vi. Solve the equation system of the current contact iteration
    [dHat_stiffMtxLM,~,~,~] = solve_FEMLinearSystem(analysis,uSaved,...
        uDotSaved,uDDotSaved,strMsh,forceVct,bodyForces,parameters,dHat_stiffMtxLM,uDot,uDDot,...
        massMtx,dampMtx,stiffMtxLM,resVecLM,computeProblemMatricesSteadyState,...
        DOFNumbering,freeDOFs_iterate,homDOFs_iterate,inhomDOFs,valuesInhomDOFs,uMeshALE,...
        solve_LinearSystem,propStrDynamics,propNLinearAnalysis,propGaussInt,strcat(tab,'\t'),outMsg);
           
    %% 7vii. Evaluate convergence conditions
    
    % Check whether the vector of active contact nodes has not changed
    isCnd_DOFs = isequal(active_old_nodes,active_nodes);
    
    % Check whether all Lagrange Multipliers DOFs are negative
    isCnd_LM = max(dHat_stiffMtxLM(noDOFs+1:length(resVecLM))) <= 0;
    
    %% 7viii. Update the iteration counter
    counterContact = counterContact + 1;
end

%% 8. Get the displacement solution
dHat = dHat_stiffMtxLM(1:noDOFs);

%% 9. Get the Lagrange Multipliers solution
lambdaHat = dHat_stiffMtxLM(noDOFs + 1:noDOFsTotal);

%% 10. Get the index of the active Lagrange Multipliers DOFs
allContactNodes = repmat(propContact.nodeIDs,contactSegments.number,1);
nodeIDs_active = allContactNodes(lambdaHat < 0);

%% 11. Return only the negative Lagrange Multipliers
lambdaHat = lambdaHat(lambdaHat < 0);

%% 12. Print info
% if strcmp(outMsg,'outputEnabled')
%     % energy of the structure
%     energy = dHat'*K*dHat;
%     % output
%     fprintf('\n');
%     fprintf('Output information...\n');
%     fprintf('\t %d active nodes found.\n',length(DOFs_LM.active_nodes));
%     fprintf('\t #DOF: %d \n\t Energy norm of the structure: %4.2f\n',noDOFs,energy);
%     if (counterContact >= propContact.maxIter)
%         fprintf('\t Max number of iterations of has been reached!\n');
%     end
% end

end