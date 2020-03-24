function [dHat,lambdaHat,nodeIDs_active,FComplete,minElSize] = ...
    solveSignoriniLagrange_2...
    (analysis,strMsh,homDOFs,inhomDOFs,valuesInhomDOFs,NBC,bodyForces,...
    parameters,segmentsContact,computeStiffMtxLoadVct,solve_LinearSystem,...
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
% to a plate in membrane action analysis for the given mesh and geometry 
% together with its Dirichlet and Neumann boundary conditions and the 
% contact constraints for multiple rigids walls using applying the Lagrange
% multiplier method
% 
%                  Input :
%               analysis : .type : The analysis type
%                 strMsh : Nodes and elements in the mesh
%                homDOFs : The global numbering of the nodes where 
%                          homogeneous Dirichlet boundary conditions are 
%                          applied 
%              inhomDOFs : The global numbering of the nodes where 
%                          inhomogeneous Dirichlet boundary conditions are 
%                          applied
%        valuesInhomDOFs : Prescribed values on the nodes where 
%                          inhomogeneous Dirichlet boundary conditions are 
%                          applied
%                    NBC :    .nodes : The nodes where Neumann boundary 
%                                      conditions are applied
%                          .loadType : The type of the load for each 
%                                      Neumann node
%                         .fctHandle : The function handle for each Neumann 
%                                      node for the computation of the load 
%                                      vector (these functions are unde the 
%                                      folder load)
%             bodyForces : Function handle to body force vector computation
%             parameters : Problem specific technical parameters
%        segmentsContact : Data sturcture containing information about the 
%                          boundaries of the rigid wall :
%                           .points : a list of 2x2 matrices containing
%                                      end points of the segment(s)
%                            .number : total number of segments
%                           .normals : normal vector of each segment
% computeStiffMtxLoadVct : Function handle to the computation of the
%                          stiffness matrix and load vector
%     solve_LinearSystem : Function handle to the computation of the 
%                          solution to the linear equation system
%    propNLinearAnalysis :     .scheme : The employed nonlinear scheme
%                           .tolerance : The residual tolerance
%                             .maxIter : The maximum number of the 
%                                        nonlinear iterations
%            propContact : Data structure containing the contact properties
%                              .nodeIds : global numbering of contact nodes
%                        .numberOfNodes : number of nodes
%           propGaussInt : On the numerical integration (quadrature)
%                                  .type : 'default', 'user'
%                            .domainNoGP : Number of Gauss Points for the 
%                                          domain integration
%                          .boundaryNoGP : Number of Gauss Points for the
%                                          boundary integration
%               caseName : The name of the case in the inputGiD case folder
%             isUnitTest : Flag on whether the case is a unit test case
%           pathToOutput : Define the path to where to write out the 
%                          results
%                 outMsg : On outputting information
%
%                 Output :
%                   dHat : The nodal displacement field
%              lambdaHat : The valid (negative) Lagrange Multipliers
%         nodeIDs_active : The IDs of the active (contact) nodes
%              FComplete : The complete force vector
%              minElSize : The minimum element area size in the mesh
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
    fprintf('_____________________________________________________________\n');
    fprintf('#############################################################\n');
    fprintf('Compute the displacement field for a plate in membrane action\n');
    fprintf('problem subject to rigid contact boundary has been initiated\n');
    fprintf('_____________________________________________________________\n\n');

    % start measuring computational time
    tic;
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
noTimeStep = 1;

% Connectivity arrays for the DOFs into the resulting vectors
DOF4Output = [1:2:noDOFs-1
              2:2:noDOFs];

% Prescribed DOFs
prescribedDOFs = sort(horzcat(homDOFs, inhomDOFs));

% Tabulation for the output in the command window
tab = '\t';

% Initialize output array
dHat_stiffMtx = zeros(noDOFs,1);

% Initialize counter of contact iterations
counterContact = 1;

% Initialize array of active nodes
nodesActive = [];

% Inialize convergence conditions
isCnd_DOFs = false;
isCnd_LM = false;

% Title for the output file
title = 'Contact analysis for a plate in membrane action';

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

%% 3. Compute initial gap function
propContact = computeGapFunction(strMsh,propContact,segmentsContact);
gapFunction = reshape(propContact.gap,[],1);

%% 4. Compute external force vector
F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,propGaussInt,outMsg);

%% 5. Compute the master stiffness matrix of the structure
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'>> Computing the stiffness matrix of the system\n'));
end
% [K,F,minElSize] = computeStiffMtxLoadVct...
%     (analysis,dHat_stiffMtx,uSaved,uDot,uDotSaved,DOFNumbering,strMsh,F,...
%     loadFactor,bodyForces,propStrDynamics,parameters,propGaussInt);
K = computeStiffnessMatrixPlateInMembraneActionLinear...
    (strMsh,parameters,analysis);

%% 3. Reduce the system according to the given constraints
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Reducing the system according to the constraints... \n');
end

stiffMtxLM = K;
resVecLM = F;

freeDOFs = 1:noDOFs;
% Remove constrained DOFs (homDBCs)
freeDOFs_iterate = freeDOFs;
homDOFs_iterate = homDOFs;


% Remove constrained DOFs (homDBCs)
% stiffMtxLM(:,homDBC) = [];
% stiffMtxLM(homDBC,:) = [];
% resVecLM(homDBC) = [];

%% 4. Solve the system according to Lagrange multipliers

% Initialize the displacement vector
dHat_stiffMtxLM = zeros(noDOFs,1);

% Initial values for the iteration
isCndMain = true;
nodesActive = [];
numberOfActiveNodes = 0;
it = 1;

% Iterate until no more invalid Lagrange multipliers AND no new active nodes
% are added in the pool AND max number of iterations in not reached
while(isCndMain && it <= propContact.maxIter)    

    %% 4.1 Build the expanded (complete) displacement vector
    %displacement_exp = buildFullDisplacement(noDOFs,homDBC,dHat_stiffMtxLM);

    %% 4.2 Detect active nodes
    nodesActiveSaved = findActiveLagrangeMultipliersContact2D(strMsh,dHat_stiffMtxLM,segmentsContact,propContact);

    % Evaluate convergence condition
    if(isequaln(nodesActiveSaved,nodesActive) && it~=1)
        isCndMain = false;
    end

    %% 4.3 Rebuild system if new active nodes found
    if(~isempty(nodesActiveSaved) && isCndMain)

        % Update the active nodes
        nodesActive = nodesActiveSaved;
        numberOfActiveNodes = length(nodesActive);

        % Build constraint matrix C and force vector F
        C = buildConstraintMatrix(noDOFs, propContact,nodesActive, segmentsContact);
        resVecLM = buildRHS(F,propContact,nodesActive,segmentsContact);

        % Build master system matrix
        stiffMtxLM = [K  C
                      C' zeros(size(C,2))];
        clear C;

        % Reduce the system according to the BCs
        freeDOFs_iterate = 1:length(resVecLM);
        homDOFs_iterate = homDOFs;
        
        
        % Reduce the system according to the BCs
        %stiffMtxLM(:,homDBC) = [];
        %stiffMtxLM(homDBC,:) = [];
        %resVecLM(homDBC) = [];
        
    end
    
    %% 4.4 Relax the system until ONLY valid Lagrange multipliers are computed
    
    % Initialize Lagrange condition
    isCndLagrange = true;
    
    while(isCndLagrange && it <= propContact.maxIter)

        isCndLagrange = false;

        %% 4.4.1 Compute the displacement and Lagrange multipliers
        [dHat_stiffMtxLM,FComplete,~,~] = solve_FEMLinearSystem ...
            (analysis,uSaved,uDotSaved,uDDotSaved,strMsh,forceVct,bodyForces,...
            parameters,dHat_stiffMtxLM,uDot,uDDot,massMtx,dampMtx,stiffMtxLM,...
            resVecLM,computeProblemMatricesSteadyState,DOFNumbering,...
            freeDOFs_iterate,homDOFs_iterate,inhomDOFs,valuesInhomDOFs,...
            uMeshALE,solve_LinearSystem,propStrDynamics,propNLinearAnalysis,...
            propGaussInt,strcat(tab,'\t'),outMsg);
        
        %% 4.4.2 Detect and delete non-valid Lagrange multipliers and related rows
        
        % Get the number of DOFs in reduced system
        nDOFs_red = length(dHat_stiffMtxLM);
        
        % Lagrange multipliers indices
        lmDOFsIndices = nDOFs_red-numberOfActiveNodes+1 : nDOFs_red;
        
        % Check if all Lagrange multipliers are valid
        if max(dHat_stiffMtxLM(lmDOFsIndices)) > 0
            isCndLagrange = true;
            isCndMain = true;
        end
        
        % Find the indices of only non-compressive Lagrange Multipliers
        lmDOFsIndices = lmDOFsIndices(dHat_stiffMtxLM(lmDOFsIndices)>=0);
       
        
        
        homDOFs_iterate = horzcat(homDOFs, lmDOFsIndices);
        constrained_DOFs = unique(horzcat(homDOFs_iterate, inhomDOFs));
        
        freeDOFs_iterate = 1:length(resVecLM);
        freeDOFs_iterate(ismember(freeDOFs_iterate,constrained_DOFs)) = [];
        
        
        
        
        
        
        
        % Delete non-valid Lagrange multipliers and related rows
        %stiffMtxLM(:,lmDOFsIndices) = [];
        %stiffMtxLM(lmDOFsIndices,:) = [];
        %resVecLM(lmDOFsIndices) = [];
        
        % Delete active nodes related to non-valid Lagrange multipliers
        nodesActive(lmDOFsIndices - nDOFs_red+numberOfActiveNodes) = [];
        
        % Update the number of active nodes
        numberOfActiveNodes = length(nodesActive);
        
        % update the iteration counter
        it = it + 1;

    end % end of inner while loop

end % end of main while loop

%% 5. Get the values for the displacement and the Lagrange multipliers

% Build full displacement vector
dHat = dHat_stiffMtxLM(1:noDOFs);

% Create an 1D vector of all contact nodes
allContactNodes = repmat(propContact.nodeIDs,segments.number,1);

% Find the active nodes where lagrange multipliers apply
lagrange.active_nodes = allContactNodes(nodesActive);

% Keep only lagrange multipliers of the active nodes
lagrange.multipliers = ...
    dHat_stiffMtxLM(nDOFs_red-numberOfActiveNodes+1 : nDOFs_red);

%% 6. Print info
if strcmp(outMsg,'outputEnabled')
    % energy of the structure
    energy = displacement'*K*displacement;
    % output
    fprintf('\n');
    fprintf('Output informations...\n');
    fprintf('\t Constraints solved in %d iterations. A total of %d equations were solved. \n',it,equations_counter);
    fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
    fprintf('\t #DOF: %d \n\t Energy norm of the structure: %4.2f\n',noDOFs,energy);
    if it >= maxIteration
        fprintf('\t Max number of iterations of has been reached !! Not Converged !!\n');
    end
end

end