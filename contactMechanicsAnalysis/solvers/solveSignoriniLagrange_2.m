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
%% Function layout
%
% 0. Read input
%
% 1. Remove fully constrained nodes 
%
% 2. Compute initial gap function
%
% 3. Compute external force vector
%
% 4. Compute the master stiffness matrix of the structure
%
% 5. Initialize the system
%
% 6. Reduce the initial system according to the given constraints
%
% 7. Loop over all contact iterations
% ->
%    7i. Assign active DOFs to the ones from previous contact iteration
%
%   7ii. Determine active contact nodes
%
%  7iii. Evaluate main convergence condition - compare active nodes
%
%   7iv. Rebuild the expanded system if new active nodes are found
%
%    7v. Relax the system until ONLY valid Lagrange multipliers are computed
%       ->
%       7v.1 Compute the displacement and Lagrange multipliers
%
%       7v.2 Detect and delete non-valid Lagrange multipliers and DOFs
%
%       7v.3. Update the iteration counter
%       <-
% <-
%
% 8. Get the displacement solution
%
% 9. Write out the results into a file if the case is not a unit test
%
% 10. Get the Lagrange Multipliers solution
%
% 11. Get the index of the active Lagrange Multipliers DOFs
%
% 12. Appendix
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

% Initialize output array
dHat_stiffMtx = zeros(noDOFs,1);

% Tabulation for the output in the command window
tab = '\t';

% Initialize counter of contact iterations
counterContact = 1;

% Initialize array of active nodes
nodesActive = [];

% Initialize the displacement vector
dHat_stiffMtxLM = zeros(noDOFs,1);

% Inialize convergence conditions
isCnd_main = true;

% Title for the output file
title = 'Contact analysis for a plate in membrane action';

%% 1. Remove fully constrained nodes 
fullyConstrainedNodes = false(propContact.numberOfNodes,1);
for i = 1:length(propContact.nodeIDs)
    DOFs = 2*propContact.nodeIDs(i,1)-1 : 2*propContact.nodeIDs(i,1);
    if ismember(DOFs(1),prescribedDOFs) && ismember(DOFs(2),prescribedDOFs)
        fullyConstrainedNodes(i,1) = true; 
    end
end
propContact.nodeIDs(fullyConstrainedNodes) = [];
propContact.numberOfNodes = length(propContact.nodeIDs);

%% 2. Compute initial gap function
propContact = computeGapFunction(strMsh,propContact,segmentsContact);

%% 3. Compute external force vector
F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,propGaussInt,outMsg);

%% 4. Compute the master stiffness matrix of the structure
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'>> Computing the stiffness matrix of the system\n'));
end
% [K,F,minElSize] = computeStiffMtxLoadVct...
%     (analysis,dHat_stiffMtx,uSaved,uDot,uDotSaved,DOFNumbering,strMsh,F,...
%     loadFactor,bodyForces,propStrDynamics,parameters,propGaussInt);
K = computeStiffnessMatrixPlateInMembraneActionLinear...
    (strMsh,parameters,analysis);

%% 5. Initialize the system

% Assign intial values for the stiffness matrix K and force vector F
stiffMtxLM = K;
resVecLM = F;

% Total DOFs in initial expanded system
noDOFsTotal = length(resVecLM);

%% 6. Reduce the initial system according to the given constraints
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'>> Reducing the system according to the constraints\n'));
end

% Remove the prescribed DOFs from the set of free DOFs
homDOFs_iterate = homDOFs;

% Find the free DOFs of the current iteration
freeDOFs_iterate = 1:noDOFs;
freeDOFs_iterate(ismember(freeDOFs_iterate,prescribedDOFs)) = [];

%% 7. Loop over all contact iterations
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'>> Loop over all contact iterations\n'));
end
while(isCnd_main && counterContact <= propContact.maxIter) 
    %% 7i. Assign active DOFs to the ones from previous contact iteration
    nodesActiveSaved = nodesActive;
    
    %% 7ii. Determine active contact nodes
    nodesActive = findActiveLagrangeMultipliersContact2D(strMsh,dHat_stiffMtxLM,segmentsContact,propContact);

    %% 7iii. Evaluate main convergence condition - compare active nodes
    if(isequaln(nodesActive,nodesActiveSaved) && counterContact~=1)
        isCnd_main = false;
    end

    %% 7iv. Rebuild the expanded system if new active nodes are found
    if(~isempty(nodesActive) && isCnd_main)
        
        % Create the Constraint matrix C and resulting force vector
        C = buildConstraintMatrix(noDOFs,propContact,nodesActive,segmentsContact);
        resVecLM = buildRightHandSide(F,propContact,nodesActive,segmentsContact);
        
        % Create the expanded system of equations corresponding
        stiffMtxLM = [K  C
                      C' zeros(size(C,2))];
        clear C;
        
        % Get total number of DOFs in expanded system
        noDOFsTotal = length(resVecLM);
        
        % Initialize the displacement vector
        dHat_stiffMtxLM = zeros(noDOFsTotal,1);

        % Assign homDOFs to current iteration
        homDOFs_iterate = homDOFs;
        
        % Find the free DOFs of the current iteration
        freeDOFs_iterate = 1:noDOFsTotal;
        freeDOFs_iterate(ismember(freeDOFs_iterate,prescribedDOFs)) = [];
    end
    
    %% 7v. Relax the system until ONLY valid Lagrange multipliers are computed
    
    % Initialize Lagrange condition
    isCnd_lagrange = true;
    
    while(isCnd_lagrange && counterContact <= propContact.maxIter)

        isCnd_lagrange = false;
        
        %% 7v.1 Compute the displacement and Lagrange multipliers
        [dHat_stiffMtxLM,FComplete,~,~] = solve_FEMLinearSystem ...
            (analysis,uSaved,uDotSaved,uDDotSaved,strMsh,forceVct,bodyForces,...
            parameters,dHat_stiffMtxLM,uDot,uDDot,massMtx,dampMtx,stiffMtxLM,...
            resVecLM,computeProblemMatricesSteadyState,DOFNumbering,...
            freeDOFs_iterate,homDOFs_iterate,inhomDOFs,valuesInhomDOFs,...
            uMeshALE,solve_LinearSystem,propStrDynamics,propNLinearAnalysis,...
            propGaussInt,strcat(tab,'\t'),outMsg);
        
        %% 7v.2 Detect and delete non-valid Lagrange multipliers and DOFs
        if(~isempty(nodesActive))
            
            % Lagrange multipliers indices
            lmDOFsIndices = noDOFs+1 : noDOFsTotal;
            
            % Evaluate inner convergence condition - compare Lagrange Multipliers
            if max(dHat_stiffMtxLM(lmDOFsIndices)) > 0
                isCnd_lagrange = true;
                isCnd_main = true;
            end
            
            % Find the indices of only non-compressive Lagrange Multipliers
            lmDOFsIndices = lmDOFsIndices(dHat_stiffMtxLM(lmDOFsIndices) >= 0);

            % Collect the DOFs of the homDBC conditions and the active 
            % lagrange multipliers in one array
            homDOFs_iterate = horzcat(homDOFs, lmDOFsIndices);
            
            % Collect all constrained DOFs into one array
            constrained_DOFs = unique(horzcat(homDOFs_iterate, inhomDOFs));

            % Find the free DOFs of the current iteration
            freeDOFs_iterate = 1:noDOFsTotal;
            freeDOFs_iterate(ismember(freeDOFs_iterate,constrained_DOFs)) = [];
            
        end
        
        %% 7v.3. Update the iteration counter
        counterContact = counterContact + 1;
        
    end % End of inner while loop

end % End of main while loop

if (counterContact >= propContact.maxIter)
    warning('\t Max number of iterations has been reached!\n');
end

%% 8. Get the displacement solution
dHat = dHat_stiffMtxLM(1:noDOFs);

%% 9. Write out the results into a file if the case is not a unit test
if ~isUnitTest
    if strcmp(outMsg,'outputEnabled')
        fprintf('>> Writting out the results to "%s"\n',strcat(pathToOutput,caseName,'/'));
    end
    writeOutputFEMPlateInMembraneActionToVTK(analysis,propNLinearAnalysis,...
        propStrDynamics,strMsh,parameters,dHat,uDot,uDDot,DOF4Output,caseName,...
        pathToOutput,title,noTimeStep);
end

%% 10. Get the Lagrange Multipliers solution
lambdaHat = dHat_stiffMtxLM(noDOFs+1 : noDOFsTotal);

%% 11. Get the index of the active Lagrange Multipliers DOFs
allContactNodes = repmat(propContact.nodeIDs,segmentsContact.number,1);
nodeIDs_active = allContactNodes(nodesActive);

%% 12. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nContact analysis took %.2d seconds \n',computationalTime);
    fprintf('Number of active (contact) nodes %d\n\n',length(nodeIDs_active));
    fprintf('_____________________Linear Analysis Ended___________________\n');
    fprintf('#############################################################\n\n\n');
end

end