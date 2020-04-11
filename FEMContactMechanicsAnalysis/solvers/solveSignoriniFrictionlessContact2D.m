function [dHat, lambdaHat, nodeIDs_active, numInter, FComplete, minElSize] = ...
    solveSignoriniFrictionlessContact2D...
    (propAnalysis, strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, ...
    computeBodyForces, parameters, segmentsContact, computeStiffMtxLoadVct, ...
    solve_LinearSystem, propNLinearAnalysis, propContact, propGaussInt, ...
    caseName, pathToOutput, isUnitTest, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Fabien Pean
%                  Marko Leskovar
%
%% Function documentation
%
% Returns the displacement field and the Lagrange multipliers corresponding 
% to the contact analysis of a plate in membrane action without friction
% provided contact constraints for multiple rigids walls using an iterative 
% approach for solving the sequential quadratic programming problem (SQP).
% 
%                  Input :
%           propAnalysis : Structure containing general information on the 
%                          analysis,
%                              .type : The analysis type
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
%                propNBC : Structure containing information on the Neumann
%                          boundary conditions,
%                             .nodes : The nodes where Neumann boundary 
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
%               numInter : Number of contact iterations
%              FComplete : The complete force vector
%              minElSize : The minimum element area size in the mesh
%
% Function layout
%
% 0. Read input
%
% 1. Remove fully constrained nodes
%
% 2. Get number of Lagrange Multipliers DOFs, total number of DOFs, total DOF numbering and initialize solution vector
%
% 3. Compute initial gap function
%
% 4. Compute external force vector
%
% 5. Compute the master stiffness matrix of the structure
%
% 6. Create the expanded system of equations corresponding to the Lagrange Multipliers method
%
% 7. Loop over all contact iterations
% ->
%    7i. Print progress message on the contact iterations
%
%   7ii. Assign inactive DOFs to the ones from previous contact iteration
%
%  7iii. Determine active contact nodes
%
%   7iv. Collect the DOFs of the homogeneous Dirichlet boundary conditions and the active contact nodes of the current contact iteration in one array
%
%    7v. Collect all constrained DOFs into one array
%
%   7vi. Find the free DOFs of the current iteration
%
%  7vii. Solve the equation system of the current contact iteration
%
% 7viii. Evaluate convergence conditions
%
%   7ix. Print progress message on the contact conditions
%
%    7x. Update the iteration counter
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
% 12. Return only the negative Lagrange Multipliers
%
% 13. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________\n');
    fprintf('#############################################################\n');
    fprintf('Compute the displacement field for a plate in membrane action\n');
    fprintf('problem subject to rigid contact boundary has been initiated\n\n');
    if ~isfield(propContact, 'numNodes')
        if ~isinteger(propContact.numNodes)
            error('Variable propContact.numNodes needs to be defined and be an integer');
        end
    end
    fprintf('Number of potential contact nodes : %d\n', propContact.numNodes);
    if ~isfield(segmentsContact, 'numSegments')
        if ~isinteger(segmentsContact.numSegments)
            error('Variable segmentsContact.numSegments needs to be defined and be an integer');
        end
    end
    fprintf('Number of rigid contact segments : %d\n', segmentsContact.numSegments);
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
precompStiffMtx = 'undefined';
precomResVct = 'undefined';

% Steady-state analysis
t = 0;
noTimeStep = 1;

% Connectivity arrays for the DOFs into the resulting vectors
DOF4Output = [1:2:noDOFs-1
              2:2:noDOFs];

% Prescribed DOFs
prescribedDOFs = sort(horzcat(homDOFs, inhomDOFs));

% Tabulation for the output in the command window
tab = '';

% Initialize output array
dHat_stiffMtx = zeros(noDOFs, 1);

% Initialize counter of contact iterations
counterContact = 1;

% Initialize array of active nodes
homDOFsLM = [];

% Inialize convergence conditions
isCnd_DOFs = false;
isCnd_LM = false;

% Title for the output file
title = 'Contact analysis for a plate in membrane action';

% Initialize array containing information about which potentially contact
% nodes are enabled for each rigid segment
propContact.activeNodesToSegment = ...
    false(propContact.numNodes, segmentsContact.numSegments);

% Array of logical strings
LogicalStr = {'false', 'true'};

%% 1. Remove fully constrained nodes
fullyConstrainedNodes = false(propContact.numNodes, 1);
for i = 1:length(propContact.nodeIDs)
    DOFs = 2*propContact.nodeIDs(i, 1) - 1 : 2*propContact.nodeIDs(i, 1);
    if ismember(DOFs(1), prescribedDOFs) && ismember(DOFs(2), prescribedDOFs)
        fullyConstrainedNodes(i,1) = true; 
    end
end
propContact.nodeIDs(fullyConstrainedNodes) = [];
propContact.numNodes = length(propContact.nodeIDs);

%% 2. Get number of Lagrange Multipliers DOFs, total number of DOFs, total DOF numbering and initialize solution vector
numDOFsLM = propContact.numNodes*segmentsContact.numSegments;
numDOFsTotal = noDOFs + numDOFsLM;
freeDOFs = 1:numDOFsTotal;
dHat_stiffMtxLM = zeros(numDOFsTotal, 1);

%% 3. Compute initial gap function
gapFunctionInitial = computeInitialGapFunction...
    (strMsh, propContact, segmentsContact);

%% 4. Compute external force vector
F = computeLoadVctFEMPlateInMembraneAction...
    (strMsh, propAnalysis, propNBC, t, propGaussInt, '');

%% 5. Compute the master stiffness matrix of the structure
if strcmp(outMsg, 'outputEnabled')
    fprintf(strcat(tab, '>> Computing the stiffness matrix of the system\n'));
end
[K, F, minElSize] = computeStiffMtxLoadVct ...
    (propAnalysis, dHat_stiffMtx, uSaved, uDot, uDotSaved, precompStiffMtx, ...
    precomResVct, DOFNumbering, strMsh, F, loadFactor, computeBodyForces, ...
    propStrDynamics, parameters, propGaussInt);

%% 6. Create the expanded system of equations corresponding to the Lagrange Multipliers method
if strcmp(outMsg, 'outputEnabled')
    fprintf(strcat(tab, '>> Creating the expanded system of equations\n'));
end
C = computeLagrangeMultipliersMatrixContact2D ...
    (noDOFs, numDOFsLM, propContact, segmentsContact);
stiffMtxLM = [K  C
              C' zeros(size(C,2))];
clear C;
resVecLM = [F
            - gapFunctionInitial];

%% 7. Loop over all contact iterations
if strcmp(outMsg, 'outputEnabled')
    fprintf(strcat(tab, '>> Loop over all contact iterations\n\n'));
end
while counterContact <= propContact.maxIter && ~(isCnd_DOFs && isCnd_LM)
    %% Debugging
%     graph.index = 1;
%     graph.visualization.geometry = 'current';
%     graph.index = plot_currentConfigurationFEMPlateInMembraneAction...
%         (strMsh,homDOFs,segmentsContact,dHat_stiffMtxLM,graph);
%     close(1);
    %% 7i. Print progress message on the contact iterations
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(strcat(tab, '\t'),'>> Contact iteration %d\n'), ...
            counterContact);
    end
    
    %% 7ii. Assign inactive DOFs to the ones from previous contact iteration
    homDOFsLM_saved = homDOFsLM;
    
    %% 7iii. Determine active contact nodes
    homDOFsLM = findInactiveLagrangeMultipliersContact2D...
        (homDOFsLM, strMsh, noDOFs, dHat_stiffMtxLM, gapFunctionInitial, ...
        segmentsContact, propContact);
    
    %% 7iv. Collect the DOFs of the homogeneous Dirichlet boundary conditions and the active contact nodes of the current contact iteration in one array
    homDOFs_iterate = horzcat(homDOFs, homDOFsLM);
      
    %% 7v. Collect all constrained DOFs into one array
    constrained_DOFs = unique(horzcat(homDOFs_iterate, inhomDOFs));
    
    %% 7vi. Find the free DOFs of the current iteration
    freeDOFs_iterate = freeDOFs;
    freeDOFs_iterate(ismember(freeDOFs_iterate, constrained_DOFs)) = [];
    
    %% 7vii. Solve the equation system of the current contact iteration
    [dHat_stiffMtxLM, FComplete, ~, ~] = solve_FEMLinearSystem ...
        (propAnalysis, uSaved, uDotSaved, uDDotSaved, strMsh, forceVct, ...
        computeBodyForces, parameters, dHat_stiffMtxLM, uDot, uDDot, ...
        massMtx, dampMtx, stiffMtxLM, resVecLM, ...
        computeProblemMatricesSteadyState, DOFNumbering, ...
        freeDOFs_iterate, homDOFs_iterate, inhomDOFs, valuesInhomDOFs, ...
        uMeshALE, solve_LinearSystem, propStrDynamics, propNLinearAnalysis, ...
        propGaussInt, strcat(tab, '\t'), outMsg);
    
    %% 7viii. Evaluate convergence conditions

    % Check whether the vector of active contact nodes has not changed
    isCnd_DOFs = isequal(homDOFsLM_saved, homDOFsLM);
    
    % Check whether all Lagrange Multipliers DOFs are negative
    isCnd_LM = min(dHat_stiffMtxLM(noDOFs + 1:length(resVecLM))) >= 0;
    
    %% 7ix. Print progress message on the contact conditions
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(strcat(tab, '\t'),'>> Convergence criteria isCnd_DOFs(%s) and isCnd_LM(%s)\n\n'), ...
            LogicalStr{isCnd_DOFs + 1}, LogicalStr{isCnd_LM + 1});
    end
    
    %% 7x. Update the iteration counter
    counterContact = counterContact + 1;
end
numInter = counterContact - 1;
if (numInter >= propContact.maxIter)
    warning('\t Max number of contact iterations has been reached!\n');
end

%% 8. Get the displacement solution
dHat = dHat_stiffMtxLM(1:noDOFs);

%% 9. Write out the results into a file if the case is not a unit test
if ~isUnitTest
    if strcmp(outMsg,'outputEnabled')
        fprintf('>> Writting out the results to "%s"\n',strcat(pathToOutput, caseName, '/'));
    end
    writeOutputFEMPlateInMembraneActionToVTK...
        (propAnalysis, propNLinearAnalysis, propStrDynamics, strMsh, parameters, ...
        dHat, uDot, uDDot, DOF4Output, caseName, pathToOutput, title, ...
        noTimeStep);
end

%% 10. Get the Lagrange Multipliers solution
lambdaHat = dHat_stiffMtxLM(noDOFs + 1:numDOFsTotal);

%% 11. Get the index of the active Lagrange Multipliers DOFs
allContactNodes = repmat(propContact.nodeIDs, segmentsContact.numSegments,1);
nodeIDs_active = allContactNodes(lambdaHat > 0);

%% 12. Return only the negative Lagrange Multipliers
lambdaHat = lambdaHat(lambdaHat > 0);

%% 13. Appendix
if strcmp(outMsg, 'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nContact analysis took %.2d seconds \n', computationalTime);
    fprintf('Number of active (contact) nodes : %d\n\n' , length(nodeIDs_active));
    fprintf('____________________Contact Analysis Ended___________________\n');
    fprintf('############################################################\n\n\n');
end

end
