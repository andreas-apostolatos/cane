function [eigenmodeShapes,naturalFrequencies,BSplinePatches,dHat] = ...
    solve_DDMEigenmodeAnalysisIGAMembrane...
    (BSplinePatches,connections,propCoupling,solve_LinearSystem,noEig,...
    propNLinearAnalysis,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the natural frequencies and the coresponding eigenmode shapes for
% the geometrically linear multipatch isogeometric membrane problem.
%
%                Input :
%       BSplinePatches : Structure containing all the information regarding 
%                        the connections between the multipatches
%          connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%                        .lambda : The NURBS discretization of the 
%                                  interface traction force field for each 
%                                  defined interface
%                            .mu : The NURBS discretization of the 
%                                  interface traction moment field for each 
%                                  defined interface
%         propCoupling : Properties of the multipatch coupling
%                           .method : Coupling method
%                           .alphaD : penalty factor for the displacement
%                                     coupling
%                           .alphaR : penalty factor for the rotation
%                                     coupling
%                             .intC : On the integration of the coupling
%                                     interface
%    solve_LinearSystem : Function handle to the solution of linear 
%                         equation systems
%                 noEig : Number of eigenmode shapes to be returned
%   propNLinearAnalysis : Properties of the nonlinear method for
%                            tackling the nonlinear problem :
%                                  .method : 'newtonRapshon'
%                             .noLoadSteps : Number of load steps
%                                     .eps : Residual tolerance
%                                 .maxIter : Maximum number of nonlinear
%                                              iterations
%                outMsg : Enables outputting information onto the command 
%                         window when chosen 'outputEnabled'
%
%                Output :
%       eigenmodeShapes : The eigenmode shapes into a matrix
%    naturalFrequencies : The natural frequencies of the discrete system in 
%                         descending order
%        BSplinePatches : The updated array of the B-Spline patches
%                  dHat : The displacement field when solving the problem 
%                         for bringing it in equilibrium with its internal 
%                         stresses
%
% Function layout :
%
% 0. Read input
%
% 1. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
%
% 2. Compute the total number of DOFs for the multipatch structure
%
% 3. Assign a number of DOFs for each Lagrange Multipliers field corresponding to the multipatch coupling
%
% 4. Compute the total number of DOFs for the coupled multipatch system
%
% 5. Create an element freedom table for each patch and each Lagrange Multipliers fields
%
% 6. Create a DOF numbering for each patch and each Lagrange Multipliers field
%
% 7. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions
%
% 8. Compute the constant problem matrices according to the chosen method
%
% 9. Compute an empty load vector for each patch
%
% 10. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
%
% 11. Initialize the solution field
%
% 12. Perform a static solve to bring the structure in equilibrium with its internal forces
%
% 13. Compute the static external loads
% ->
%     13i. Initialize the load vector and the tangent stiffness matrix
%
%    13ii. Get the Neumann boundary conditions of the current patch
%
%   13iii. Check if there is a non-conservative loading associated with the current patch
%
%    13iv. Loop over the Neumann boundary conditions of the current patch
%    ->
%          13iv.1. Get the function handle for the load vector computation
%
%          13iv.2. Compute the load vector and the tangent matrix resulting from the application of follower loads
%
%          13iv.3. If the loading is not conservative add the contribution to the non-conservative load vector
%
%          13iv.4. Add The compute external load vector into the B-Spline array
%    <-
% <-
%
% 13. Compute the linear stiffness matrix of the membrane problem
%
% 14. Compute the mass matrix of the membrane problem
%
% 15. Solve the generalized eigenvalue problem to get the eigenmode shapes
%
% 16. Apply the boundary conditions onto the system
%
% 17. Appendix
% 
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_______________________________________________________________\n');
    fprintf('###############################################################\n');
    fprintf('Modal analysis for the domain decomposition of the isogeometric\n');
    fprintf('membrane using the ');
    if strcmp(propCoupling.method,'lagrangeMultipliers')
        if isfield(propCoupling,'alphaD')
            if propCoupling.alphaD ~= 0
                fprintf('augmented Lagrange Multipliers\n');
                fprintf('method has been initiated\n\n');
                couplingMethod = 'augmented Lagrange Multipliers';
            else
                fprintf('Lagrange Multipliers method\n');
                fprintf('has been initiated\n\n');
                couplingMethod = 'Lagrange Multipliers';
            end
        end
    elseif strcmp(propCoupling.method,'penalty')
        fprintf('Penalty method has been \n');
        fprintf('initiated\n\n');
        couplingMethod = 'Penalty';
    elseif strcmp(propCoupling.method,'nitsche')
        fprintf('Nitsche method has been \n');
        fprintf('initiated\n\n');
        couplingMethod = 'Nitsche';
    else
        fprintf('\n');
        error('Choose coupling method in propCoupling.method');
    end
    fprintf('Number of eigenfrequencies requested : %d',noEig)
    fprintf('\n\n');
    fprintf('Coupling properties \n');
    fprintf('------------------- \n\n');
    if strcmp(propCoupling.method,'penalty')
        if isfield(propCoupling,'alphaD')
            for iConnections = 1:connections.No
                idI = connections.xiEtaCoup(iConnections,1);
                idJ = connections.xiEtaCoup(iConnections,2);
                fprintf('Penalty factor for the dislacement field across gamma^(%d,%d): %d\n',idI,idJ,propCoupling.alphaD(iConnections,1));
            end
        else
            fprintf('\n');
            error('Choose penalty factor for the displacement field in propCoupling.alphaD');
        end
    elseif strcmp(propCoupling.method,'lagrangeMultipliers')
        if isfield(propCoupling,'alphaD')
            if propCoupling.alphaD ~= 0
                for iConnections = 1:connections.No
                    idI = connections.xiEtaCoup(iConnections,1);
                    idJ = connections.xiEtaCoup(iConnections,2);
                    fprintf('Penalty factor for the dislacement field across gamma^(%d,%d): %d\n',idI,idJ,propCoupling.alphaD(iConnections,1));
                end
            end
        end
        if ~isfield(connections,'lambda')
            fprintf('\n');
            fprintf('Choose a discretization for the Lagrange Multipliers field in connections.lambda');
        end
    elseif strcmp(propCoupling.method,'nitsche')
        if propCoupling.estimationStabilPrm == true;
            fprintf('Automatic estimation of the stabilization is enabled\n');
        else
            if isfield(propCoupling,'alphaD')
                fprintf('Stabilization factor is chosen as %d\n',propCoupling.alphaD);
            else
                 fprintf('\n');
                 error('Choose stabilization parameter in propCoupling.alphaD');
            end
        end
        if isfield(propCoupling,'gammaTilde')
            fprintf('Linear combination factor for the interface tractions is chosen as %d',propCoupling.gammaTilde);
        else
            fprintf('\n');
            error('Choose linear combination factor for the interface tractions in propCoupling.gammaTilde');
        end
    end
    fprintf('\n\n');
    isWeakDBC = false;
    for iPatches = 1:length(BSplinePatches)
        if isfield(BSplinePatches{iPatches},'weakDBC')
            if isfield(BSplinePatches{iPatches}.weakDBC,'noCnd')
                if BSplinePatches{iPatches}.weakDBC.noCnd > 0
                    isWeakDBC = true;
                    break;
                end
            end
        end
    end
    if isWeakDBC
        fprintf('Weak Dirichlet boundary conditions \n');
        fprintf('---------------------------------- \n\n');
        for iPatches = 1:length(BSplinePatches)
            BSplinePatch = BSplinePatches{iPatches};
            if isfield(BSplinePatch,'weakDBC')
                if isfield(BSplinePatch.weakDBC,'noCnd')
                    if BSplinePatch.weakDBC.noCnd > 0
                        fprintf('Weak boundary conditions using the %s method are applied \n',BSplinePatch.weakDBC.method);
                        fprintf('over %d boundaries of the B-Spline patch %d: \n',BSplinePatch.weakDBC.noCnd,iPatches);
                        if strcmp(BSplinePatch.weakDBC.method,'Nitsche')
                            if BSplinePatch.weakDBC.estimationStabilPrm == true
                                fprintf('Automatic estimation of the stabilization parameter is enabled \n');
                            end
                            if isfield(BSplinePatch.weakDBC,'computeConstMtx')
                                if isfield(BSplinePatch.weakDBC.alpha)
                                    fprintf('Manual stabilization parameter chosen as %d\n',BSplinePatch.weakDBC.alpha)
                                else
                                    error('Manual stabilization parameter weakDBC.alpha needs to be assigned\n');
                                end
                            end
                        elseif strcmp(BSplinePatch.weakDBC.method,'Penalty')
                            if isfield(BSplinePatch.weakDBC.alpha)
                                fprintf('The penalty parameter is chosen as %d',BSplinePatch.weakDBC.alpha);
                            else
                                error('The penalty parameter needs to be assigned');
                            end
                        end
                        fprintf('\n');
                    end
                end
            end
        end
    end
    fprintf('\n');
    fprintf('_______________________________________________________________\n');
    fprintf('\n');
    
    % start measuring computational time
    tic;
end

%% 0. Read input

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Tabulation of the output messages
tab = '\t';

% Flag on whether the reference is updated
isReferenceUpdated = false;

% Check if Lagrange Multipliers field for the displacement coupling is
% enforced
isLagrangeMultipliersDisplacementsEnabled = false;
if isfield(connections,'lambda') && strcmp(propCoupling.method,'lagrangeMultipliers') 
    isLagrangeMultipliersDisplacementsEnabled = true;
    if isempty(connections.lambda)
        error('Lagrange multipliers discretization for the displacement coupling exists but its empty');
    end
elseif ~isfield(connections,'lambda') && strcmp(propCoupling.method,'lagrangeMultipliers') 
    error('Lagrange Multipliers method for the coupling is chosen but field connections.lambda is undefined');
end

% Check if Lagrange Multipliers field for the rotational coupling is
% enforced
isLagrangeMultipliersRotationsEnabled = false;
if isfield(connections,'mu') && strcmp(propCoupling.method,'lagrangeMultipliers')
    isLagrangeMultipliersRotationsEnabled = true;
    if isempty(connections.mu)
        error('Lagrange multipliers discretization for the rotational coupling exists but its empty');
    end
end

% Assign the computation of the tangent stiffness matrix according to the selected method
if strcmp(propCoupling.method,'penalty') || strcmp(propCoupling.method,'lagrangeMultipliers')
    computeTanStiffMtxAndResVct = @computeTangentStiffMtxResVctDDMPenaltyIGAMembrane;
elseif strcmp(propCoupling.method,'nitsche')
    computeTanStiffMtxAndResVct = @computeTangentStiffMtxResVctDDMNitscheIGAMembrane;
else
    error('Choose coupling method in propCoupling.method');
end

% Number of patches
noPatches = length(BSplinePatches);

% On the application of weak Dirichlet boundary conditions
isWeakDBC = zeros(noPatches,1);
for iPatches = 1:noPatches
    if ~isempty(BSplinePatches{iPatches}.weakDBC)
        if isfield(BSplinePatches{iPatches}.weakDBC,'noCnd')
            if BSplinePatches{iPatches}.weakDBC.noCnd > 0
                isWeakDBC(iPatches,1) = true;
            end
        end
    end
end

% Check if there exist a nonconservative loading
isConservative = true(noPatches,1);
for iPatches = 1:noPatches
    for iNBC = 1:BSplinePatches{iPatches}.NBC.noCnd
        if BSplinePatches{iPatches}.NBC.isFollower(iNBC,1)
            isConservative(iPatches,1) = false;
            break;
        end
    end
end

% Counter of the nonlinear iterations
counterNonlinearIterations = 1;

% Initialize the dummy arrays
dHatDot = 'undefined';
dHatDDot = 'undefined';
dHatDDotSaved = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
noTimeStep = 'undefined';
noPatch = 'undefined';
plot_IGANLinear = 'undefined';
graph = 'undefined';

% Set co-simulation with Empire to false
isCosimulationWithEmpire = false;

% Initialize tangent matrix contribution from follower load
tanMtxLoad = struct([]);

% Initialize the time
t = 0;

% Define function handle for the computation of the updated geometry
computeUpdatedGeometry = @computeUpdatedGeometryIGAThinStructureMultipatches;

% Get the number of weak Dirichlet boundary conditions
noWeakDBCCnd = 0;
for iPatches = 1:noPatches
    if ~isempty(BSplinePatches{iPatches}.weakDBC)
        if isfield(BSplinePatches{iPatches}.weakDBC,'noCnd')
            noWeakDBCCnd = noWeakDBCCnd + ...
                BSplinePatches{iPatches}.weakDBC.noCnd;
        end
    end
end

% Static linear analysis
propTransientAnalysis.timeDependence = 'steadyState';
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;

% Load factor definition
loadFactor = 1;

% Initialize auxiliary variables
sizePrevious = 0;

%% 1. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.noDOFs = 3*BSplinePatches{iPatches}.noCPs;
    noDOFsPatchLM = 0;
    if isWeakDBC(iPatches,1)
        if strcmp(BSplinePatches{iPatches}.weakDBC.method,'lagrangeMultipliers')
            for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
                noDOFsPatchLMCnd = 3*length(BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.CP(:,1));
                noDOFsPatchLM = noDOFsPatchLM + noDOFsPatchLMCnd;
                if iCnd == 1
                    index = 3*BSplinePatches{iPatches}.noCPs;
                else
                    index = BSplinePatches{iPatches}.weakDBC.lambda{iCnd-1}.EFT(length(BSplinePatches{iPatches}.weakDBC.lambda{iCnd-1}.EFT));
                end
                BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.EFT = index + 1:index + noDOFsPatchLMCnd;
            end
        end
    end
    BSplinePatches{iPatches}.noDOFs = BSplinePatches{iPatches}.noDOFs + noDOFsPatchLM;
end

%% 2. Compute the total number of DOFs for the multipatch structure
noDOFsPatches = 0;
for iPatches = 1:noPatches
    noDOFsPatches = noDOFsPatches + BSplinePatches{iPatches}.noDOFs;
end

%% 3. Assign a number of DOFs for each Lagrange Multipliers field corresponding to the multipatch coupling
if strcmp(couplingMethod,'Lagrange Multipliers') || strcmp(couplingMethod,'augmented Lagrange Multipliers')
    noDOFsLagrangeMultipliers = 0;
    for iConnections = 1:connections.No
        noDOFsLagrangeMultipliers = noDOFsLagrangeMultipliers + ...
            3*connections.lambda{iConnections}.noCPs;
        if isLagrangeMultipliersRotationsEnabled
            noDOFsLagrangeMultipliers = noDOFsLagrangeMultipliers + ...
            2*connections.mu{iConnections}.noCPs;
        end
    end
end

%% 4. Compute the total number of DOFs for the coupled multipatch system
if strcmp(couplingMethod,'Lagrange Multipliers') || strcmp(couplingMethod,'augmented Lagrange Multipliers')
    noDOFs = noDOFsPatches + noDOFsLagrangeMultipliers;
else
    noDOFs = noDOFsPatches ;
end
if noEig > noDOFs
    error('The number of  requested eigenfrequencies can be up to %',noDOFs);
end

%% 5. Create an element freedom table for each patch and each Lagrange Multipliers fields

% For the patches :
% _________________

for iPatches = 1:noPatches
    noDOFsPatch = BSplinePatches{iPatches}.noDOFs;
    if iPatches == 1
        BSplinePatches{iPatches}.EFTPatches = 1:noDOFsPatch;
    else
        BSplinePatches{iPatches}.EFTPatches = ...
            BSplinePatches{iPatches-1}.EFTPatches(length(BSplinePatches{iPatches-1}.EFTPatches)) + ...
            1:BSplinePatches{iPatches-1}.EFTPatches(length(BSplinePatches{iPatches-1}.EFTPatches)) + ...
            noDOFsPatch;
    end
end

% For the Lagrange Multipliers fields :
% _____________________________________

if strcmp(couplingMethod,'Lagrange Multipliers') || strcmp(couplingMethod,'augmented Lagrange Multipliers')
    for iConnections = 1:connections.No
        if iConnections == 1
            index = noDOFsPatches;
            connections.lambda{iConnections}.EFTLagrangeMultipliers = ...
                index+1:index+3*connections.lambda{iConnections}.noCPs;
        else
            if isLagrangeMultipliersRotationsEnabled
                index = connections.mu{iConnections-1}.EFTLagrangeMultipliers(length(connections.mu{iConnections-1}.EFTLagrangeMultipliers));
            else
                index = connections.lambda{iConnections-1}.EFTLagrangeMultipliers(length(connections.lambda{iConnections-1}.EFTLagrangeMultipliers));
            end
            connections.lambda{iConnections}.EFTLagrangeMultipliers = index + 1:index+3*connections.lambda{iConnections}.noCPs;
        end
        if isLagrangeMultipliersRotationsEnabled
            index = connections.lambda{iConnections}.EFTLagrangeMultipliers(length(connections.lambda{iConnections}.EFTLagrangeMultipliers));
            connections.mu{iConnections}.EFTLagrangeMultipliers = index + 1:index + 2*connections.mu{iConnections}.noCPs;
        end
    end
end

%% 6. Create a DOF numbering for each patch and each Lagrange Multipliers field

% For the patches :
% _________________

for iPatches = 1:noPatches
    % Get the number of Control Points in xi-direction
    nxi = length(BSplinePatches{iPatches}.CP(:,1,1));
    
    % Get the number of Control Points in eta-direction
    neta = length(BSplinePatches{iPatches}.CP(1,:,1));
    
    % Initialize the DOF numbering array
    BSplinePatches{iPatches}.DOFNumbering = zeros(nxi,neta,3);
    
    % Compute the entries of the DOF numbering array
    k = 1;
    for cpj = 1:neta
        for cpi = 1:nxi
            BSplinePatches{iPatches}.DOFNumbering(cpi,cpj,1) = k;
            BSplinePatches{iPatches}.DOFNumbering(cpi,cpj,2) = k + 1;
            BSplinePatches{iPatches}.DOFNumbering(cpi,cpj,3) = k + 2;

            % Update counter
            k = k + 3;
        end
    end
    
    % Create a DOF numbering for each Lagrange Multipliers field within the
    % patch
    if isWeakDBC(iPatches,1)
        if strcmp(BSplinePatches{iPatches}.weakDBC.method,'lagrangeMultipliers')
            for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
                % Get the number of Control Points
                noXiLambda = length(BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.Xi);

                % Initialize the field of the DOF numbering
                BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.DOFNumbering = zeros(noXiLambda,3);

                 % Compute the entries of the DOF numbering array
                k = 1;
                for cpi = 1:noXiLambda
                    BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.DOFNumbering(cpi,1) = k;
                    BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.DOFNumbering(cpi,2) = k + 1;
                    BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.DOFNumbering(cpi,3) = k + 2;

                    % Update counter
                    k = k + 3;
                end
            end
        end
    end
end

% For the Lagrange Multipliers field :
% ____________________________________

if strcmp(couplingMethod,'Lagrange Multipliers') || strcmp(couplingMethod,'augmented Lagrange Multipliers')
    for iConnections = 1:connections.No
        % Get the number of Control Points in xi-direction
        if isLagrangeMultipliersDisplacementsEnabled
            nxiLambda = length(connections.lambda{iConnections}.CP(:,1));
        end
        if isLagrangeMultipliersRotationsEnabled
            nxiMu = length(connections.mu{iConnections}.CP(:,1));
        end

        % Lagrange multipliers for the traction forces :
        % ______________________________________________

        if isLagrangeMultipliersDisplacementsEnabled
            % Initialize the field
            connections.lambda{iConnections}.DOFNumbering = zeros(nxiLambda,3);

            % Compute the entries of the DOF numbering array
            k = 1;
            for cpi = 1:nxiLambda
                connections.lambda{iConnections}.DOFNumbering(cpi,1) = k;
                connections.lambda{iConnections}.DOFNumbering(cpi,2) = k + 1;
                connections.lambda{iConnections}.DOFNumbering(cpi,3) = k + 2;

                % Update counter
                k = k + 3;
            end
        end
        % Lagrange multipliers for the traction moments :
        % _______________________________________________

        if isLagrangeMultipliersRotationsEnabled
            % Initialize the field
            connections.mu{iConnections}.DOFNumbering = zeros(nxiMu,2);

            % Compute the entries of the DOF numbering array
            k = 1;
            for cpi = 1:nxiMu
                connections.mu{iConnections}.DOFNumbering(cpi,1) = k;
                connections.mu{iConnections}.DOFNumbering(cpi,2) = k + 1;

                % Update counter
                k = k + 2;
            end
        end
    end
end

%% 7. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions
if strcmp(outMsg,'outputEnabled')
    message = 'Computing the constant problem matrices for each patch for the application of weak Dirichlet boundary conditions\n';
    fprintf([tab,'>>',' ',message]);
end
for iPatches = 1:noPatches
    if isWeakDBC(iPatches,1)
        isConstMtxWeakDBC = true;
        if isfield(BSplinePatches{iPatches}.weakDBC,'method')
            if strcmp(BSplinePatches{iPatches}.weakDBC.method,'penalty') || ...
                    (strcmp(BSplinePatches{iPatches}.weakDBC.method,'nitsche') && ...
                    ~BSplinePatches{iPatches}.weakDBC.estimationStabilPrm)
                computeWeakDBCConstantProblemMatrices = @computeWeakDBCMtxPenaltyIGAMembrane;
            elseif strcmp(BSplinePatches{iPatches}.weakDBC.method,'lagrangeMultipliers')
                computeWeakDBCConstantProblemMatrices = @computeWeakDBCMtxLagrangeMultipliersIGAMembrane;
            elseif strcmp(BSplinePatches{iPatches}.weakDBC.method,'nitsche') && ...
                    BSplinePatches{iPatches}.weakDBC.estimationStabilPrm
                isConstMtxWeakDBC = false;
            else
                error('Define a valid method for the application of weak Dirichlet boundary conditions in BSplinePatches{iPatches}.weakDBC.method');
            end
            if isConstMtxWeakDBC
                BSplinePatches{iPatches}.KConstant = ...
                    computeWeakDBCConstantProblemMatrices ...
                    (BSplinePatches{iPatches},connections,...
                    BSplinePatches{iPatches}.noDOFs,propCoupling);
            else
                BSplinePatches{iPatches}.KConstant = 'undefined';
            end
        else
            error('Define variable method for the application of weak Dirichlet boundary conditions in BSplinePatches{iPatches}.weakDBC');
        end
    else
        BSplinePatches{iPatches}.KConstant = 'undefined';
    end
end

%% 8. Compute the constant problem matrices according to the chosen method
if strcmp(outMsg,'outputEnabled')
    message = 'Computing the constant problem matrices for the multipatch coupling\n';
    fprintf([tab,'>>',' ',message]);
end
if strcmp(couplingMethod,'Penalty') || ...
        (strcmp(couplingMethod,'Nitsche') && ~propCoupling.estimationStabilPrm) || ...
        (~(strcmp(couplingMethod,'Lagrange Multipliers') || strcmp(couplingMethod,'augmented Lagrange Multipliers')) && isfield(propCoupling,'alphaR'))
    if (strcmp(couplingMethod,'Nitsche') && ~propCoupling.estimationStabilPrm)
        if isfield(propCoupling,'alphaD')
            error('For automatic estimation of the stabilization parameter corresponding to the multipatch coupling no displacement penalization needs to be specified in couplingMethod.alphaD');
        end
    end
    constStiff = computeConstantMtxForDDMPenaltyIGAThinStructure...
        (BSplinePatches,connections,noDOFs,propCoupling);
elseif strcmp(couplingMethod,'Lagrange Multipliers') || strcmp(couplingMethod,'augmented Lagrange Multipliers')
    constStiff = computeConstantMtxForDDMLagrangeMultipliersIGAThinStructure...
        (BSplinePatches,connections,noDOFs,propCoupling);
else
    constStiff = 'undefined';
end

%% 9. Compute an empty load vector for each patch
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.FGamma = zeros(3*BSplinePatches{iPatches}.noCPs,1);
end

%% 10. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
homDOFs = [];
inhomDOFs = [];
valuesInhomDOFs = [];
masterDOFs = [];
slaveDOFs = [];
for iPatches = 1:noPatches
    if iPatches ~= 1
        % Determine the number of the DOFs of the previous patches
        sizePrevious = sizePrevious + 3*BSplinePatches{iPatches-1}.noCPs;
        
        % Add the numbering DOFs where homogeneous Dirichlet boundary
        % conditions are applied from the patch level
        homDOFsPatch = sizePrevious + BSplinePatches{iPatches}.homDOFs;
        homDOFs = mergesorted(homDOFs,homDOFsPatch);
        
        % Add the numbering DOFs where inhomogeneous Dirichlet boundary
        % conditions are applied from the patch level
        inhomDOFsPatch = sizePrevious + BSplinePatches{iPatches}.inhomDOFs;
        inhomDOFs = mergesorted(inhomDOFs,inhomDOFsPatch);
        
        % Add the prescribed values to the DOFs where inhomogeneous 
        % Dirichlet boundary conditions are applied from the patch level
%         valuesInhomDOFsPatch = sizePrevious + BSplinePatches{counterPatches}.valuesInhomDOFs;
%         valuesInhomDOFs = mergesorted(valuesInhomDOFs,valuesInhomDOFsPatch);
        valuesInhomDOFsPatch = BSplinePatches{iPatches}.valuesInhomDOFs;
        if ~isempty(valuesInhomDOFsPatch)
            valuesInhomDOFs = horzcat(valuesInhomDOFs,valuesInhomDOFsPatch);
        end
        
        % Add the numbering of the master DOFs in the patch level
        masterDOFsPatch = sizePrevious + BSplinePatches{iPatches}.masterDOFs;
        masterDOFs = mergesorted(masterDOFs,masterDOFsPatch);
        
        % Add the numbering of the master DOFs in the patch level
        slaveDOFsPatch = sizePrevious + BSplinePatches{iPatches}.slaveDOFs;
        slaveDOFs = mergesorted(slaveDOFs,slaveDOFsPatch);
    else
        homDOFs = BSplinePatches{iPatches}.homDOFs;
        inhomDOFs = BSplinePatches{iPatches}.inhomDOFs;
        valuesInhomDOFs = BSplinePatches{iPatches}.valuesInhomDOFs;
        masterDOFs = BSplinePatches{iPatches}.masterDOFs;
        slaveDOFs = BSplinePatches{iPatches}.slaveDOFs;
    end
end

% Vector of the global numbering of the unconstrained DOFs
freeDOFs = 1:noDOFs;
freeDOFs(ismember(freeDOFs,homDOFs)) = [];
freeDOFs(ismember(freeDOFs,inhomDOFs)) = [];

%% 11. Initialize the solution field
dHat = zeros(noDOFs,1);

%% 12. Perform a static solve to bring the structure in equilibrium with its internal forces
if strcmp(outMsg,'outputEnabled')
    message = 'Performing steady-state analysis to bring the structure in equilibrium\n';
    fprintf([tab,'>>',' ',message]);
end
[dHat,~,~,~,~,~,~,~,BSplinePatches,propCoupling,~] = ...
    solve_IGANLinearSystem...
    (analysis,dHatSaved,dHatDotSaved,dHatDDotSaved,BSplinePatches,connections,...
    dHat,dHatDot,dHatDDot,constStiff,massMtx,dampMtx,computeTanStiffMtxAndResVct,...
    computeUpdatedGeometry,freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,...
    masterDOFs,slaveDOFs,solve_LinearSystem,t,propCoupling,propTransientAnalysis,...
    propNLinearAnalysis,plot_IGANLinear,isReferenceUpdated,isCosimulationWithEmpire,...
    strcat(tab,'\t'),graph,outMsg);

%% 13. Compute the static external loads
for iPatches = 1:noPatches
    %% 13i. Initialize the load vector and the tangent stiffness matrix
    FGamma = zeros(BSplinePatches{iPatches}.noDOFs,1);
    tanMtxLoad{iPatches} = zeros(BSplinePatches{iPatches}.noDOFs);
    BSplinePatches{iPatches}.FGamma = FGamma;
    
    %% 13ii. Get the Neumann boundary conditions of the current patch
    NBC = BSplinePatches{iPatches}.NBC;
    
    %% 13iii. Check if there is a non-conservative loading associated with the current patch
    if ~isConservative(iPatches,1)
        BSplinePatches{iPatches}.FNonConservative = ...
            zeros(BSplinePatches{iPatches}.noDOFs,1);
    end

    %% 13iv. Loop over the Neumann boundary conditions of the current patch
    for iNBC = 1:NBC.noCnd
        %% 13iv.1. Get the function handle for the load vector computation
        funcHandle = str2func(NBC.computeLoadVct{iNBC});
        
        %% 13iv.2. Compute the load vector and the tangent matrix resulting from the application of follower loads
        if ~(propTransientAnalysis.isStaticStep && NBC.isTimeDependent(iNBC,1))
            [FGamma,tanMtxLoadPatch] = funcHandle(FGamma,BSplinePatches{iPatches},...
                NBC.xiLoadExtension{iNBC},NBC.etaLoadExtension{iNBC},...
                NBC.loadAmplitude{iNBC},NBC.loadDirection{iNBC},...
                NBC.isFollower(iNBC,1),t,BSplinePatches{iPatches}.int,'');
            if NBC.isFollower(iNBC,1)
                tanMtxLoad{iPatches} = tanMtxLoad{iPatches} + tanMtxLoadPatch;
            end
        end
        
        %% 13iv.3. If the loading is not conservative add the contribution to the non-conservative load vector
        if NBC.isFollower(iNBC,1)
            BSplinePatches{iPatches}.FNonConservative = BSplinePatches{iPatches}.FNonConservative + ...
                FGamma;
        end
        
        %% 13iv.4. Add The compute external load vector into the B-Spline array
        BSplinePatches{iPatches}.FGamma = BSplinePatches{iPatches}.FGamma +...
            FGamma;
    end
end

%% 14. Compute the linear stiffness matrix of the membrane problem
if strcmp(outMsg,'outputEnabled')
    message = 'Computing the linear stiffness matrix of the structure\n';
    fprintf([tab,'>>',' ',message]);
end
stiffMtxLinear = computeTanStiffMtxAndResVct(constStiff,tanMtxLoad,dHat,dHatSaved,...
    dHatDot,dHatDotSaved,BSplinePatches,connections,propCoupling,loadFactor,...
    noPatch,noTimeStep,counterNonlinearIterations,noWeakDBCCnd,...
    t,propTransientAnalysis,isReferenceUpdated,strcat(tab,'\t'),outMsg);

%% 15. Compute the mass matrix of the membrane problem
if strcmp(outMsg,'outputEnabled')
    message = 'Computing the mass matrix of the structure\n';
    fprintf([tab,'>>',' ',message]);
end
massMtx = computeIGAMassMtxThinStructure(BSplinePatches,noDOFs);

%% 16. Solve the generalized eigenvalue problem to get the eigenmode shapes
if strcmp(outMsg,'outputEnabled')
    message = 'Solve the generalized eigenvalue problem to get the eigenmode shapes\n';
    fprintf([tab,'>>',' ',message]);
end
[eigenmodeShapesDBC,naturalFrequencies] = ...
    eigs(stiffMtxLinear(freeDOFs,freeDOFs),massMtx(freeDOFs,freeDOFs),noEig,'sm');
naturalFrequencies = naturalFrequencies*ones(length(naturalFrequencies),1);
naturalFrequencies = sqrt(naturalFrequencies);
naturalFrequencies = naturalFrequencies./(2*pi);
[noModeShapes,m] = size(eigenmodeShapesDBC);

%% 17. Apply the boundary conditions onto the system
if strcmp(outMsg,'outputEnabled')
    message = 'Applying the Dirichlet boundary conditions at each mode shape\n\n';
    fprintf([tab,'>>',' ',message]);
end
eigenmodeShapes = zeros(noModeShapes,m);
for iModes = 1:noEig
    eigenmodeShapes(freeDOFs,iModes) = eigenmodeShapesDBC(:,iModes);
    eigenmodeShapes(homDOFs,iModes) = 0;
    eigenmodeShapes(inhomDOFs,iModes) = 0;
end

%% 18. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Modal analysis took %.2d seconds \n\n',computationalTime);
    fprintf('______________________Modal Analysis Ended_____________________\n');
    fprintf('##############################################################\n\n\n');
end

end
