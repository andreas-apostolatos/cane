function [BSplinePatches, CPHistoryMultipatch, propCoupling, resHistory, ...
    isConverged, noIter] = ...
    solve_DDMFormFindingIGAMembrane ...
    (BSplinePatches, connections, propCoupling, propFormFinding, ...
    solve_LinearSystem, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-Spline multipatch array with the updated Control Point
% coordinates performing a form-finding analysis.
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
%             propFormFinding : Properties on the form finding analyse :
%                                .tolerance : Tolerance on the displacement 
%                                             difference between the 
%                                             different form finding
%                                             iterations
%                                .maxNoIter : Maximum number of iterations 
%                                             for convergence of the form 
%                                             finding iterations
%                                .minNoIter : Minimum number of
%                                             form-finding iterations
%         solve_LinearSystem : Function handle to the equation system
%                              solver
%                     outMsg : Enables outputting information onto the
%                              command window when chosen 'outputEnabled'
%
%                     Output :
%             BSplinePatches : The updated B-Spline patch array with the
%                              form-finding results
%        CPHistoryMultipatch : Array containing the control point history
%                              of each patch within the form-finding
%                              iterations
%               propCoupling : Updated array of the coupling properties
%                              with the values of the stabilization factor
%                              in case the Nitsche method for the
%                              multipatch coupling or application of weak
%                              Dirichlet boundary conditions is enabled
%                 resHistory : The history of the error for each form
%                              finding iteration
%                isConverged : Flag on whether the form-finding iterations
%                              have converged
%                     noIter : Number of iterations
%
% Function layout :
%
% 0. Read input
%
% 1. Assign the computation of the tangent stiffness matrix according to the selected method
%
% 2. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
%
% 3. Compute the total number of DOFs for the multipatch structure
%
% 4. Create a freedom table for each patch in the multipatch geometry
%
% 5. Assign a number of DOFs for each Lagrange Multipliers field corresponding to the multipatch coupling
%
% 6. Find the global numbering of the DOFs which are associated to displacements and not to the Lagrange Multipliers
%
% 7. Create a freedom table for each Lagrange Multipliers field corresponding to each patch connection
%
% 8. Compute the total number of DOFs for the coupled multipatch system
%
% 9. Loop over all patches and change parameters corresponding to the form-finding analysis
% ->
%    9i. Change the patch parameters and save the original values
%
%   9ii. Change the load parameters for each patch
%
%  9iii. Change the parameters of the cable elements
% <-
%
% 10. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
%
% 11. Find the master/slave and the domain DOFs
%
% 12. Create a DOF numbering for each patch and each Lagrange Multipliers field
%
% 13. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions
%
% 14. Compute the constant problem matrices according to the chosen method
%
% 15. Initialize the displacement field
%
% 16. Loop over all the form-finding iterations
% ->
%     16i. Solve the static linear problem
%
%    16ii. Check if the solution is not NaN
%
%   16iii. Update the Control Point coordinates of the B-Spline patch
%
%    16iv. Compute the difference between the displacements from both iteration steps
%
%     16v. Compute the residual and check convergence
%
%    16vi. Save the displacement field of the previous form-finding step
%
%   16vii. Update form-finding iteration counter and the pseudotime
% <-
%
% 17. Re-assign arrays
%
% 18. Appendix
%
%% Function main body
if ~strcmp(propCoupling.method, 'penalty') && ~strcmp(propCoupling.method, 'lagrangeMultipliers') && ...
        ~strcmp(propCoupling.method, 'mortar') && ~strcmp(propCoupling.method, 'nitsche')
    error('Wrong coupling method for the multipatch structure in propCoupling.method was selected');
end
if strcmp(outMsg, 'outputEnabled')
    fprintf('______________________________________________________________\n');
    fprintf('##############################################################\n');
    fprintf('Form-finding analysis for the domain decomposition of the\n');
    fprintf('isogeometric membrane using the ');
    if strcmp(propCoupling.method, 'lagrangeMultipliers')
        if isfield(propCoupling, 'alphaD')
            if propCoupling.alphaD ~= 0
                fprintf('augmented Lagrange Multipliers\n');
                fprintf('method has been initiated\n\n');
            else
                fprintf('Lagrange Multipliers method\n');
                fprintf('has been initiated\n\n');
            end
        end
    elseif strcmp(propCoupling.method, 'mortar')
        fprintf('mortar method has been \n');
        fprintf('initiated\n\n');
    elseif strcmp(propCoupling.method, 'penalty')
        fprintf('Penalty method has been \n');
        fprintf('initiated\n\n');
    elseif strcmp(propCoupling.method, 'nitsche')
        fprintf('Nitsche method has been \n');
        fprintf('initiated\n\n');
    else
        fprintf('\n');
        error('Choose coupling method in propCoupling.method');
    end
    fprintf('Form-finding properties \n');
    fprintf('----------------------- \n\n');
    if isfield(propFormFinding, 'tolerance')
        fprintf('Tolerance = %d\n', propFormFinding.tolerance);
    else
        error('Define form-finding tolerance in propFormFinding.tolerance');
    end
    if isfield(propFormFinding, 'minNoIter')
        fprintf('Minimum number of iterations = %d\n', propFormFinding.minNoIter);
    else
        error('Define the minimum number of form-finding iterations in propFormFinding.minNoIter');
    end
    if isfield(propFormFinding, 'maxNoIter')
        fprintf('Maximum number of iterations = %d\n', propFormFinding.maxNoIter);
    else
        error('Define the maximum number of form-finding iterations in propFormFinding.maxNoIter');
    end
    fprintf('\n\n');
    fprintf('Coupling properties \n');
    fprintf('------------------- \n\n');
    if strcmp(propCoupling.method, 'penalty')
        if ~isfield(propCoupling, 'alphaD')
            error('Choose penalty factors for the displacement field in propCoupling.alphaD');
        end
        if size(propCoupling.alphaD) ~= connections.No
            error('The number of the penalty factors in propCoupling.alphaD must be the same as for the number of connections connections.No');
        end
        for iPatches = 1:connections.No
            fprintf('Displacement penalty factor for patch %d equal to alpha = %d \n', ...
                iPatches, propCoupling.alphaD(iPatches, 1));
        end
    elseif strcmp(propCoupling.method,'lagrangeMultipliers')
        if ~isfield(propCoupling, 'alphaD')
            error('Choose penalty factors for the displacement field in propCoupling.alphaD');
        end
        if size(propCoupling.alphaD) ~= connections.No
            error('The number of the penalty factors in propCoupling.alphaD must be the same as for the number of connections connections.No');
        end
        for iPatches = 1:connections.No
            fprintf('Displacement penalty factor for patch %d equal to alpha = %d \n', ...
                iPatches, propCoupling.alphaD(iPatches, 1));
        end
        if ~isfield(connections, 'lambda')
            fprintf('\n');
            fprintf('Choose a discretization for the Lagrange Multipliers field in connections.lambda');
        end
    elseif strcmp(propCoupling.method, 'mortar')
        if ~isfield(propCoupling, 'isSlaveSideCoarser')
            error('Define if the coarse side is the slave in propCoupling.isSlaveSideCoarser by true or false');
        else
            if propCoupling.isSlaveSideCoarser
                fprintf('The coarser side is chosen as slave and the finer as master \n');
            else
                fprintf('The finer side is chosen as slave and the coarser as master \n');
            end
        end
        if ~isfield(propCoupling,'computeRearrangedProblemMtrcs')
            error('Define the function handle to the computation of the re-arranged system corresponding to the mortar method in propCoupling.computeRearrangedProblemMtrcs');
        else
            fprintf('The function handle to the computation of the re-arranged system is \n%s\n', ...
                func2str(propCoupling.computeRearrangedProblemMtrcs));
        end
    elseif strcmp(propCoupling.method, 'nitsche')
        if propCoupling.estimationStabilPrm == true
            fprintf('Automatic estimation of the stabilization is enabled\n');
        else
            if isfield(propCoupling, 'alphaD')
                fprintf('Stabilization factor is chosen as %d\n', propCoupling.alphaD);
            else
                 fprintf('\n');
                 error('Choose stabilization parameter in propCoupling.alphaD');
            end
        end
        if isfield(propCoupling, 'gammaTilde')
            fprintf('Linear combination factor for the interface tractions is chosen as %d', ...
                propCoupling.gammaTilde);
        else
            fprintf('\n');
            error('Choose linear combination factor for the interface tractions in propCoupling.gammaTilde');
        end
    end
    fprintf('\n\n');
    isWeakDBC = false;
    for iPatches = 1:length(BSplinePatches)
        if isfield(BSplinePatches{iPatches}, 'weakDBC')
            if isfield(BSplinePatches{iPatches}.weakDBC, 'noCnd')
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
            if isfield(BSplinePatch, 'weakDBC')
                if isfield(BSplinePatch.weakDBC, 'noCnd')
                    if BSplinePatch.weakDBC.noCnd > 0
                        fprintf('Weak boundary conditions using the %s method are applied \n', BSplinePatch.weakDBC.method);
                        fprintf('over %d boundaries of the B-Spline patch %d: \n', BSplinePatch.weakDBC.noCnd, iPatches);
                        if strcmp(BSplinePatch.weakDBC.method, 'Nitsche')
                            if BSplinePatch.weakDBC.estimationStabilPrm == true
                                fprintf('Automatic estimation of the stabilization parameter is enabled \n');
                            end
                            if isfield(BSplinePatch.weakDBC, 'computeConstMtx')
                                if isfield(BSplinePatch.weakDBC.alpha)
                                    fprintf('Manual stabilization parameter chosen as %d\n', BSplinePatch.weakDBC.alpha)
                                else
                                    error('Manual stabilization parameter weakDBC.alpha needs to be assigned\n');
                                end
                            end
                        elseif strcmp(BSplinePatch.weakDBC.method, 'Penalty')
                            if isfield(BSplinePatch.weakDBC.alpha)
                                fprintf('The penalty parameter is chosen as %d', BSplinePatch.weakDBC.alpha);
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
    fprintf('______________________________________________________________\n');
    fprintf('\n');
    tic;
end

%% 0. Read input

% Analysis type
analysis.type = 'isogeometricMembraneNonlinear';

% Tabulation of the output messages
tab = '\t';

% Number of patches
noPatches = length(BSplinePatches);

% Check input
for iPatches = 1:noPatches
    if isfield(BSplinePatches{iPatches}.weakDBC, 'computeTangMtxResVct') && ...
            (strcmp(BSplinePatches{iPatches}.weakDBC.method, 'penalty') || ...
            strcmp(BSplinePatches{iPatches}.weakDBC.method, 'lagrangeMultipliers'))
        error('No field weakDBC.%s needs to be defined for the method %s corresponding to the application of weak Dirichlet boundary conditions', ...
            BSplinePatches{iPatches}.weakDBC.computeTangMtxResVct, BSplinePatches{iPatches}.weakDBC.method);
    end
end

% Form-find material and analysis parameter
EYoungFF = 0.0;
poissonRatioFF = 0.0;
densityFF = 0.0;
isFollowerFF = true;

% Check if Lagrange Multipliers field for the displacement coupling is
% enforced
isLagrangeMultipliersDisplacementsEnabled = false;
if (isfield(connections, 'lambda') && strcmp(propCoupling.method, 'lagrangeMultipliers') ) || strcmp(propCoupling.method, 'mortar')
    isLagrangeMultipliersDisplacementsEnabled = true;
    if isfield(connections, 'lambda')
        if isempty(connections.lambda) && strcmp(propCoupling.method, 'lagrangeMultipliers')
            error('Lagrange multipliers discretization for the displacement coupling exists but its empty');
        end
    end
elseif ~isfield(connections, 'lambda') && strcmp(propCoupling.method, 'lagrangeMultipliers')
    error('Lagrange Multipliers method for the coupling is chosen but field connections.lambda is undefined');
end

% Check if Lagrange Multipliers field for the rotational coupling is
% enforced
isLagrangeMultipliersRotationsEnabled = false;
if isfield(connections, 'mu') && strcmp(propCoupling.method, 'lagrangeMultipliers')
    isLagrangeMultipliersRotationsEnabled = true;
    if isempty(connections.mu)
        error('Lagrange multipliers discretization for the rotational coupling exists but its empty');
    end
end

% On the application of weak Dirichlet boundary conditions
isWeakDBC = zeros(noPatches, 1);
for iPatches = 1:noPatches
    if ~isempty(BSplinePatches{iPatches}.weakDBC)
        if isfield(BSplinePatches{iPatches}.weakDBC, 'noCnd')
            if BSplinePatches{iPatches}.weakDBC.noCnd > 0
                isWeakDBC(iPatches, 1) = true;
            end
        end
    end
end

% Assign the the flag corresponding to whether the reference geometry is
% updated
isReferenceUpdated = true;

% Initialize convergence flag to false
isConverged = false;

% Initialize the history of the Control Points through the form finding
% iterations and the residual history
CPHistoryMultipatch = struct([]);
for iPatches = 1:noPatches
    CPHistoryMultipatch{iPatches}.CPHistory{1} = BSplinePatches{iPatches}.CP;
end
resHistory = zeros(propFormFinding.maxNoIter, 1);

% Initialize form finding iteration counter
counterFoFiIter = 1;

% Initialize the dummy arrays
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
dHatDot = 'undefined';
dHatDDot = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
updateDirichletBCs = 'undefined';
propIDBC = 'undefined';

% Set flag for co-simulation with EMPIRE to false
isCosimulationWithEmpire = false;

% Set graph properties
propGraph.index = 1;
propGraph.postprocConfig = 'current';
propGraph.resultant = 'displacement';
propGraph.component = '2norm';

% On plotting the deformed configuration through the nonlinear iterations
plot_IGANLinear = '';

% Pseudo-transient analysis only for storing the coupling data
propTransientAnalysis.timeDependence = 'pseudotransient';
propTransientAnalysis.TStart = 0;
propTransientAnalysis.TEnd = propFormFinding.maxNoIter;
propTransientAnalysis.dt = (propTransientAnalysis.TEnd - propTransientAnalysis.TStart)/...
    propFormFinding.maxNoIter;
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;
t = propTransientAnalysis.TStart;

% Nonlinear analysis
propNLinearAnalysis.method = 'newtonRapshon';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 1e-3;
propNLinearAnalysis.maxIter = 2;

% Initialize auxiliary arrays
parameters = struct([]);
NBC = struct([]);
parametersCablesSaved = struct([]);

%% 1. Assign the computation of the tangent stiffness matrix according to the selected method
if strcmp(propCoupling.method, 'penalty') || strcmp(propCoupling.method, 'lagrangeMultipliers') || ...
        strcmp(propCoupling.method, 'mortar')
    computeTanStiffMtxAndResVct = @computeTangentStiffMtxResVctDDMPenaltyIGAMembrane;
elseif strcmp(propCoupling.method, 'nitsche')
    computeTanStiffMtxAndResVct = @computeTangentStiffMtxResVctDDMNitscheIGAMembrane;
else
    error('Choose coupling method in propCoupling.method');
end

%% 2. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.noDOFs = 3*BSplinePatches{iPatches}.noCPs;
    noDOFsPatchLM = 0;
    if isWeakDBC(iPatches,1)
        if strcmp(BSplinePatches{iPatches}.weakDBC.method, 'lagrangeMultipliers')
            for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
                noDOFsPatchLMCnd = 3*length(BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.CP(:,1));
                noDOFsPatchLM = noDOFsPatchLM + noDOFsPatchLMCnd;
                if iCnd == 1
                    index = 3*BSplinePatches{iPatches}.noCPs;
                else
                    index = BSplinePatches{iPatches}.weakDBC.lambda{iCnd - 1}.EFT(length(BSplinePatches{iPatches}.weakDBC.lambda{iCnd - 1}.EFT));
                end
                BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.EFT = index + 1:index + noDOFsPatchLMCnd;
            end
        end
    end
    BSplinePatches{iPatches}.noDOFs = BSplinePatches{iPatches}.noDOFs + noDOFsPatchLM;
end

%% 3. Compute the total number of DOFs for the multipatch structure
noDOFsPatches = 0;
for iPatches = 1:noPatches
    noDOFsPatches = noDOFsPatches + BSplinePatches{iPatches}.noDOFs;
end

%% 4. Create a freedom table for each patch in the multipatch geometry
for iPatches = 1:noPatches
    noDOFsPatch = BSplinePatches{iPatches}.noDOFs;
    if iPatches == 1
        BSplinePatches{iPatches}.EFTPatches = 1:noDOFsPatch;
    else
        BSplinePatches{iPatches}.EFTPatches = ...
            BSplinePatches{iPatches - 1}.EFTPatches(length(BSplinePatches{iPatches - 1}.EFTPatches)) + ...
            1:BSplinePatches{iPatches - 1}.EFTPatches(length(BSplinePatches{iPatches - 1}.EFTPatches)) + ...
            noDOFsPatch;
    end
end

%% 5. Assign a number of DOFs for each Lagrange Multipliers field corresponding to the multipatch coupling
if strcmp(propCoupling.method, 'lagrangeMultipliers')
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

%% 6. Find the global numbering of the DOFs which are associated to displacements and not to the Lagrange Multipliers
EFTFreeOfLM = [];
for iPatches = 1:noPatches
    if iPatches == 1
        noCPStart = 0;
    else
        noCPStart = BSplinePatches{iPatches - 1}.EFTPatches(end);
    end
    EFTPatchFreeofLM = noCPStart + 1:noCPStart + 3*BSplinePatches{iPatches}.noCPs;
    EFTFreeOfLM = [EFTFreeOfLM EFTPatchFreeofLM];
end

%% 7. Create a freedom table for each Lagrange Multipliers field corresponding to each patch connection
if strcmp(propCoupling.method, 'lagrangeMultipliers')
    for iConnections = 1:connections.No
        if iConnections == 1
            index = noDOFsPatches;
            connections.lambda{iConnections}.EFTLagrangeMultipliers = ...
                index + 1:index + 3*connections.lambda{iConnections}.noCPs;
        else
            if isLagrangeMultipliersRotationsEnabled
                index = connections.mu{iConnections-1}.EFTLagrangeMultipliers(length(connections.mu{iConnections-1}.EFTLagrangeMultipliers));
            else
                index = connections.lambda{iConnections - 1}.EFTLagrangeMultipliers(length(connections.lambda{iConnections-1}.EFTLagrangeMultipliers));
            end
            connections.lambda{iConnections}.EFTLagrangeMultipliers = index + 1:index+3*connections.lambda{iConnections}.noCPs;
        end
        if isLagrangeMultipliersRotationsEnabled
            index = connections.lambda{iConnections}.EFTLagrangeMultipliers(length(connections.lambda{iConnections}.EFTLagrangeMultipliers));
            connections.mu{iConnections}.EFTLagrangeMultipliers = index + 1:index + 2*connections.mu{iConnections}.noCPs;
        end
    end
end

%% 8. Compute the total number of DOFs for the coupled multipatch system
if strcmp(propCoupling.method, 'lagrangeMultipliers')
    noDOFs = noDOFsPatches + noDOFsLagrangeMultipliers;
else
    noDOFs = noDOFsPatches ;
end

%% 9. Loop over all patches and change parameters corresponding to the form-finding analysis
for iPatches = 1:noPatches
    %% 9i. Change the patch parameters and save the original values
    parameters{iPatches} = BSplinePatches{iPatches}.parameters;
    BSplinePatches{iPatches}.parameters.E = EYoungFF;
    BSplinePatches{iPatches}.parameters.nue = poissonRatioFF;
    BSplinePatches{iPatches}.parameters.rho = densityFF;
    
    %% 9ii. Change the load parameters for each patch
    if isfield(BSplinePatches{iPatches}.NBC, 'noCnd')
        for iCnd = 1:BSplinePatches{iPatches}.NBC.noCnd
            if isfield(BSplinePatches{iPatches}.NBC, 'isFollower')
                NBC{iPatches} = BSplinePatches{iPatches}.NBC;
                BSplinePatches{iPatches}.NBC.isFollower(iCnd, 1) = isFollowerFF;
            else
                error('Neumann condition %d of patch %d does not contain variable isFollower', ...
                    iCnd, iPatches);
            end
        end
    else
        error('Patch %d does not contain Neumann boundary conditions', iPatches);
    end
    
    %% 9iii. Change the parameters of the cable elements
    if BSplinePatches{iPatches}.cables.No > 0
        parametersCablesSaved{iPatches}.parameters = ...
            BSplinePatches{iPatches}.cables.parameters;
    else
        parametersCablesSaved{iPatches}.parameters = [];
    end
    for iCables = 1:BSplinePatches{iPatches}.cables.No
        BSplinePatches{iPatches}.cables.parameters{iCables}.E = 0.0;
        BSplinePatches{iPatches}.cables.parameters{iCables}.areaCS = ...
            parametersCablesSaved{iPatches}.parameters{iCables}.areaCS;
        BSplinePatches{iPatches}.cables.parameters{iCables}.rho = 0.0;
        BSplinePatches{iPatches}.cables.parameters{iCables}.prestress = ...
            parametersCablesSaved{iPatches}.parameters{iCables}.prestress;
    end
end

%% 10. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
sizePrevious = 0;
homDOFs = [];
inhomDOFs = [];
valuesInhomDOFs = [];
masterDOFs = [];
slaveDOFs = [];
for iPatches = 1:noPatches
    if iPatches ~= 1
        % Determine the number of the DOFs of the previous patches
        sizePrevious = sizePrevious + BSplinePatches{iPatches - 1}.noDOFs;
        
        % Add the numbering DOFs where homogeneous Dirichlet boundary
        % conditions are applied from the patch level
        homDOFsPatch = sizePrevious + BSplinePatches{iPatches}.homDOFs;
        homDOFs = mergesorted(homDOFs, homDOFsPatch);
        
        % Add the numbering DOFs where inhomogeneous Dirichlet boundary
        % conditions are applied from the patch level
        inhomDOFsPatch = sizePrevious + BSplinePatches{iPatches}.inhomDOFs;
        inhomDOFs = mergesorted(inhomDOFs, inhomDOFsPatch);
        
        % Add the prescribed values to the DOFs where inhomogeneous 
        % Dirichlet boundary conditions are applied from the patch level
%         valuesInhomDOFsPatch = sizePrevious + BSplinePatches{counterPatches}.valuesInhomDOFs;
%         valuesInhomDOFs = mergesorted(valuesInhomDOFs,valuesInhomDOFsPatch);
        valuesInhomDOFsPatch = BSplinePatches{iPatches}.valuesInhomDOFs;
        if ~isempty(valuesInhomDOFsPatch)
            valuesInhomDOFs = horzcat(valuesInhomDOFs, valuesInhomDOFsPatch);
        end
        
        % Add the numbering of the master DOFs in the patch level
        masterDOFsPatch = sizePrevious + BSplinePatches{iPatches}.masterDOFs;
        masterDOFs = mergesorted(masterDOFs, masterDOFsPatch);
        
        % Add the numbering of the master DOFs in the patch level
        slaveDOFsPatch = sizePrevious + BSplinePatches{iPatches}.slaveDOFs;
        slaveDOFs = mergesorted(slaveDOFs, slaveDOFsPatch);
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
freeDOFs(ismember(freeDOFs, homDOFs)) = [];
freeDOFs(ismember(freeDOFs, inhomDOFs)) = [];

%% 11. Find the master/slave and the domain DOFs
fixedDOFs = mergesorted(homDOFs, inhomDOFs);
if strcmp(propCoupling.method, 'mortar')
    [BSplinePatches,connections] = findDOFsDDMMortarIGAMembrane ...
        (BSplinePatches, connections, propCoupling, fixedDOFs);
end

%% 12. Create a DOF numbering for each patch and each Lagrange Multipliers field

% For the patches :
% _________________

for iPatches = 1:noPatches
    % Get the number of Control Points in xi-direction
    nxi = length(BSplinePatches{iPatches}.CP(:, 1, 1));
    
    % Get the number of Control Points in eta-direction
    neta = length(BSplinePatches{iPatches}.CP(1, :, 1));
    
    % Initialize the DOF numbering array
    BSplinePatches{iPatches}.DOFNumbering = zeros(nxi, neta, 3);
    
    % Compute the entries of the DOF numbering array
    k = 1;
    for cpj = 1:neta
        for cpi = 1:nxi
            BSplinePatches{iPatches}.DOFNumbering(cpi, cpj, 1) = k;
            BSplinePatches{iPatches}.DOFNumbering(cpi, cpj, 2) = k + 1;
            BSplinePatches{iPatches}.DOFNumbering(cpi, cpj, 3) = k + 2;

            % Update counter
            k = k + 3;
        end
    end
    
    % Create a DOF numbering for each Lagrange Multipliers field within the
    % patch
    if isWeakDBC(iPatches, 1)
        if strcmp(BSplinePatches{iPatches}.weakDBC.method, 'lagrangeMultipliers')
            for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
                % Get the number of Control Points
                noXiLambda = length(BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.Xi);

                % Initialize the field of the DOF numbering
                BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.DOFNumbering = zeros(noXiLambda, 3);

                % Compute the entries of the DOF numbering array
                k = 1;
                for cpi = 1:noXiLambda
                    BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 1) = k;
                    BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 2) = k + 1;
                    BSplinePatches{iPatches}.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 3) = k + 2;

                    % Update counter
                    k = k + 3;
                end
            end
        end
    end
end

% For the Lagrange Multipliers field :
% ____________________________________

if strcmp(propCoupling.method, 'lagrangeMultipliers')
    for iConnections = 1:connections.No
        % Get the number of Control Points in xi-direction
        if isLagrangeMultipliersDisplacementsEnabled
            noXiLambda = length(connections.lambda{iConnections}.CP(:, 1));
        end
        if isLagrangeMultipliersRotationsEnabled
            nxiMu = length(connections.mu{iConnections}.CP(:, 1));
        end

        % Lagrange multipliers for the traction forces :
        % ______________________________________________

        if isLagrangeMultipliersDisplacementsEnabled
            % Initialize the field
            connections.lambda{iConnections}.DOFNumbering = zeros(noXiLambda, 3);

            % Compute the entries of the DOF numbering array
            k = 1;
            for cpi = 1:noXiLambda
                connections.lambda{iConnections}.DOFNumbering(cpi, 1) = k;
                connections.lambda{iConnections}.DOFNumbering(cpi, 2) = k + 1;
                connections.lambda{iConnections}.DOFNumbering(cpi, 3) = k + 2;

                % Update counter
                k = k + 3;
            end
        end
        % Lagrange multipliers for the traction moments :
        % _______________________________________________

        if isLagrangeMultipliersRotationsEnabled
            % Initialize the field
            connections.mu{iConnections}.DOFNumbering = zeros(nxiMu, 2);

            % Compute the entries of the DOF numbering array
            k = 1;
            for cpi = 1:nxiMu
                connections.mu{iConnections}.DOFNumbering(cpi, 1) = k;
                connections.mu{iConnections}.DOFNumbering(cpi, 2) = k + 1;

                % Update counter
                k = k + 2;
            end
        end
    end
end

%% 13. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions
for iPatches = 1:noPatches
    if isWeakDBC(iPatches,1)
        isConstMtxWeakDBC = true;
        if isfield(BSplinePatches{iPatches}.weakDBC, 'method')
            if strcmp(BSplinePatches{iPatches}.weakDBC.method, 'penalty') || ...
                    (strcmp(BSplinePatches{iPatches}.weakDBC.method, 'nitsche') && ...
                    ~BSplinePatches{iPatches}.weakDBC.estimationStabilPrm)
                computeWeakDBCConstantProblemMatrices = @computeWeakDBCMtxPenaltyIGAMembrane;
            elseif strcmp(BSplinePatches{iPatches}.weakDBC.method, 'lagrangeMultipliers')
                computeWeakDBCConstantProblemMatrices = @computeWeakDBCMtxLagrangeMultipliersIGAMembrane;
            elseif strcmp(BSplinePatches{iPatches}.weakDBC.method, 'nitsche') && ...
                    BSplinePatches{iPatches}.weakDBC.estimationStabilPrm
                isConstMtxWeakDBC = false;
            else
                error('Define a valid method for the application of weak Dirichlet boundary conditions in BSplinePatches{iPatches}.weakDBC.method');
            end
            if isConstMtxWeakDBC
                BSplinePatches{iPatches}.KConstant = ...
                    computeWeakDBCConstantProblemMatrices ...
                    (BSplinePatches{iPatches}, connections,...
                    BSplinePatches{iPatches}.noDOFs, propCoupling);
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

%% 14. Compute the constant problem matrices according to the chosen method
if strcmp(propCoupling.method, 'penalty') || ...
        (strcmp(propCoupling.method, 'nitsche') && ~propCoupling.estimationStabilPrm) || ...
        (~strcmp(propCoupling.method, 'lagrangeMultipliers') && isfield(propCoupling, 'alphaR'))
    if (strcmp(propCoupling.method, 'nitsche') && ~propCoupling.estimationStabilPrm)
        if isfield(propCoupling, 'alphaD')
            error('For automatic estimation of the stabilization parameter corresponding to the multipatch coupling no displacement penalization needs to be specified in couplingMethod.alphaD');
        end
    end
    computeConstantProblemMatrices = @computeConstantMtxForDDMPenaltyIGAThinStructure;
elseif strcmp(propCoupling.method, 'lagrangeMultipliers')
    computeConstantProblemMatrices = @computeConstantMtxForDDMLagrangeMultipliersIGAThinStructure;
elseif strcmp(propCoupling.method, 'mortar')
    computeConstantProblemMatrices = @computeConstantMtxForDDMMortarIGAThinStructure;
else
    computeConstantProblemMatrices = 'undefined';
end
if isa(computeConstantProblemMatrices, 'function_handle')
    KConstant = computeConstantProblemMatrices(BSplinePatches,connections,noDOFs,propCoupling);
else
    KConstant = 'undefined';
end

%% 15. Initialize the displacement field
dHatPrevious = zeros(noDOFs, 1);

%% 16. Loop over all the form-finding iterations
if strcmp(outMsg, 'outputEnabled')
    msgPNR = sprintf(strcat(tab, '\tLooping over the form finding iterations\n',...
        tab,'\t----------------------------------------\n\n'));
    fprintf(msgPNR);
end
while ~isConverged && counterFoFiIter <= propFormFinding.maxNoIter
    %% 16i. Solve the static linear problem
    [dHat, ~, rH, ~, ~, ~, ~, ~, BSplinePatches, propCoupling, ~] = ...
        solve_IGANLinearSystem ...
        (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
        connections, zeros(noDOFs,1), dHatDot, dHatDDot, KConstant, ...
        massMtx, dampMtx, computeTanStiffMtxAndResVct, ...
        @computeUpdatedGeometryIGAThinStructureMultipatches, freeDOFs, ...
        homDOFs, inhomDOFs, valuesInhomDOFs, updateDirichletBCs, masterDOFs, ...
        slaveDOFs, solve_LinearSystem, t, propCoupling, propTransientAnalysis, ...
        propNLinearAnalysis, propIDBC, plot_IGANLinear, isReferenceUpdated, ...
        isCosimulationWithEmpire, strcat('\t', tab), propGraph, '');
    if length(find(rH)) > 2
        warning('More than 1 iterations needed for convergence');
    end
    
    %% 16ii. Check if the solution is not NaN
    if sum(isnan(dHat))
        if strcmp(outMsg, 'outputEnabled')
            warning('NaN solution has been obtained');
        end
        break;
    end
    
    %% Debugging
%     graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
%         (BSplinePatches,dHat,graph,'outputEnabled');

    %% 16iii. Update the Control Point coordinates of the B-Spline patch
    for iPatches = 1:noPatches
        BSplinePatches{iPatches}.CP = BSplinePatches{iPatches}.CPd;
        CPHistoryMultipatch{iPatches}.CPHistory{counterFoFiIter + 1} = ...
            BSplinePatches{iPatches}.CP;
    end
    
    %% 16iv. Compute the difference between the displacements from both iteration steps
	delta_dHat = dHat - dHatPrevious;
    
    %% 16v. Compute the residual and check convergence
    resHistory(counterFoFiIter,1) = norm(delta_dHat(EFTFreeOfLM));
    if strcmp(outMsg, 'outputEnabled')
        msgNR = sprintf(strcat(tab,'\t||delta_dHat|| = %d at form-finding iteration No. %d \n'), ...
            resHistory(counterFoFiIter, 1), counterFoFiIter);
        fprintf(msgNR);
    end
    if counterFoFiIter > 1
        if resHistory(counterFoFiIter, 1) > resHistory(counterFoFiIter - 1, 1)
            if strcmp(outMsg, 'outputEnabled')
                warning('The error in the form-finding iterations is increasing\n');
            end
        end
    end
    if resHistory(counterFoFiIter, 1) < propFormFinding.tolerance
        if counterFoFiIter > propFormFinding.minNoIter
            if strcmp(outMsg, 'outputEnabled')
                fprintf(strcat(tab,' \tForm-finding iterations converged!\n\n'));
            end
            isConverged = true;
            break;
        end
    end

    %% 16vi. Save the displacement field of the previous form-finding step
    dHatPrevious = dHat;
    
    %% 16vii. Update form-finding iteration counter and the pseudotime
    counterFoFiIter = counterFoFiIter + 1;
    t = t + propTransientAnalysis.dt;
end

%% 17. Re-assign arrays
noIter = counterFoFiIter - 1;
if noIter == propFormFinding.maxNoIter
    if strcmp(outMsg, 'outputEnabled')
        warning(strcat(tab, ' \tForm-finding iterations did not converge up to tolerance %f!\n\n'), ...
            propFormFinding.tolerance);
    end
end
for iPatches = 1:noPatches
    BSplinePatches{iPatches} = rmfield(BSplinePatches{iPatches}, 'parameters');
    BSplinePatches{iPatches}.parameters = parameters{iPatches};
    BSplinePatches{iPatches} = rmfield(BSplinePatches{iPatches}, 'NBC');
    BSplinePatches{iPatches}.NBC = NBC{iPatches};
    for iCables = 1:BSplinePatches{iPatches}.cables.No
        BSplinePatches{iPatches}.cables.parameters{iCables}.E = parametersCablesSaved{iPatches}.parameters{iCables}.E;
        BSplinePatches{iPatches}.cables.parameters{iCables}.areaCS = parametersCablesSaved{iPatches}.parameters{iCables}.areaCS;
        BSplinePatches{iPatches}.cables.parameters{iCables}.rho = parametersCablesSaved{iPatches}.parameters{iCables}.rho;
        BSplinePatches{iPatches}.cables.parameters{iCables}.prestress = parametersCablesSaved{iPatches}.parameters{iCables}.prestress;
    end
end

%% 18. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Form finding analysis took %.2d seconds \n\n', computationalTime);
    fprintf('_________________Form Finding Analysis Ended__________________\n');
    fprintf('##############################################################\n\n\n');
end

end
