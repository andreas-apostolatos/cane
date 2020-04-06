function [dHat, CPHistory, resHistory, hasConverged, FComplete, rankD, ...
    condK, minEig, BSplinePatches, propCoupling, minElAreaSize] = ...
    solve_DDMIGAMembraneMultipatchesNLinear ...
    (BSplinePatches, connections, propCoupling, propNLinearAnalysis, ...
    solve_LinearSystem, plot_IGANLinear, graph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement field corresponding to the geometrically
% nonlinear analysis over the multipatch isogeometric membrane analysis for
% the chosen coupling method.
%
%                Input :
%       BSplinePatches : Structure containing all the information regarding 
%                        the connections between the multipatches
%          connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                 ...      ...    ...   ...  ...   ...]
%         propCoupling : Properties of the multipatch coupling
%                           .alphaD : penalty factor for the displacement
%                                     coupling
%                           .alphaR : penalty factor for the rotation
%                                     coupling
%                             .intC : On the integration of the coupling
%                                     interface 
% propNLinearAnalysis : Structure on the non-linear analysis settings :
%                                    #Load Steps/Stages#
%                              .noLoadSteps : Selected number of load steps 
%                                             (stages of the loading)
%                                .tolerance : Tolerance for the Newton 
%                                             iterations on the 2-norm
%                                  .maxIter : Maximum number of the Newton 
%                                             iteration
%   solve_LinearSystem : Function handle to the linear equation system
%                        solver
%      plot_IGANLinear : Function handle to the plotting of the current
%                        configuration throughout the nonlinear iterations
%                graph : On the graphics
%               outMsg : Whether or not to output message on refinement 
%                        progress
%                       'outputEnabled' : enables output information
%   
%               Output :
%                 dHat : Array containing the displacement field of each 
%                        patch in the coupled system
%            CpHistory : Array containing the deformation history of the
%                        patches in the coupled system
%           resHistory : The history of the residual force magnitude
%                        throughout the non-linea analysis and the Newton 
%                        iterations
%           minElASize : Array containing the minimum element area size for
%                        each patch
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
% 6. Create a freedom table for each Lagrange Multipliers field corresponding to each patch connection
%
% 7. Compute the total number of DOFs for the coupled multipatch system
%
% 8. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
%
% 9. Find the mortar (master/slave) DOFs
%
% 10. Create a DOF numbering for each patch and each Lagrange Multipliers field
%
% 11. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions
%
% 12. Compute the constant problem matrices according to the chosen method
%
% 13. Solve the geometrically nonlinear problem
%
% 14. Appendix
%
%% Function main body
if ~strcmp(propCoupling.method,'penalty') && ~strcmp(propCoupling.method,'lagrangeMultipliers') && ...
        ~strcmp(propCoupling.method,'mortar') && ~strcmp(propCoupling.method,'nitsche')
    error('Undefined coupling method for the multipatch structure in propCoupling.method was selected');
end
if strcmp(outMsg,'outputEnabled')
    fprintf('______________________________________________________________\n');
    fprintf('##############################################################\n');
    fprintf('Geometrically nonlinear analysis for the domain decomposition\n');
    fprintf('of the isogeometric membrane using the ');
    if strcmp(propCoupling.method,'lagrangeMultipliers')
        if isfield(propCoupling,'alphaD')
            if propCoupling.alphaD ~= 0
                fprintf('augmented Lagrange Multipliers\n');
                fprintf('method has been initiated\n\n');
            else
                fprintf('Lagrange Multipliers method\n');
                fprintf('has been initiated\n\n');
            end
        end
    elseif strcmp(propCoupling.method,'penalty')
        fprintf('Penalty method has been \n');
        fprintf('initiated\n\n');
	elseif strcmp(propCoupling.method,'mortar')
        fprintf('mortar method has been \n');
        fprintf('initiated\n\n');
    elseif strcmp(propCoupling.method,'nitsche')
        fprintf('Nitsche method has been \n');
        fprintf('initiated\n\n');
    else
        fprintf('\n');
        error('Choose coupling method in propCoupling.method');
    end
    fprintf('\n\n');
    fprintf('Coupling properties \n');
    fprintf('------------------- \n\n');
    if strcmp(propCoupling.method,'penalty')
        if ~isfield(propCoupling,'alphaD')
            error('Choose penalty factors for the displacement field in propCoupling.alphaD');
        end
        if size(propCoupling.alphaD) ~= connections.No
            error('The number of the penalty factors in propCoupling.alphaD must be the same as for the number of connections connections.No');
        end
        for iPatches = 1:connections.No
            fprintf('Displacement penalty factor for patch %d equal to alpha = %d \n',iPatches,propCoupling.alphaD(iPatches,1));
        end
    elseif strcmp(propCoupling.method,'lagrangeMultipliers')
        if ~isfield(propCoupling,'alphaD')
            error('Choose penalty factors for the displacement field in propCoupling.alphaD');
        end
        if size(propCoupling.alphaD) ~= connections.No
            error('The number of the penalty factors in propCoupling.alphaD must be the same as for the number of connections connections.No');
        end
        for iPatches = 1:connections.No
            fprintf('Displacement penalty factor for patch %d equal to alpha = %d \n',iPatches,propCoupling.alphaD(iPatches,1));
        end
        if ~isfield(connections,'lambda')
            fprintf('\n');
            fprintf('Choose a discretization for the Lagrange Multipliers field in connections.lambda');
        end
    elseif strcmp(propCoupling.method,'mortar')
        if ~isfield(propCoupling,'isSlaveSideCoarser')
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
            fprintf('The function handle to the computation of the re-arranged system is \n%s\n',func2str(propCoupling.computeRearrangedProblemMtrcs));
        end
    elseif strcmp(propCoupling.method,'nitsche')
        if propCoupling.estimationStabilPrm == true
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

% Flag on whether the reference configuration is updated
isReferenceUpdated = false;

% Check input
for iPatches = 1:noPatches
    if isfield(BSplinePatches{iPatches}.weakDBC,'computeTangMtxResVct') && ...
            (strcmp(BSplinePatches{iPatches}.weakDBC.method,'penalty') || ...
            strcmp(BSplinePatches{iPatches}.weakDBC.method,'lagrangeMultipliers'))
        error('No field weakDBC.%s needs to be defined for the method %s corresponding to the application of weak Dirichlet boundary conditions', ...
            BSplinePatches{iPatches}.weakDBC.computeTangMtxResVct,BSplinePatches{iPatches}.weakDBC.method);
    end
end

% Check if Lagrange Multipliers field for the displacement coupling is
% enforced
isLagrangeMultipliersDisplacementsEnabled = false;
if (isfield(connections,'lambda') && (strcmp(propCoupling.method,'lagrangeMultipliers')) || strcmp(propCoupling.method,'mortar'))
    isLagrangeMultipliersDisplacementsEnabled = true;
    if isfield(connections,'lambda')
        if isempty(connections.lambda) && strcmp(propCoupling.method,'lagrangeMultipliers')
            error('Lagrange multipliers discretization for the displacement coupling exists but its empty');
        end
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

% Set the co-simulation flag with EMPIRE to false
isCosimulationWithEmpire = false;

% Pseudo-transient analysis only for storing the coupling data
propTransientAnalysis.timeDependence = 'steadyState';
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;
t = 'undefined';

%% 1. Assign the computation of the tangent stiffness matrix according to the selected method
if strcmp(propCoupling.method,'penalty') || strcmp(propCoupling.method,'lagrangeMultipliers') || ...
        strcmp(propCoupling.method,'mortar')
    computeTanStiffMtxAndResVct = @computeTangentStiffMtxResVctDDMPenaltyIGAMembrane;
elseif strcmp(propCoupling.method,'nitsche')
    computeTanStiffMtxAndResVct = @computeTangentStiffMtxResVctDDMNitscheIGAMembrane;
else
    error('Choose coupling method in propCoupling.method');
end

%% 2. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
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
            BSplinePatches{iPatches-1}.EFTPatches(length(BSplinePatches{iPatches-1}.EFTPatches)) + ...
            1:BSplinePatches{iPatches-1}.EFTPatches(length(BSplinePatches{iPatches-1}.EFTPatches)) + ...
            noDOFsPatch;
    end
end

%% 5. Assign a number of DOFs for each Lagrange Multipliers field corresponding to the multipatch coupling
if strcmp(propCoupling.method,'lagrangeMultipliers')
    noDOFsLagrangeMultipliers = 0;
    for iConnOuter = 1:connections.No
        noDOFsLagrangeMultipliers = noDOFsLagrangeMultipliers + ...
            3*connections.lambda{iConnOuter}.noCPs;
        if isLagrangeMultipliersRotationsEnabled
            noDOFsLagrangeMultipliers = noDOFsLagrangeMultipliers + ...
            2*connections.mu{iConnOuter}.noCPs;
        end
    end
end

%% 6. Create a freedom table for each Lagrange Multipliers field corresponding to each patch connection
if strcmp(propCoupling.method,'lagrangeMultipliers')
    for iConnOuter = 1:connections.No
        if iConnOuter == 1
            index = noDOFsPatches;
            connections.lambda{iConnOuter}.EFTLagrangeMultipliers = ...
                index + 1:index + 3*connections.lambda{iConnOuter}.noCPs;
        else
            if isLagrangeMultipliersRotationsEnabled
                index = connections.mu{iConnOuter-1}.EFTLagrangeMultipliers(length(connections.mu{iConnOuter-1}.EFTLagrangeMultipliers));
            else
                index = connections.lambda{iConnOuter-1}.EFTLagrangeMultipliers(length(connections.lambda{iConnOuter-1}.EFTLagrangeMultipliers));
            end
            connections.lambda{iConnOuter}.EFTLagrangeMultipliers = index + 1:index+3*connections.lambda{iConnOuter}.noCPs;
        end
        if isLagrangeMultipliersRotationsEnabled
            index = connections.lambda{iConnOuter}.EFTLagrangeMultipliers(length(connections.lambda{iConnOuter}.EFTLagrangeMultipliers));
            connections.mu{iConnOuter}.EFTLagrangeMultipliers = index + 1:index + 2*connections.mu{iConnOuter}.noCPs;
        end
    end
end

%% 7. Compute the total number of DOFs for the coupled multipatch system
if strcmp(propCoupling.method,'lagrangeMultipliers')
    noDOFs = noDOFsPatches + noDOFsLagrangeMultipliers;
else
    noDOFs = noDOFsPatches ;
end
dHat = zeros(noDOFs,1);

%% 8. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
sizePrevious = 0;
homDOFs = [];
inhomDOFs = [];
valuesInhomDOFs = [];
masterDOFs = [];
slaveDOFs = [];
for iPatches = 1:noPatches
    if iPatches ~= 1
        % Determine the number of the DOFs of the previous patches
        sizePrevious = sizePrevious + BSplinePatches{iPatches-1}.noDOFs;
        
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
if strcmp(propCoupling.method,'lagrangeMultipliers')
    freeDOFs = 1:noDOFs;
else
    freeDOFs = 1:noDOFsPatches;
end
freeDOFs(ismember(freeDOFs,homDOFs)) = [];
freeDOFs(ismember(freeDOFs,inhomDOFs)) = [];

%% 9. Find the mortar (master/slave) DOFs
fixedDOFs = mergesorted(homDOFs,inhomDOFs);
if strcmp(propCoupling.method,'mortar')
    [BSplinePatches,connections] = findDOFsDDMMortarIGAMembrane...
        (BSplinePatches,connections,propCoupling,fixedDOFs);
end

%% 10. Create a DOF numbering for each patch and each Lagrange Multipliers field

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


for iConnOuter = 1:connections.No
    % Get the number of Control Points in xi-direction
    if isLagrangeMultipliersDisplacementsEnabled
        noXiLambda = length(connections.lambda{iConnOuter}.CP(:,1));
    end
    if isLagrangeMultipliersRotationsEnabled
        nxiMu = length(connections.mu{iConnOuter}.CP(:,1));
    end
    
    % Lagrange multipliers for the traction forces :
    % ______________________________________________
    
    if isLagrangeMultipliersDisplacementsEnabled
        % Initialize the field
        connections.lambda{iConnOuter}.DOFNumbering = zeros(noXiLambda,3);

        % Compute the entries of the DOF numbering array
        k = 1;
        for cpi = 1:noXiLambda
            connections.lambda{iConnOuter}.DOFNumbering(cpi,1) = k;
            connections.lambda{iConnOuter}.DOFNumbering(cpi,2) = k + 1;
            connections.lambda{iConnOuter}.DOFNumbering(cpi,3) = k + 2;

            % Update counter
            k = k + 3;
        end
    end
    % Lagrange multipliers for the traction moments :
    % _______________________________________________
    
    if isLagrangeMultipliersRotationsEnabled
        % Initialize the field
        connections.mu{iConnOuter}.DOFNumbering = zeros(nxiMu,2);

        % Compute the entries of the DOF numbering array
        k = 1;
        for cpi = 1:nxiMu
            connections.mu{iConnOuter}.DOFNumbering(cpi,1) = k;
            connections.mu{iConnOuter}.DOFNumbering(cpi,2) = k + 1;

            % Update counter
            k = k + 2;
        end
    end
end

%% 11. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions
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

%% 12. Compute the constant problem matrices according to the chosen method
if strcmp(propCoupling.method,'penalty') || (strcmp(propCoupling.method,'nitsche') && ~propCoupling.estimationStabilPrm)
    computeConstantProblemMatrices = @computeConstantMtxForDDMPenaltyIGAThinStructure;
elseif strcmp(propCoupling.method,'lagrangeMultipliers')
    computeConstantProblemMatrices = @computeConstantMtxForDDMLagrangeMultipliersIGAThinStructure;
elseif strcmp(propCoupling.method,'mortar')    
    computeConstantProblemMatrices = @computeConstantMtxForDDMMortarIGAThinStructure;
else
    computeConstantProblemMatrices = 'undefined';
end
if isa(computeConstantProblemMatrices,'function_handle')
    KConstant = computeConstantProblemMatrices(BSplinePatches,connections,noDOFs,propCoupling);
else
    KConstant = 'undefined';
end

%% 13. Solve the geometrically nonlinear problem
[dHat, CPHistory, resHistory, hasConverged, FComplete, rankD, condK, ...
    minEig, BSplinePatches, propCoupling, minElAreaSize] = ...
    solve_IGANLinearSystem ...
    (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
    connections, dHat, dHatDot, dHatDDot, KConstant, massMtx, dampMtx, ...
    computeTanStiffMtxAndResVct, ...
    @computeUpdatedGeometryIGAThinStructureMultipatches, freeDOFs, ...
    homDOFs, inhomDOFs, valuesInhomDOFs, updateDirichletBCs, masterDOFs, ...
    slaveDOFs, solve_LinearSystem, t, propCoupling, ...
    propTransientAnalysis, propNLinearAnalysis, propIDBC, ...
    plot_IGANLinear, isReferenceUpdated, isCosimulationWithEmpire, ...
    strcat('\t',tab), graph, outMsg);

%% 14. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Static geometrically nonlinear analysis took %.2d seconds \n\n', computationalTime);
    fprintf('_____________Nonlinear Analysis Ended________________\n');
    fprintf('#####################################################\n\n\n');
    fprintf('\n');
end

end
