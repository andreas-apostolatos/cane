function [dHatHistory,resHistory,BSplinePatches,propCoupling,minElAreaSize] = ...
    solve_DDMIGAMembraneMultipatchesNLinearTransient...
    (BSplinePatches,connections,computeInitCnds,propCoupling,...
    propNLinearAnalysis,propStrDynamics,propPostproc,solve_LinearSystem,...
    propOutput,pathToOutput,caseName,propEmpireCoSimulation,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the history of the displacement field corresponding to the
% solution of the nonlinear transient isogeometric membrane problem
% modelled with multiple patches coupled with either the Penalty, the 
% Lagrange Multipliers, or the Nitsche method.
%
%                    Input :
%           BSplinePatches : Structure containing all the information 
%                            regarding the connections between the 
%                            multipatches
%              connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%                        .lambda : Cell array containing the discretization
%                                  of the Lagrange Multipliers fields
%                                  corresponding to the traction forces
%                            .mu : Cell array containing the discretization
%                                  of the Lagrange Multipliers fields
%                                  corresponding to the traction moments
%          computeInitCnds : Function handle to the computation of initial 
%                            conditions
%             propCoupling : Properties of the multipatch coupling
%                              .method : The coupling method
%                .estimationStabilPrm : Flag for the automatic estimation 
%                                       of the stabilization parameter 
%                                       when Nitsche is chosen as a 
%                                       coupling method
%                             .alphaD : penalty factor for the displacement
%                                       coupling
%                             .alphaR : penalty factor for the rotation
%                                       coupling
%                               .intC : On the integration of the coupling
%                                       interface 
%      propNLinearAnalysis : Structure on the non-linear analysis settings:
%                              .noLoadSteps : Selected number of load steps 
%                                             (stages of the loading)
%                                .tolerance : Tolerance for the Newton 
%                                             iterations on the 2-norm
%                                  .maxIter : Maximum number of the Newton 
%                                             iteration
%          propStrDynamics : On the properties of the structural dynamics 
%                            analysis:
%                           .method : Time integration method
%                        .isDamping : Flag on whether damping is enabled
%                           .TStart : Starting time of the simulation
%                             .TEnd : End time of the simulation
%                      .noTimeSteps : Number of time steps
%                               .dt : Time step size
%     .computeProblemMtrcsTransient : Function handle to the computation of 
%                                     the problem matrices according to the 
%                                     defined time integration scheme given 
%                                     the system matrix of the steady-state
%                                     system as well as the mass matrix of 
%                                     the problem
%                .computeUpdatedVct : Function handle to the computation of 
%                                     the updated values of the discretized
%             propPostproc : On the properties of the postprocessing
%                               .resultant : Array of strings containing 
%                                            the name of each resultant to 
%                                            be written out
%                        .computeResultant : Array of strings representing 
%                                            the function handle for the 
%                                            computation of the desirable 
%                                            resultant
%       solve_LinearSystem : Function handle to the linear equation system
%                            solver
%              propOutput : Properties for writting out the results
%                              .writeOutput : Function handle to outputting 
%                                             the results in a folder
%                           .writeFrequency : Frequency for writting out
%                                             the results
%                            .saveFrequency : Frequency for saving temorary
%                                             results at the time step
%                                             level
%             pathToOutput : Path to the folder where the results are 
%                            written into
%                 caseName : The name of the case in the inputGiD case 
%                            folder
%   propEmpireCoSimulation : Properties for the co-simulation with EMPIRE
%                           .isCoSimulation : Flag on whether co-simulation 
%                                             with EMPIRE is assumed
%                         .isInterfaceLayer : Flag on whether the matlab
%                                             client is used as an
%                                             interface layer
%             .isGaussPointProvidedDirichlet : Flag on whether the Gauss
%                                              Point data along the
%                                              Dirichlet boundaries are
%                                              provided
%             .isGaussPointProvidedInterface : Flag on whether the Gauss
%                                              Point data along the
%                                              interface boundaries are
%                                              provided
%                          .propIntDirichlet : Integration along the
%                                              Dirichlet boundaries
%                                              properties,
%                                                .type : 'default', 'user'
%                                               .noGPs : Number of Gauss
%                                                        Points
%                          .propIntInterface : Integration along the
%                                              interface boundaries
%                                              properties,
%                                                .type : 'default', 'user'
%                                               .noGPs : Number of Gauss
%                                                        Points
%                              .strMatlabXml : Name of the xml file for 
%                                              connection to Empire
%                   outMsg : Whether or not to output message on refinement 
%                            progress
%                              'outputEnabled' : enables output information
%   
%               Output :
%                 dHat : Array containing the displacement field of each 
%                        patch in the coupled system
%           resHistory : The history of the residual force magnitude
%                        throughout the non-linea analysis and the Newton 
%                        iterations
%       BSplinePatches : The updated B-Spline patches array
%         propCoupling : The updated structure of the coupling properties
%           minElASize : Array containing the minimum element area size for
%                        each patch
%
% Function Layout :
%
% 0. Read input
%
% 1. Assign the function handle to the computation of the constant matrix of the system corresponding to the multipatch coupling
%
% 2. Assign the computation of the tangent stiffness matrix according to the selected method
%
% 3. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
%
% 4. Compute the total number of DOFs for the multipatch structure
%
% 5. Create a freedom table for each patch in the multipatch geometry
%
% 6. Assign a number of DOFs for each Lagrange Multipliers field corresponding to the multipatch coupling
%
% 7. Create a freedom table for each Lagrange Multipliers field corresponding to each patch connection
%
% 8. Compute the total number of DOFs for the coupled multipatch system
%
% 9. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary and master-slave conditions are applied to account for the coupled system
%
% 10. Find the mortar (master/slave) DOFs
%
% 11. Create a DOF numbering for each patch and each Lagrange Multipliers field
%
% 12. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions and update the number of DOFs per patch if it is enriched with Lagrange Multipliers
%
% 13. Solve the transient problem
%
% 14. Appendix
%
%% Function main body
if strcmp(propCoupling.method,'penalty')
    couplingMethod = 'Penalty';
elseif strcmp(propCoupling.method,'lagrangeMultipliers')
    if isfield(propCoupling,'alphaD')
        if propCoupling.alphaD ~= 0
            couplingMethod = 'augmented Lagrange Multipliers';
        else
            couplingMethod = 'Lagrange Multipliers';
        end
    else
        couplingMethod = 'Lagrange Multipliers';
    end
elseif strcmp(propCoupling.method,'mortar')
    couplingMethod = 'Mortar';   
elseif strcmp(propCoupling.method,'nitsche')
    couplingMethod = 'Nitsche';
else
    error('Wrong coupling method for the multipatch structure in propCoupling.method was selected');
end
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________\n');
    fprintf('#############################################################\n');
    fprintf('Transient geometrically nonlinear analysis for a multipatch \n');
    fprintf('isogeometric membrane structure has been initiated using the\n');
    fprintf('%s method \n\n',couplingMethod);
    fprintf('Coupling properties\n');
    fprintf('-------------------\n\n');
    if strcmp(couplingMethod,'Penalty') || strcmp(couplingMethod,'augmented Lagrange Multipliers')
        if ~isfield(propCoupling,'alphaD')
            error('Choose penalty factors for the displacement field in propCoupling.alphaD');
        end
        if size(propCoupling.alphaD) ~= connections.No
            error('The number of the penalty factors in propCoupling.alphaD must be the same as for the number of connections connections.No');
        end
        for iPatches = 1:connections.No
            fprintf('Displacement penalty factor for patch %d equal to alpha = %d \n',iPatches,propCoupling.alphaD(iPatches,1));
        end
        if isfield(propCoupling,'alphaR')
            fprintf('Penalty factor for the rotations equals %d\n',propCoupling.alphaR);
        end
    elseif strcmp(propCoupling.method,'lagrangeMultipliers')
        for iConnections = 1:length(connections.lambda)
            if isfield(connections,'lambda')
                fprintf('Number of elements for lambda of connection %d: %d\n',iConnections,length(unique(connections.lambda{iConnections}.Xi)-1));
            else
                error('Define Lagrange Multipliers discretization for the patch coupling\n');
            end
            if isfield(connections,'mu')
                fprintf('Number of elements for mu of connection %d: %d',iConnections,length(unique(connections.mu{iConnections}.Xi)-1));
            end
        end
        fprintf('\n');
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
    elseif strcmp(propCoupling.method,'Nitsche')
        if propCoupling.estimationStabilPrm
            fprintf('Automatic estimation of the stabilization parameter has been selected');
        else
            fprintf('Manual assignment of the stabilization parameter has been selected');
            if isfield(propCoupling,'alphaD')
                fprintf('Stabilization parameter for the displacements equals %d\n',propCoupling.alphaD);
            else
                error('Choose stabilization parameter for the displacement continuity\n');
            end 
            if isfield(propCoupling,'alphaR')
                fprintf('Penalty factor for the rotations equals %d\n',propCoupling.alphaR);
            end
            if isfield(propCoupling,'gammaTilde')
                fprintf('Linear combination factor for the interface traction field chosen equal to %d\n',propCoupling.gammaTilde);
            else
                error('Choose linear combination factor the interface traction field in propCoupling.gammaTilde');
            end
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
    if propEmpireCoSimulation.isCoSimulation
        fprintf('Co-simulation through Empire \n');
        fprintf('---------------------------- \n\n');
        if isempty(propEmpireCoSimulation.strMatlabXml)
            error('No xml file for the Matlab code is provided');
        else
            fprintf('Matlab xml file for connection to Empire: %s.xml\n',propEmpireCoSimulation.strMatlabXml);
        end
        fprintf('\n');
    end
    fprintf('_____________________________________________________________\n');
    fprintf('\n');
    
    % start measuring computational time
    tic;
end

%% 0. Read input

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Flag on whether the reference geometry is updated
isReferenceUpdated = false;

% Title for the output file
if ~strcmp(propNLinearAnalysis.method,'undefined')
    title = 'Geometrically nonlinear transient isogeometric membrane analysis';
else
    title = 'Geometrically linear transient isogeometric membrane analysis';
end

% Function handle to the computation of the mass matrix
computeMassMtx = @computeIGAMassMtxThinStructure;

% Function handle to the computation of the updated deformed geometry 
computeUpdatedGeometry = @computeUpdatedGeometryIGAThinStructureMultipatches;

% Initialize dummy variables
nodesALE = 'undefined';
computeUpdatedMesh = 'undefined';

% Tabulation for the output message
tab = '\t';

% Check if Lagrange Multipliers field for the displacement coupling is
% enforced
isLagrangeMultipliersDisplacementsEnabled = false;
if (isfield(connections,'lambda') && strcmp(propCoupling.method,'lagrangeMultipliers') ) || strcmp(propCoupling.method,'mortar')
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
if strcmp(propCoupling.method,'lagrangeMultipliers') || strcmp(propCoupling.method,'mortar')
    isLagrangeMultipliersRotationsEnabled = false;
    if isfield(connections,'mu')
        isLagrangeMultipliersRotationsEnabled = true;
        if isempty(connections.mu)
            error('Lagrange multipliers discretization for the rotational coupling exists but its empty');
        end
    end
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

%% 1. Assign the function handle to the computation of the constant matrix of the system corresponding to the multipatch coupling
if strcmp(propCoupling.method,'penalty') || ...
        (strcmp(propCoupling.method,'nitsche') && ~propCoupling.estimationStabilPrm) || ...
        (~strcmp(propCoupling.method,'lagrangeMultipliers') && isfield(propCoupling,'alphaR'))
    if (strcmp(propCoupling.method,'nitsche') && ~propCoupling.estimationStabilPrm)
        if isfield(propCoupling,'alphaD')
            error('For automatic estimation of the stabilization parameter corresponding to the multipatch coupling no displacement penalization needs to be specified in couplingMethod.alphaD');
        end
    end
    computeConstantProblemMatrices = @computeConstantMtxForDDMPenaltyIGAThinStructure;
elseif strcmp(propCoupling.method,'lagrangeMultipliers')
    computeConstantProblemMatrices = @computeConstantMtxForDDMLagrangeMultipliersIGAThinStructure;
elseif strcmp(propCoupling.method,'mortar')
    computeConstantProblemMatrices = @computeConstantMtxForDDMMortarIGAThinStructure;
else
    computeConstantProblemMatrices = 'undefined';
end

%% 2. Assign the computation of the tangent stiffness matrix according to the selected method
if strcmp(propCoupling.method,'penalty') || strcmp(propCoupling.method,'lagrangeMultipliers') || ...
       strcmp(propCoupling.method,'mortar') 
    computeProblemMtrcsSteadyState = @computeTangentStiffMtxResVctDDMPenaltyIGAMembrane;
elseif strcmp(propCoupling.method,'nitsche')
    computeProblemMtrcsSteadyState = @computeTangentStiffMtxResVctDDMNitscheIGAMembrane;
else
    error('Choose coupling method in propCoupling.method');
end

%% 3. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.noDOFs = 3*BSplinePatches{iPatches}.noCPs;
    if propEmpireCoSimulation.isCoSimulation
        BSplinePatches{iPatches}.noDOFsEmpire = BSplinePatches{iPatches}.noDOFs;
    end
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

%% 4. Compute the total number of DOFs for the multipatch structure
noDOFsPatches = 0;
for iPatches = 1:noPatches
    noDOFsPatches = noDOFsPatches + BSplinePatches{iPatches}.noDOFs;
end

%% 5. Create a freedom table for each patch in the multipatch geometry
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

%% 6. Assign a number of DOFs for each Lagrange Multipliers field corresponding to the multipatch coupling
if strcmp(propCoupling.method,'lagrangeMultipliers')
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

%% 7. Create a freedom table for each Lagrange Multipliers field corresponding to each patch connection
if strcmp(propCoupling.method,'lagrangeMultipliers')
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

%% 8. Compute the total number of DOFs for the coupled multipatch system
if strcmp(propCoupling.method,'lagrangeMultipliers')
    noDOFs = noDOFsPatches + noDOFsLagrangeMultipliers;
else
    noDOFs = noDOFsPatches ;
end

%% 9. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary and master-slave conditions are applied to account for the coupled system
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
freeDOFs = 1:noDOFs;
freeDOFs(ismember(freeDOFs,homDOFs)) = [];
freeDOFs(ismember(freeDOFs,inhomDOFs)) = [];

%% 10. Find the mortar (master/slave) DOFs
fixedDOFs = mergesorted(homDOFs,inhomDOFs);
if strcmp(propCoupling.method,'mortar')
    [BSplinePatches,connections] = findDOFsDDMMortarIGAMembrane...
        (BSplinePatches,connections,propCoupling,fixedDOFs);
end

%% 11. Create a DOF numbering for each patch and each Lagrange Multipliers field

% For the patches :
% _________________

for iPatches = 1:noPatches
    % Get the number of Control Points in xi-direction
    nxi = length(BSplinePatches{iPatches}.CP(:,1,1));
    
    % Get the number of Control Points in eta-direction
    neta = length(BSplinePatches{iPatches}.CP(1,:,1));
    
    % Initialize the DOF numbering array
    BSplinePatches{iPatches}.DOFNumbering = zeros(nxi,neta,3);
    
    % Initialize counter
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

% For the Lagrange Multipliers fields :
% _____________________________________

if strcmp(propCoupling.method,'lagrangeMultipliers') || strcmp(propCoupling.method,'mortar')
    for iConnections = 1:connections.No
        % Get the number of Control Points in xi-direction
        if isLagrangeMultipliersDisplacementsEnabled
            noXiLambda = length(connections.lambda{iConnections}.CP(:,1));
        end
        if isLagrangeMultipliersRotationsEnabled
            nxiMu = length(connections.mu{iConnections}.CP(:,1));
        end

        % Lagrange multipliers for the traction forces :
        % ______________________________________________

        if isLagrangeMultipliersDisplacementsEnabled
            % Initialize the field
            connections.lambda{iConnections}.DOFNumbering = zeros(noXiLambda,3);

            % Compute the entries of the DOF numbering array
            k = 1;
            for cpi = 1:noXiLambda
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

%% 12. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions and update the number of DOFs per patch if it is enriched with Lagrange Multipliers
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

%% 13. Solve the transient problem
[dHatHistory,resHistory,BSplinePatches,propCoupling,minElAreaSize] = ...
    solve_IGATransientAnalysis...
    (analysis,BSplinePatches,connections,freeDOFs,homDOFs,inhomDOFs,...
    valuesInhomDOFs,masterDOFs,slaveDOFs,nodesALE,computeInitCnds,...
    @solve_IGANLinearSystem,computeConstantProblemMatrices,...
    computeMassMtx,computeProblemMtrcsSteadyState,computeUpdatedMesh,...
    computeUpdatedGeometry,solve_LinearSystem,propCoupling,...
    propStrDynamics,propNLinearAnalysis,propPostproc,caseName,...
    pathToOutput,title,propOutput,isReferenceUpdated,...
    propEmpireCoSimulation,tab,outMsg);

%% 14. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;
    
    fprintf('Transient geometrically nonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('____________Transient Nonlinear Analysis Ended______________\n');
    fprintf('#############################################################\n\n\n');
    fprintf('\n');
end

end
