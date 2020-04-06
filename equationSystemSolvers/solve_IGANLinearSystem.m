function [u, CPHistory, resHistory, isConverged, FComplete, rankD, condK, ...
    minEig, BSplinePatches, propCoupling, minElAreaSize] = ...
    solve_IGANLinearSystem...
    (analysis, uSaved, uDotSaved, uDDotSaved, BSplinePatches, connections, ...
    u, uDot, uDDot, constMtx, massMtx, dampMtx, computeNLinearMtrcsSteadyState, ...
    computeUpdatedGeometry, freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateDirichletBCs, masterDOFs, slaveDOFs, solve_LinearSystem, t, ...
    propCoupling, propTransientAnalysis, propNLinearAnalysis, propIDBC, ...
    plot_IGANLinear, isReferenceUpdated, isCosimulationWithEmpire, tab, ...
    graph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the solution to a nonlinear system which correspond to the
% isogeometric discretization of the underlying field.
%
%                          Input :
%                       analysis : .type : The analysis type
%                         uSaved : The discrete solution field of the 
%                                  previous time step
%                      uDotSaved : The time derivative of the discrete 
%                                  solution field of the previous time step
%                     uDDotSaved : The second order time derivative of the 
%                                  discrete solution field of the previous 
%                                  time
%                 BSplinePatches : An array of structures {patch1,patch2,…}
%                                  each of the patch structures containing 
%                                  the following information:
%                                 .p,.q: Polynomial degrees
%                              .Xi,.Eta: knot vectors
%                                   .CP: Control Points coordinates and 
%                                        weights
%                              .isNURBS: Flag on whether the basis is a 
%                                        NURBS or a B-Spline
%                                 .NBC : Structure containing information 
%                                        on the application of the Neumann 
%                                        boundary conditions
%                                       .noCnd : Number of Neumann boundary 
%                                                conditions
%                             .xiLoadExtension : Cell array {.noCnd} 
%                                                containing the load 
%                                                extensions in the xi-
%                                                direction
%                            .etaLoadExtension : Cell array {.noCnd} 
%                                                containing the load 
%                                                extensions in the eta-
%                                                direction
%                               .loadAmplitude : Array (1,.noCnd) 
%                                                containing the load 
%                                                amplitudes
%                               .loadDirection : Array (1,.noCnd) 
%                                                containing the load 
%                                                directions
%                              .computeLoadVct : Cell array {.noCnd} 
%                                                containing the function 
%                                                name for the computation 
%                                                of the load vector
%                              .isConservative : Array (1,.noCnd) of flags 
%                                                indicating whether the 
%                                                load is conservative or 
%                                                not
%                             .isTimeDependent : Array (1,.noCnd) 
%                                                specifying whether each of
%                                                the loading  conditions is 
%                                                time dependent
%                                .DOFNumbering : Numbering of the DOFs 
%                                                sorted into a 3D array
%                                  .parameters : material parameters of the 
%                                                isogeometric structure
%                                         .int : On the numerical 
%                                                integration
%                                                     .type : 'default' or 
%                                                             'user'
%                                                    .xiNGP : No. of GPs 
%                                                             along xi-
%                                                             direction 
%                                                             for stiffness 
%                                                             entries
%                                                   .etaNGP : No. of GPs 
%                                                             along eta-
%                                                             direction for 
%                                                             stiffness 
%                                                             entries
%                                             .xiNGPForLoad : No. of GPs 
%                                                             along xi-
%                                                             direction for 
%                                                             load
%                                            .etaNGPForLoad : No. of GPs 
%                                                             along eta-
%                                                             direction for 
%                                                             load entries
%                                              .noGPForLoad : No. of GPs 
%                                                             along 
%                                                             boundary
%                                  connections : Define the connection 
%                                                between the patches:
%                                                       .No : Number of 
%                                                             connections
%                                                .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                                              ...      ...      ...  ...   ...  ...]
%                              u : Initial guess for the primary field
%                           uDot : Initial guess for the time derivative of 
%                                  the primary field
%                          uDDot : Initial guess for the second time
%                                  derivative of the primary field
%                       constMtx : The matrices which stay constant 
%                                  throughout the nonlinear computations
%                        massMtx : The mass matrix of the problem
%                        dampMtx : The damping matrix of the problem
% computeNLinearMtrcsSteadyState : Function handle for the computation of 
%                                  the matrices and vectors corresponding
%                                  to the steady-state problem
%         computeUpdatedGeometry : Function handle to the compute updated 
%                                  geometry for the system
%                       freeDOFs : The global numbering of the uncostrained 
%                                  DOFs
%                        homDOFs : The global numbering of the DOFs where
%                                  homogeneous Dirichlet boundary 
%                                  conditions are applied
%                      inhomDOFs : The global numbering of the DOFs where
%                                  inhomogeneous Dirichlet boundary 
%                                  conditions are applied
%                valuesInhomDOFs : The vector containing the magnitude of 
%                                  the prescribed values on the DOFs with
%                                  inhomogeneous Dirichlet boundary 
%                                  conditions
%             updateDirichletBCs : Function handle to the computation of
%                                  the updated prescribed values of the
%                                  inhomogeneous Dirichlet boundary
%                                  conditions
%                     masterDOFs : The global numbering of the DOFs which 
%                                  drive the master-slave relations
%                      slaveDOFs : The global numbering of the DOFs which 
%                                  are forced to be equal to the master 
%                                  DOFs
%             solve_LinearSystem : Function handle to the solver for the 
%                                  linear equation system
%                              t : The current time instance of the 
%                                  transient simulation
%                   propCoupling : Properties of the multipatch coupling
%                                   .alphaD : penalty factor for the 
%                                             displacement coupling
%                                   .alphaR : penalty factor for the 
%                                             rotation coupling
%                                     .intC : On the integration of the 
%                                             coupling interface
%          propTransientAnalysis : Transient analysis parameters:
%                                 .isStaticStep : Flag on whether the 
%                                                 current linear system
%                                                 corresponds to a static
%                                                 step before proceeding to
%                                                 the time loop
%                                       .method : Time integration method
%                                    .alphaBeta : Bossak parameter
%                                        .gamma : Bossak parameter
%                                       .TStart : Start time of the 
%                                                 simulation
%                                         .TEnd : End time of the simulation
%                                           .nT : Number of time steps
%                                           .dt : Time step
%                 .computeProblemMtrcsTransient : Function handle to the
%                                                 computation of the
%                                                 problem matrices
%                                                 according to the defined
%                                                 time integration method
%                                                 given the system matrix 
%                                                 of the steady-state
%                                                 system as well as the 
%                                                 mass matrix of the
%                                                 problem
%                            .computeUpdatedVct : Function handle to the 
%                                                 computation of the 
%                                                 updated values of the 
%                                                 discretized
%            propNLinearAnalysis : Nonlinear analysis parameters
%                                  .noLoadSteps : The number of the 
%                                                 nonlinear steps
%                                       .method : The nonlinear solution 
%                                                 method
%                                          .eps : The residual tolerance
%                                      .maxIter : The maximum number of
%                                                 nonlinear iterations
%                       propIDBC : Structure containing information on the 
%                                  inhomogeneous Dirichlet boundary 
%                                  conditions,
%                              .numCnd : Number of segements where those
%                                        conditions are applied
%                              .xiSpan : .noConditions x 2 array containing
%                                        the knot span extension of the
%                                        segments where those conditions
%                                        are applied in xi-direction
%                             .etaSpan : .noConditions x 2 array containing
%                                        the knot span extension of the
%                                        segments where those conditions
%                                        are applied in eta-direction
%                 .prescribedDirection : .noConditions x 1 array containing
%                                        the direction of the load
%                                        application for each condition
%                  .isUniqueOnBoundary : .noConditions x 1 array containing
%                                        flags on whether each of those
%                                        conditions are unique over their
%                                        application boundary
%                     .prescribedValue : Array of size .noConditions which
%                                        contains handles to functions
%                                        which determine the prescribed
%                                        values at each segment
%                          .isDominant : Flag on whether the inhomogeneous
%                                        Dirichlet boundary conditions are
%                                        dominant over the homogeneous ones
%                                        or notions x 2 array containing
%                                        the knot span extension of the
%                                        segments where those conditions
%                                        are applied in eta-direction
%                plot_IGANLinear : Function handle to the plotting of the 
%                                  current configuration throughout the 
%                                  nonlinear iterations
%             isReferenceUpdated : Flag on whether the reference geometry 
%                                  is updated
%       isCosimulationWithEmpire : Flag on whether co-simulation through
%                                  Empire is assumed
%                            tab : Tabulation for the output messages
%                          graph : On plotting graphs
%                         outMsg : On printing information during analysis 
%                                  in the command window
%
%                         Output :
%                              u : The converged discrete solution vector
%                      CPHistory : The Control Point coordinates history
%                     resHistory : Array containing the evolution of the 
%                                  residual throughout the nonlinear 
%                                  iterations and the load steps
%                    isConverged : Flag on whether the nonlinear iterations 
%                                  have converged
%                      FComplete : Dummy output for this function
%                          rankD : Dummy output for this function
%                          condK : Dummy output for this function
%                         minEig : Dummy output for this function
%                 BSplinePatches : The updated array of B-Spline patches
%                   propCoupling : Updated structure of the coupling
%                                  properties
%                  minElAreaSize : The minimum element area over the
%                                  isogeometric mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Update the time-dependent inhomogeneous Dirichlet boundary conditions
%
% 2. Initialize the load vector on the multipatch structure
% ->
%    2i. Initialize the load vector of the patch
%
%   2ii. Get the Neumann boundary conditions of the current patch
%
%  2iii. Check if there is a non-conservative loading associated with the current patch
%
%   2iv. Loop over the Neumann boundary conditions of the current patch
%   ->
%        2iv.1. Initialize the load vector for the current condition
%
%        2iv.2. Get the function handle for the load vector computation  
%
%        2iv.3. Compute the load vector and the tangent matrix resulting from the application of follower loads
%
%        2iv.4. If the loading is not conservative add the contribution to the non-conservative load vector
%
%        2iv.5. Add The compute external load vector into the B-Spline array
%   <-
% <-
% 
% 3. Loop over all the fixed point (coupling) iterations if external coupling through Empire is considered
% ->
%    3i. Initialize the force vector to be received by Empire
%
%   3ii. Receive field from Empire
%
%  3iii. Distribute the received forces from Empire to each patch
%
%   3iv. Loop over all load steps
%   ->
%        3iv.1. Compute the load factor
%
%        3iv.2. Loop over all the nonlinear iterations
%        ->
%               3iv.2i. Compute the tangent stiffness matrix and the residual vector for the steady-state problem
%
%              3iv.2ii. Compute the tangent stiffness matrix and the residual vector for the transient problem
%
%             3iv.2iii. Re-arrangement of the equation system for specific treatment of condition enforcement
%
%              3iv.2iv. Get the array of inhomogeneous Dirichlet boundary conditions
%
%               3iv.2v. Compute the right-hand side (RHS) residual vector in equation
%
%              3iv.2vi. Check condition for convergence on the residual vector
%
%             3iv.2vii. Solve the linearized equation system
%
%            3iv.2viii. Re-assemble to the complete vector of unknowns
%
%              3iv.2ix. Update the patch geometry
%
%               3iv.2x. Loop over the patches and update the external load vector if a nonconservative load is encountered
%               ->
%                       3iv.2x.1. Get the Neumann boundary conditions of the current patch
%
%                       3iv.2x.2. Erase the nonconservative part from the load application and the corresponding tangent matrix in order to update it
%
%                       3iv.2x.3. Loop over the Neumann boundary conditions of the current patch
%                       ->
%                                 3iv.2x.3i. Check if the load application is conservative or not and if the current linear system corresponds to the static step and the load is transient and continue if it yes
%
%                                3iv.2x.3ii. Get the function handle for the computation of the load vector at the current patch
%
%                               3iv.2x.3iii. Compute the non-conservative load vector corresponding to the next iteration step and the tangent matrix resulting from the application of follower loads
%                       <-
%
%                       3iv.2x.4. Add the nonconservative part to the load vector
%               <-
%
%              3iv.2xi. Save the results of the Control Point deformation
%        <-
%
%        3iv.3. If the nonlinear iterations have not converged break the loop
%
%        3iv.4. Plot the current configuration
%   <-
%
%    3v. Send the converged solution to Empire
%
%   3vi. Check convergence of the fixed point (coupling) iterations
%
%  3vii. Subtract the received from Empire forces from each patch 
% <-
%
% 4. Close figures
%
%% Function main body

%% 0. Read input

% Define a tolerance
eps = 1e-14;

% Print out message on the nonlinear method under use
if strcmp(outMsg, 'outputEnabled')
    fprintf(strcat(tab, 'Nonlinear method : %s\n\n'), propNLinearAnalysis.method);
end

% Get number of patches
noPatches = length(BSplinePatches);
if noPatches == 1
    noPatch = 1;
else
    noPatch = 'undefined';
end

% Get the number of weak Dirichlet boundary conditions
noWeakDBCCnd = 0;
for iPatches = 1:noPatches
    if isfield(BSplinePatches{iPatches}, 'weakDBC')
        if isfield(BSplinePatches{iPatches}.weakDBC, 'noCnd')
            noWeakDBCCnd = noWeakDBCCnd + ...
                BSplinePatches{iPatches}.weakDBC.noCnd;
        end
    end
end

% Get the time step number if the analysis is transient
if ~ischar(propTransientAnalysis)
    if isfield(propTransientAnalysis, 'timeDependence')
        if strcmp(propTransientAnalysis.timeDependence, 'transient') || ...
                strcmp(propTransientAnalysis.timeDependence, 'pseudotransient')
            noTimeStep = (t - propTransientAnalysis.TStart)/propTransientAnalysis.dt + 1;
            if mod(noTimeStep,1) ~= 0
                if abs(noTimeStep - round(noTimeStep)) > eps*noTimeStep/propTransientAnalysis.dt
                    error('Time step number was computed as non-integer')
                end
                noTimeStep = round(noTimeStep);
            end
        end
    end
end
if ~exist('noTimeStep', 'var')
    noTimeStep = 'undefined';
end

% Get the coupling type in case of co-simulation using Empire
if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
    couplingType = EMPIRE_API_getUserDefinedText('couplingType');
    if strcmp(couplingType, 'looseCoupling')
        isIterativeCoupling = false;
    elseif strcmp(couplingType, 'iterativeCoupling')
        isIterativeCoupling = true;
    else
        error('Coupling type in propEmpireCoSimulation.strMatlabXml is %s but it has to be defined as either looseCoupling or iterativeCoupli', ...
            couplingType);
    end
end

% Flag on the fixed point iterations
isFixedPointConvergent = false;

% Total number of DOFs when coupling through Empire is considered
if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
    noDOFsEmpire = 0;
    for iPatches = 1:noPatches
        noDOFsEmpire = noDOFsEmpire + BSplinePatches{iPatches}.noDOFsEmpire;
    end
end

% Get the global numbering of the DOFs which are sent to Empire
if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
    noDOFsPrevious = 0;
    EFTEmpire = zeros(noDOFsEmpire, 1);
    for iPatches = 1:noPatches
        if iPatches ~= 1
            noDOFsPrevious = noDOFsPrevious + BSplinePatches{iPatches - 1}.noDOFsEmpire;
        end
        EFTEmpire(noDOFsPrevious + 1:noDOFsPrevious + BSplinePatches{iPatches}.noDOFsEmpire, 1) = ...
            BSplinePatches{iPatches}.EFTPatches(1,1:BSplinePatches{iPatches}.noDOFsEmpire);
    end
end

% Dummy output
FComplete = 'undefined';
rankD = 'undefined';
condK = 'undefined';
minEig = 'undefined';

% Initialize the displacement history array for each patch
CPHistory = struct([]);
for iPatches = 1:noPatches
    CPHistory{iPatches} = zeros(length(BSplinePatches{iPatches}.CP(:, 1, 1)), ...
        length(BSplinePatches{iPatches}.CP(1, :, 1)), ...
        length(BSplinePatches{iPatches}.CP(1, 1, :)), ...
        propNLinearAnalysis.noLoadSteps + 1);
    CPHistory{iPatches}(:, :, :, 1) = BSplinePatches{iPatches}.CP;
end

% Check if there exist a nonconservative loading
tanMtxLoad = struct([]);
isConservative = true(noPatches, 1);
for iPatches = 1:noPatches
    for iNBC = 1:BSplinePatches{iPatches}.NBC.noCnd
        if BSplinePatches{iPatches}.NBC.isFollower(iNBC, 1)
            isConservative(iPatches, 1) = false;
            break;
        end
    end
end

% Initialize the tangent matrix due to follower loading for each patch
for iPatches = 1:noPatches
    if ~isConservative(iPatches, 1)
        tanMtxLoad{iPatches} = zeros(BSplinePatches{iPatches}.noDOFs);
    else
        tanMtxLoad{iPatches} = 'undefined';
    end
end

% Initialize the residual history array
resHistory = zeros(propNLinearAnalysis.maxIter, propNLinearAnalysis.noLoadSteps);

%% 1. Update the time-dependent inhomogeneous Dirichlet boundary conditions
if ~ischar(propIDBC)
    if isa(updateDirichletBCs, 'function_handle')
        [homDOFs, inhomDOFs, valuesInhomDOFs] = ...
            updateDirichletBCs(BSplinePatches, homDOFs, propIDBC, t);
    else
        error('Variable updateInhomDOFs should define a function handle');
    end
end

%% 2. Initialize the load vector on the multipatch structure
for iPatches = 1:noPatches
    %% 2i. Initialize the load vector of the patch
    BSplinePatches{iPatches}.FGamma = zeros(BSplinePatches{iPatches}.noDOFs, 1);
    
    %% 2ii. Get the Neumann boundary conditions of the current patch
    NBC = BSplinePatches{iPatches}.NBC;
    
    %% 2iii. Check if there is a non-conservative loading associated with the current patch
    if ~isConservative(iPatches, 1)
        BSplinePatches{iPatches}.FNonConservative = ...
            zeros(BSplinePatches{iPatches}.noDOFs, 1);
    end

    %% 2iv. Loop over the Neumann boundary conditions of the current patch
    for iNBC = 1:NBC.noCnd
        %% 2iv.1. Initialize the load vector for the current condition
        FGamma = zeros(BSplinePatches{iPatches}.noDOFs, 1);
        
        %% 2iv.2. Get the function handle for the load vector computation
        funcHandle = str2func(NBC.computeLoadVct{iNBC});
        
        %% 2iv.3. Compute the load vector and the tangent matrix resulting from the application of follower loads
        if ~(propTransientAnalysis.isStaticStep && NBC.isTimeDependent(iNBC, 1))
            [FGamma, tanMtxLoadPatch] = funcHandle ...
                (FGamma,BSplinePatches{iPatches}, ...
                NBC.xiLoadExtension{iNBC}, NBC.etaLoadExtension{iNBC}, ...
                NBC.loadAmplitude{iNBC}, NBC.loadDirection{iNBC}, ...
                NBC.isFollower(iNBC, 1), t, BSplinePatches{iPatches}.int, ...
                '');
            if NBC.isFollower(iNBC, 1)
                tanMtxLoad{iPatches} = tanMtxLoad{iPatches} + tanMtxLoadPatch;
            end
        end
        
        %% 2iv.4. If the loading is not conservative add the contribution to the non-conservative load vector
        if NBC.isFollower(iNBC, 1)
            BSplinePatches{iPatches}.FNonConservative = BSplinePatches{iPatches}.FNonConservative + ...
                FGamma;
        end
        
        %% 2iv.5. Add The compute external load vector into the B-Spline array
        BSplinePatches{iPatches}.FGamma = BSplinePatches{iPatches}.FGamma + ...
            FGamma;
    end
end

%% 3. Loop over all the fixed point (coupling) iterations if external coupling through Empire is considered
while ~isFixedPointConvergent
    %% 3i. Initialize the force vector to be received by Empire
    if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
        forceVctEmpire = zeros(1, noDOFsEmpire);
    end
    
    %% 3ii. Receive field from Empire
    if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
        fprintf(strcat(tab, 'Receiving field from Empire\n\n'));
        EMPIRE_API_recvDataField('defaultField', noDOFsEmpire, forceVctEmpire);
    end
    
    %% 3iii. Distribute the received forces from Empire to each patch
    if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
        for iPatches = 1:noPatches
            BSplinePatches{iPatches}.FGamma(1:BSplinePatches{iPatches}.noDOFsEmpire, 1) = ...
                BSplinePatches{iPatches}.FGamma(1:BSplinePatches{iPatches}.noDOFsEmpire, 1) + ...
                forceVctEmpire(1, BSplinePatches{iPatches}.EFTPatches(1, 1:BSplinePatches{iPatches}.noDOFsEmpire))';
        end
    end
    
    %% 3iv. Loop over all load steps
    for iLoadStep = 1:propNLinearAnalysis.noLoadSteps
        if strcmp(outMsg, 'outputEnabled')
            fprintf(strcat(tab, 'Load step %d/%d \n'), iLoadStep, ...
                propNLinearAnalysis.noLoadSteps);
            fprintf(strcat(tab, '----------------\n'));
            fprintf('\n');
        end
        %% 3iv.1. Compute the load factor
        loadFactor = iLoadStep/propNLinearAnalysis.noLoadSteps;

        %% 3iv.2. Loop over all the nonlinear iterations
        if strcmp(outMsg,'outputEnabled')
            msgPNR = sprintf(strcat(tab, '\tLooping over the nonlinear iterations\n', ...
                tab, '\t-------------------------------------\n\n'));
            fprintf(msgPNR);
        end
        for iNLinearIter = 1:propNLinearAnalysis.maxIter
            %% 3iv.2i. Compute the tangent stiffness matrix and the residual vector for the steady-state problem
            [tanMtx, resVct, BSplinePatches, propCoupling, minElAreaSize] = ...
                computeNLinearMtrcsSteadyState ...
                (constMtx, tanMtxLoad, u, uSaved, uDot, uDotSaved, ...
                BSplinePatches, connections, propCoupling, loadFactor, ...
                noPatch, noTimeStep, iNLinearIter, noWeakDBCCnd, t, ...
                propTransientAnalysis, isReferenceUpdated, ...
                strcat(tab, '\t'), outMsg);

            %% 3iv.2ii. Compute the tangent stiffness matrix and the residual vector for the transient problem
            if isa(propTransientAnalysis.computeProblemMtrcsTransient,'function_handle') && ...
                ~propTransientAnalysis.isStaticStep
                [tanMtx,resVct] = propTransientAnalysis.computeProblemMtrcsTransient ...
                    (u, uSaved, uDot, uDotSaved, uDDot, uDDotSaved, ...
                    massMtx, dampMtx, tanMtx, resVct, propTransientAnalysis);
            elseif ~isa(propTransientAnalysis.computeProblemMtrcsTransient, 'function_handle') && ...
                    ~(strcmp(propTransientAnalysis.timeDependence, 'steadyState') || ...
                    strcmp(propTransientAnalysis.timeDependence, 'pseudotransient'))
                error('Variable propTransientAnalysis.computeProblemMtrcsTransient is undefined but the simulation is transient');
            end

            %% 3iv.2iii. Re-arrangement of the equation system for specific treatment of condition enforcement
            if ~ischar(propCoupling)
                if isfield(propCoupling, 'method')
                    if strcmp(propCoupling.method, 'mortar')
                        if ~isfield(propCoupling, 'computeRearrangedProblemMtrcs')
                            error('For the mortar method a function handle under propCoupling.computeRearrangedProblemMtrcs needs to be specified');
                        end
                        [tanMtx, resVct] = propCoupling.computeRearrangedProblemMtrcs ...
                            (constMtx, BSplinePatches, connections, tanMtx, resVct);
                    end
                end
            end

            %% 3iv.2iv. Get the array of inhomogeneous Dirichlet boundary conditions
            if iscell(valuesInhomDOFs)
                if length(valuesInhomDOFs) ~= length(inhomDOFs)
                    error('The arrays of the inhomogeneous Dirichlet boundary conditions and their values are not the same');
                end
                valuesInhomDOFsTemp = zeros(1, length(inhomDOFs));
                for iInhomDBC = 1:length(inhomDOFs)
                    valuesInhomDOFsTemp(1, iInhomDBC) = valuesInhomDOFs{iInhomDBC}(t);
                end
                clear valuesInhomDOFs;
                valuesInhomDOFs = valuesInhomDOFsTemp;
            end

            %% 3iv.2v. Compute the right-hand side (RHS) residual vector in equation
            RHS = - resVct;
            if norm(valuesInhomDOFs) ~= 0 && iNLinearIter == 1
                RHS = RHS + tanMtx(:, inhomDOFs)*valuesInhomDOFs';
            end

            %% 3iv.2vi. Check condition for convergence on the residual vector

            % Compute the norm of the residual vector over the free DOFs
            resHistory(iNLinearIter, iLoadStep) = norm(resVct(freeDOFs));

            % Issue a message on the evolution of the residual vector
            if strcmp(outMsg, 'outputEnabled')
                msgNR = sprintf(strcat(tab, '\t||resVct|| = %d at nonlinear iteration No. %d \n'), ...
                    resHistory(iNLinearIter, iLoadStep), iNLinearIter);
                fprintf(msgNR);
            end

            % Check the convergence of the nonlinear iterations
            if resHistory(iNLinearIter, iLoadStep) <= propNLinearAnalysis.eps && iNLinearIter ~= 1
                if strcmp(outMsg, 'outputEnabled')
                    msgANR = sprintf(strcat('\n', tab, '\tNonlinear iterations converged! \n \n \n'));
                    fprintf(msgANR);
                end
                isConverged = true;
                break;
            end

            % If the nonlinear iterations do not converge after specified limit:
            if iNLinearIter == propNLinearAnalysis.maxIter
                if strcmp(outMsg, 'outputEnabled')
                    warning('Nonlinear iterations did not converge');
                end

                % Flag on the convergence of the nonlinear iterations
                isConverged = false;
            end

            %% 3iv.2vii. Solve the linearized equation system
            [deltauRed, isLinearSystemConverged] = solve_LinearSystem ...
                (tanMtx(freeDOFs, freeDOFs), RHS(freeDOFs), u(freeDOFs));
            if ~isLinearSystemConverged
                error('Linear equation solver did not converge');
            end

            %% 3iv.2viii. Re-assemble to the complete vector of unknowns
            u(freeDOFs) = u(freeDOFs) + deltauRed;
            if iNLinearIter == 1
                u(homDOFs) = 0;
                u(inhomDOFs) = valuesInhomDOFs;
            end

            %% 3iv.2ix. Update the patch geometry
            if isa(computeUpdatedGeometry, 'function_handle')
                BSplinePatches = computeUpdatedGeometry(BSplinePatches, u);
            end

            %% 3iv.2x. Loop over the patches and update the external load vector if a nonconservative load is encountered
            for iPatches = 1:noPatches
                %% 3iv.2x.1. Get the Neumann boundary conditions of the current patch
                NBC = BSplinePatches{iPatches}.NBC;

                %% 3iv.2x.2. Erase the nonconservative part from the load application and the corresponding tangent matrix in order to update it
                if ~isConservative(iPatches, 1)
                    BSplinePatches{iPatches}.FGamma = BSplinePatches{iPatches}.FGamma - ...
                        BSplinePatches{iPatches}.FNonConservative;
                    BSplinePatches{iPatches}.FNonConservative = ...
                        zeros(BSplinePatches{iPatches}.noDOFs, 1);
                    tanMtxLoad{iPatches} = zeros(BSplinePatches{iPatches}.noDOFs);
                end

                %% 3iv.2x.3. Loop over the Neumann boundary conditions of the current patch
                for iNBC = 1:NBC.noCnd
                    %% 3iv.2x.3i. Check if the load application is conservative or not and if the current linear system corresponds to the static step and the load is transient and continue if it yes
                    if ~NBC.isFollower(iNBC, 1) || ...
                            (propTransientAnalysis.isStaticStep && NBC.isTimeDependent(iNBC, 1))
                        continue; 
                    end

                    %% 3iv.2x.3ii. Get the function handle for the computation of the load vector at the current patch
                    funcHandle = str2func(NBC.computeLoadVct{iNBC});

                    %% 3iv.2x.3iii. Compute the non-conservative load vector corresponding to the next iteration step and the tangent matrix resulting from the application of follower loads
                    [BSplinePatches{iPatches}.FNonConservative, tanMtxLoadPatch] = ...
                        funcHandle ...
                        (BSplinePatches{iPatches}.FNonConservative, ...
                        BSplinePatches{iPatches}, NBC.xiLoadExtension{iNBC}, ...
                        NBC.etaLoadExtension{iNBC}, NBC.loadAmplitude{iNBC}, ...
                        NBC.loadDirection{iNBC}, NBC.isFollower(iNBC, 1), ...
                        t, BSplinePatches{iPatches}.int, '');
                    if NBC.isFollower(iNBC, 1)
                        tanMtxLoad{iPatches} = tanMtxLoad{iPatches} + tanMtxLoadPatch;
                    end
                end

                %% 3iv.2x.4. Add the nonconservative part to the load vector
                if ~isConservative(iPatches, 1)
                    BSplinePatches{iPatches}.FGamma = BSplinePatches{iPatches}.FGamma + ...
                        BSplinePatches{iPatches}.FNonConservative;
                end
            end

            %% 3iv.2xi. Save the results of the Control Point deformation
            for iPatches = 1:noPatches
                CPHistory{iPatches}(:, :, :, iLoadStep + 1) = ...
                    BSplinePatches{iPatches}.CPd;
            end
        end

        %% 3iv.3. If the nonlinear iterations have not converged break the loop
        if ~isConverged
            break;
        end

        %% 3iv.4. Plot the current configuration
        if ~ischar(plot_IGANLinear)
            figure(graph.index)
            graph.index = plot_IGANLinear(BSplinePatches, u, graph, '');
            graph.index = graph.index - 1;
        end
    end
    
    %% 3v. Send the converged solution to Empire
    if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
        fprintf(strcat(tab, 'Sending converged solution to Empire\n\n'));
        EMPIRE_API_sendDataField('defaultField', noDOFsEmpire, u(EFTEmpire));
    end
        
    %% 3vi. Check convergence of the fixed point (coupling) iterations
    if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
        if isIterativeCoupling
            convergenceSignal = EMPIRE_API_recvConvergenceSignal();
            if convergenceSignal == 1
                isFixedPointConvergent = true;
            end
            fprintf(strcat(tab, 'Fixed point (coupling) iterations converged\n\n'));
        else
            isFixedPointConvergent = true;
        end
    else
        isFixedPointConvergent = true;
    end
    
    %% 3vii. Subtract the received from Empire forces from each patch
    if isCosimulationWithEmpire && ~propTransientAnalysis.isStaticStep
        for iPatches = 1:noPatches
            BSplinePatches{iPatches}.FGamma(1:BSplinePatches{iPatches}.noDOFsEmpire, 1) = ...
                BSplinePatches{iPatches}.FGamma(1:BSplinePatches{iPatches}.noDOFsEmpire, 1) - ...
                forceVctEmpire(1, BSplinePatches{iPatches}.EFTPatches(1, 1:BSplinePatches{iPatches}.noDOFsEmpire))';
        end
    end
end

%% 4. Close figures
if ~ischar(plot_IGANLinear)
    close(figure(graph.index));
end

end
