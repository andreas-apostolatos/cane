function [eigenmodeShapes, naturalFrequencies, dHat] = ...
    solve_eigenmodeAnalysisIGAMembrane ...
    (BSplinePatch, solve_LinearSystem, noEig, propNLinearAnalysis, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the natural frequencies and the corresponding eigenmode shapes
% for the linear isogeometric membrane problem.
%
%               Input :
%        BSplinePatch : The parameters defining a single patch model
%                       .p,.q : The polynomial orders of the patch in xi-
%                               and eta- parametric directions
%                    .Xi,.Eta : The knot vectors of the patch in xi- and 
%                               eta- parametric directions
%                         .CP : The set of Control Poin coordinates and
%                               weights
%                    .isNURBS : Flag on whether the underlying basis is a
%                               NURBS or a B-Spline
%                 .parameters : Technical parameters defining the
%                               structural problem
%                    .homDOFs : The global numbering of the DOFs where
%                               homogeneous Dirichlet boundary conditions
%                               are applied
%                  .inhomDOFs : The global numbering of the DOFs where
%                               inhomogeneous Dirichlet boundary conditions
%                               are applied
%            .valuesInhomDOFs : The corresponding values to the
%                               inhomogeneous Dirichlet boundary conditions
%                        .NBC : Structure defining the Neumann boundary
%                               conditions :
%                                         .noCnd : Number of conditions
%                                   .xiExtension : Load extensions in
%                                                  xi-parametric directions
%                                 .etaExtensions : Load extensions in
%                                                  eta-parametric 
%                                                  directions
%                                 .loadAmplitude : Definition of the load
%                                                  amplitude
%                                 .loadDirection : Definition of the
%                                                  direction of the load
%                                .computeLoadVct : Function handles to the
%                                                  computation of the load 
%                                                  vector
%                                .isConservative : Flags defining whether
%                                                  the corresponding
%                                                  conditions are
%                                                  conservative or not
%                               .isTimeDependent : Flag on whether the
%                                                  corresponding condition
%                                                  is time dependent or not
%                     .cables : Structure defining the embedded cables :
%                                       .No : Number of cables
%                              .xiExtension : Cable extensions in xi-
%                                             parametric directions
%                             .etaExtension : Cable extensions in eta-
%                                             parametric directions
%                               .parameters : Technical parameters of the
%                                             cables
%                                      .int : Quadrature type
%          solve_LinearSystem : Function handle to the solution of linear
%                               equation systems
%                       noEig : Number of eigenmode shapes to be returned
%                               it can either be a positive integer or a
%                               character 'all' if all eigenvalues are
%                               requested
%         propNLinearAnalysis : Properties of the nonlinear method for
%                               tackling the nonlinear problem :
%                                    .method : 'newtonRapshon'
%                               .noLoadSteps : Number of load steps
%                                       .eps : Residual tolerance
%                                   .maxIter : Maximum number of nonlinear
%                                              iterations
%                      outMsg : Allows for printing information onto the
%                               command window when chosen 'outputEnabled'
%
%                      Output :
%             eigenmodeShapes : The eigenmode shapes into a matrix
%          naturalFrequencies : The natural frequencies of the discrete
%                               system in descending order
%                        dHat : The displacement field when solving the
%                               problem for bringing it in equilibrium with
%                               its internal stresses
%
% Function layout :
%
% 0. Read input
%
% 1. Find the total number of DOFs for the patch including the Lagrange Multipliers and make an EFT for each Lagrange Multipliers field employed
%
% 2. Compute the constant matrices of the patch corresponding to the application of weak boundary conditions
%
% 3. Compute an empty load vector for each patch
%
% 4. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
%
% 5. Initialize the solution field
%
% 6. Perform a static solve to bring the structure in equilibrium with its internal forces
%
% 7. Compute the static external loads
%    
%    7i. Initialize the load vector and the tangent stiffness matrix
%
%   7ii. Get the Neumann boundary conditions of the current patch
%
%  7iii. Check if there is a non-conservative loading associated with the current patch
%
%   7iv. Loop over the Neumann boundary conditions of the current patch
%   ->
%        7iv.1. Get the function handle for the load vector computation
%
%        7iv.2. Compute the load vector and the tangent matrix resulting from the application of follower loads
%
%        7iv.3. If the loading is not conservative add the contribution to the non-conservative load vector
%
%        7iv.4. Add The compute external load vector into the B-Spline array
%   <-
%
% 8. Compute the linear stiffness matrix of the membrane problem
%
% 9. Compute the mass matrix of the membrane problem
%
% 10. Solve the generalized eigenvalue problem to get the eigenmode shapes
%
% 11. Apply the boundary conditions onto the system
%
% 12. Appendix
%
%% Function main body
if strcmp(outMsg, 'outputEnabled')
    fprintf('_______________________________________________________________\n');
    fprintf('###############################################################\n');
    fprintf('Modal analysis for the isogeometric membrane has been initiated\n');
    fprintf('\n');
    if isnumeric(noEig)
        if mod(noEig, 1) ~= 0
            error('The number of eigenfrequencies must be an integer');
        end
        if noEig <= 0
            error('The number of eigenfrequencies must be strictly positive');
        end
        fprintf('Number of eigenfrequencies requested : %d', noEig);
    else
        error('Variable noEig must be numeric');
    end
    fprintf('\n');
    isWeakDBC = false;
    if isfield(BSplinePatch, 'weakDBC')
        if isfield(BSplinePatch.weakDBC, 'noCnd')
            if BSplinePatch.weakDBC.noCnd > 0
                isWeakDBC = true;
            end
        end
    end
	if isWeakDBC
        fprintf('Weak Dirichlet boundary conditions \n');
        fprintf('---------------------------------- \n\n');
        if isfield(BSplinePatch, 'weakDBC')
            if isfield(BSplinePatch.weakDBC, 'noCnd')
                if BSplinePatch.weakDBC.noCnd > 0
                    fprintf('Weak boundary conditions using the %s method are applied \n', ...
                        BSplinePatch.weakDBC.method);
                    if strcmp(BSplinePatch.weakDBC.method, 'Nitsche')
                        if BSplinePatch.weakDBC.estimationStabilPrm
                            fprintf('Automatic estimation of the stabilization parameter is enabled \n');
                        end
                        if isfield(BSplinePatch.weakDBC, 'computeConstMtx')
                            if isfield(BSplinePatch.weakDBC.alpha)
                                fprintf('Manual stabilization parameter chosen as %d\n', ...
                                    BSplinePatch.weakDBC.alpha)
                            else
                                error('Manual stabilization parameter weakDBC.alpha needs to be assigned\n');
                            end
                        end
                    elseif strcmp(BSplinePatch.weakDBC.method, 'Penalty')
                        if isfield(BSplinePatch.weakDBC.alpha)
                            fprintf('The penalty parameter is chosen as %d', ...
                                BSplinePatch.weakDBC.alpha);
                        else
                            error('The penalty parameter needs to be assigned');
                        end
                    end
                    fprintf('\n');
                end
            end
        end
	end
    fprintf('\n');
    fprintf('_______________________________________________________________\n');
    fprintf('\n');
    tic;
end

%% 0. Read input

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Create the multipatch structure for the B-Splines
BSplinePatches = {BSplinePatch};

% Tabulation of the output messages
tab = '\t';

% Counter of the nonlinear iterations
counterNonlinearIterations = 1;

% Initialize the dummy arrays
dHatDot = 'undefined';
dHatDDot = 'undefined';
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
noTimeStep = 'undefined';
noPatch = 'undefined';
propCoupling = 'undefined';
constStiff = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
connections = 'undefined';
masterDOFs = 'undefined';
slaveDOFs = 'undefined';
updateDirichletBCs = 'undefined';
propIDBC = 'undefined';
plot_IGANLinear = 'undefined';
graph = 'undefined';

% Co-simulation with Empire
isCosimulationWithEmpire = false;

% Initialize tangent matrix contribution from follower load
tanMtxLoad = struct([]);

% Initialize the time
t = 0;

% Assign flag on whether the reference geometry is updated
isReferenceUpdated = false;

% Define function handle for the computation of the updated geometry
computeUpdatedGeometry = @computeUpdatedGeometryIGAThinStructureMultipatches;

% Get the number of weak Dirichlet boundary conditions
noWeakDBCCnd = 0;
if isfield(BSplinePatches{1}, 'weakDBC')
    if isfield(BSplinePatches{1}.weakDBC, 'noCnd')
        noWeakDBCCnd = BSplinePatches{1}.weakDBC.noCnd;
    end
end

% Check if there exist a nonconservative loading
isConservative = true;
for iNBC = 1:BSplinePatches{1}.NBC.noCnd
    if BSplinePatches{1}.NBC.isFollower(iNBC,1)
        isConservative = false;
        break;
    end
end

% Static linear analysis
propTransientAnalysis.timeDependence = 'steadyState';
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;

% Load factor definition
loadFactor = 1;

% Function computation to the stiffness matrix of the membrane problem
if BSplinePatches{1}.noElmnts == 1
    computeTangentStiffMtxResVct = @computeTangentStiffMtxResVctIGAMembraneNLinearOutdated;
else
    computeTangentStiffMtxResVct = @computeTangentStiffMtxResVctIGAMembraneNLinear;
end

%% 1. Find the total number of DOFs for the patch including the Lagrange Multipliers and make an EFT for each Lagrange Multipliers field employed

% Assign an EFT for each Lagrange Multipliers field
noDOFsPatchLM = 0;
if isfield(BSplinePatches{1}.weakDBC, 'computeConstMtx')
    if strcmp(BSplinePatches{1}.weakDBC.computeConstMtx, ...
            'computeWeakDBCMtxLagrangeMultipliersIGAMembrane')
        for iCnd = 1:BSplinePatches{1}.weakDBC.noCnd
            noDOFsPatchLMCnd = 3*length(BSplinePatches{1}.weakDBC.lambda{iCnd}.CP(:, 1));
            noDOFsPatchLM = noDOFsPatchLM + noDOFsPatchLMCnd;
            if iCnd == 1
                index = 3*BSplinePatches{1}.noCPs;
            else
                index = BSplinePatches{1}.weakDBC.lambda{iCnd - 1}.EFT(length(BSplinePatches{1}.weakDBC.lambda{iCnd - 1}.EFT));
            end
            BSplinePatches{1}.weakDBC.lambda{iCnd}.EFT = index + 1:index + noDOFsPatchLMCnd;
        end
    end
end
BSplinePatches{1}.noDOFs = 3*BSplinePatches{1}.noCPs + noDOFsPatchLM;
if noEig > BSplinePatches{1}.noDOFs
    error('The number of  requested eigenfrequencies can be up to %', BSplinePatches{1}.noDOFs);
end

% Compute the number of DOFs
numDOFs = BSplinePatches{1}.noDOFs;

% Create the element freedom table for the BSplinePatch into the array of
% the patches
BSplinePatches{1}.EFTPatches = 1:BSplinePatches{1}.noDOFs;

% Create a DOF numbering for each Lagrange Multipliers field of the patch
if isfield(BSplinePatches{1}.weakDBC, 'computeConstMtx')
    if strcmp(BSplinePatches{1}.weakDBC.computeConstMtx, ...
            'computeWeakDBCMtxLagrangeMultipliersIGAMembrane')
        for iCnd = 1:BSplinePatches{1}.weakDBC.noCnd
            % Get the number of Control Points
            nxiLambda = length(BSplinePatches{1}.weakDBC.lambda{iCnd}.Xi);

            % Initialize the field of the DOF numbering
            BSplinePatches{1}.weakDBC.lambda{iCnd}.DOFNumbering = zeros(nxiLambda, 3);

             % Compute the entries of the DOF numbering array
            k = 1;
            for cpi = 1:nxiLambda
                BSplinePatches{1}.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 1) = k;
                BSplinePatches{1}.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 2) = k + 1;
                BSplinePatches{1}.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 3) = k + 2;

                % Update counter
                k = k + 3;
            end
        end
    end
end

%% 2. Compute the constant matrices of the patch corresponding to the application of weak boundary conditions
if strcmp(outMsg, 'outputEnabled')
    message = 'Computing the constant matrices related to the application of weak Dirichlet boundary conditions\n';
    fprintf([tab, '>>', ' ', message]);
end
if isfield(BSplinePatches{1}.weakDBC, 'computeConstMtx') && BSplinePatches{1}.weakDBC.noCnd > 0
    computeWeakDBCConstantProblemMatrices = str2func(BSplinePatches{1}.weakDBC.computeConstMtx);
    BSplinePatches{1}.KConstant = ...
        computeWeakDBCConstantProblemMatrices...
        (BSplinePatches{1}, connections, numDOFs, propCoupling);
else
    BSplinePatches{1}.KConstant = 'undefined';
end

%% 3. Compute an empty load vector for each patch
BSplinePatches{1}.FGamma = zeros(BSplinePatches{1}.noDOFs, 1);

%% 4. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
homDOFs = BSplinePatches{1}.homDOFs;
inhomDOFs = BSplinePatches{1}.inhomDOFs;
valuesInhomDOFs = BSplinePatches{1}.valuesInhomDOFs;
freeDOFs = 1:numDOFs;
freeDOFs(ismember(freeDOFs, homDOFs)) = [];
freeDOFs(ismember(freeDOFs, inhomDOFs)) = [];

%% 5. Initialize the solution field
dHat = zeros(numDOFs, 1);

%% 6. Perform a static solve to bring the structure in equilibrium with its internal forces
if strcmp(outMsg,'outputEnabled')
    message = 'Performing steady-state analysis to bring the structure in equilibrium\n';
    fprintf([tab, '>>', ' ', message]);
end
[dHat, ~, ~, ~, ~, ~, ~, ~, BSplinePatches, propCoupling, ~] = ...
    solve_IGANLinearSystem...
    (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
    connections, dHat, dHatDot, dHatDDot, constStiff, massMtx, dampMtx, ...
    computeTangentStiffMtxResVct, computeUpdatedGeometry, freeDOFs, ...
    homDOFs, inhomDOFs, valuesInhomDOFs, updateDirichletBCs, masterDOFs, ...
    slaveDOFs, solve_LinearSystem, t, propCoupling, propTransientAnalysis, ...
    propNLinearAnalysis, propIDBC, plot_IGANLinear, isReferenceUpdated, ...
    isCosimulationWithEmpire, strcat(tab, '\t'), graph, outMsg);

%% 7. Compute the static external loads
%% 7i. Initialize the load vector and the tangent stiffness matrix
FGamma = zeros(BSplinePatches{1}.noDOFs, 1);
tanMtxLoad{1} = zeros(BSplinePatches{1}.noDOFs);
BSplinePatches{1}.FGamma = FGamma;

%% 7ii. Get the Neumann boundary conditions of the current patch
NBC = BSplinePatches{1}.NBC;

%% 7iii. Check if there is a non-conservative loading associated with the current patch
if ~isConservative
    BSplinePatches{1}.FNonConservative = ...
        zeros(BSplinePatches{1}.noDOFs, 1);
end

%% 7iv. Loop over the Neumann boundary conditions of the current patch
% warning([tab 'Tangent stiffness matrix computation does not take into account the follower loads!']);
for iNBC = 1:NBC.noCnd
    %% 7iv.1. Get the function handle for the load vector computation
    funcHandle = str2func(NBC.computeLoadVct{iNBC});

    %% 7iv.2. Compute the load vector and the tangent matrix resulting from the application of follower loads
    if ~(propTransientAnalysis.isStaticStep && NBC.isTimeDependent(iNBC, 1))
        [FGamma,tanMtxLoadPatch] = funcHandle ...
            (FGamma,BSplinePatches{1}, NBC.xiLoadExtension{iNBC}, ...
            NBC.etaLoadExtension{iNBC}, NBC.loadAmplitude{iNBC}, ...
            NBC.loadDirection{iNBC}, NBC.isFollower(iNBC,1), t, ...
            BSplinePatches{1}.int, '');
        if NBC.isFollower(iNBC, 1)
            tanMtxLoad{1} = tanMtxLoad{1} + tanMtxLoadPatch;
        end
    end

    %% 7iv.3. If the loading is not conservative add the contribution to the non-conservative load vector
    if NBC.isFollower(iNBC, 1)
        BSplinePatches{1}.FNonConservative = BSplinePatches{1}.FNonConservative + ...
            FGamma;
    end

    %% 7iv.4. Add The compute external load vector into the B-Spline array
    BSplinePatches{1}.FGamma = BSplinePatches{1}.FGamma +...
        FGamma;
end

%% 8. Compute the linear stiffness matrix of the membrane problem
if strcmp(outMsg, 'outputEnabled')
    message = 'Computing the linear stiffness matrix of the structure\n';
    fprintf([tab, '>>', ' ', message]);
end
stiffMtxLinear = computeTangentStiffMtxResVct ...
    (constStiff, tanMtxLoad, dHat, dHatSaved, dHatDot, dHatDotSaved, ...
    BSplinePatches, connections, propCoupling, loadFactor, noPatch, ...
    noTimeStep, counterNonlinearIterations, noWeakDBCCnd, ...
    propTransientAnalysis, isReferenceUpdated, strcat(tab, '\t'), outMsg);

%% 9. Compute the mass matrix of the membrane problem
if strcmp(outMsg, 'outputEnabled')
    message = 'Computing the mass matrix of the structure\n';
    fprintf([tab, '>>', ' ', message]);
end
massMtx = computeIGAMassMtxThinStructure(BSplinePatches, numDOFs);

%% 10. Solve the generalized eigenvalue problem to get the eigenmode shapes
if strcmp(outMsg,'outputEnabled')
    message = 'Solving the generalized eigenvalue problem to get the eigenmode shapes\n';
    fprintf([tab, '>>', ' ', message]);
end
[eigenmodeShapesDBC, naturalFrequencies] = ...
    eigs(stiffMtxLinear(freeDOFs, freeDOFs), massMtx(freeDOFs, freeDOFs), ...
    noEig, 'sm');
naturalFrequencies = naturalFrequencies*ones(length(naturalFrequencies), 1);
naturalFrequencies = sqrt(naturalFrequencies);
naturalFrequencies = naturalFrequencies./(2*pi);
[noModeShapes,m] = size(eigenmodeShapesDBC);

%% 11. Apply the boundary conditions onto the system
if strcmp(outMsg, 'outputEnabled')
    message = 'Applying the Dirichlet boundary conditions at each mode shape\n\n';
    fprintf([tab, '>>', ' ', message]);
end
eigenmodeShapes = zeros(noModeShapes, m);
for iMode = 1:noEig
    eigenmodeShapes(freeDOFs, iMode) = eigenmodeShapesDBC(:, iMode);
    eigenmodeShapes(homDOFs, iMode) = 0;
    eigenmodeShapes(inhomDOFs, iMode) = 0;
end

%% 12. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Modal analysis took %.2d seconds \n\n', computationalTime);
    fprintf('__________________Form Finding Analysis Ended__________________\n');
    fprintf('##############################################################\n\n\n');
end

end