function [BSplinePatch, CPHistory, resHistory, isConverged, noIter] = ...
    solve_formFindingIGAMembrane ...
    (BSplinePatch,propFormFinding, solve_LinearSystem, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the B-Spline patch structure containing the updated values of the
% Control Points for the form-found geometry.
%
%           Input :
%    BSplinePatch : Structure containing information on the B-Spline patch
%                       .p,.q : The polynomial orders of the B-Spline patch
%                    .Xi,.Eta : The knot vectors in both parametric
%                               directions
%                         .CP : The set of Control Point coordinates and
%                               weights
%                    .isNURBS : Flag on whether the B-Spline patch is a
%                               NURBS or not
%                 .parameters : On the technical and geometrical parameters 
%                               of the B-Spline patch :
%                                   .prestress : On the prestress values :
%                                                 .voigtVector : The
%                                                                prestress
%                                                                values
%                                                                arranged
%                                                                in a Voigt
%                                                                vector
%             propFormFinding : Properties on the form finding analyse :
%                                .tolerance : Tolerance on the displacement 
%                                             difference between the 
%                                             different form finding
%                                             iterations
%                                .maxNoIter : Maximum number of iterations 
%                                             for convergence of the form 
%                                             finding iterations
%         solve_LinearSystem : Function handle to the equation system
%                              solver
%                     outMsg : Enables outputting information onto the
%                              command window when chosen 'outputEnabled'
%
%                      Output :
%                BSplinePatch : The updated with the new values of the
%                               Control Points B-Spline patch stucture
%                                   .CP : updated Control Point locations
%                                         corresponding to the form found
%                                         shape
%                              .weakDBC : Updated weak Dirichlet boundary
%                                         conditions array :
%                                     .automaticStabilization : History of
%                                                               the values
%                                                               for the
%                                                               stabiliza-
%                                                               tion
%                                                               factors
%                                                               through the
%                                                               iterations
%                                                               
%                   CPHistory : The history of the Control Points
%                               throughout the form finding iterations
%                  resHistory : The history of the error on the
%                               displacement field through the form-finding
%                               iterations
%                hasConverged : Flag on whether the form-finding iterations
%                               have converged up to the assigned tolerance 
%
% Function layout :
%
% 0. Read input
%
% 1. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
%
% 2. Compute the total number of DOFs for the B-Spline patch
%
% 3. Find the global DOF numbering of the DOFs where homogeneous, inhomogeneous conditions are applied
%
% 4. Create an element freedom table for the B-Spline patch
%
% 5. Find the total number of DOFs for the patch including the Lagrange Multipliers and make an EFT for each Lagrange Multipliers field employed
%
% 6. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions
%
% 7. Place the B-Spline patch into an array
%
% 8. Initialize Initialize the displacement field
%
% 9. Loop over all the form finding iterations
% ->
%    9i. Solve the static linear problem
%
%   9ii. Update the Control Point coordinates of the B-Spline patch
%
%  9iii. Compute the difference between the displacements from both iteration steps
%
%   9iv. Compute the residual and check convergence
%
%    9v. Save the displacement field of the previous form-finding step
%
%   9vi. Update form-finding iteration counter
% <-
%
% 10. Re-assign arrays
%
% 11. Appendix
%
%% Function main body              
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________________\n');
    fprintf('##################################################################\n');
    fprintf('Form Finding analysis for a membrane has been initiated \n\n');
    fprintf('Form finding iteration tolerance : %f \n', propFormFinding.tolerance);
    fprintf('Maximum number of form finding iterations: %d \n\n', propFormFinding.maxNoIter);
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
                    fprintf('Weak boundary conditions using the %s method are applied \n', BSplinePatch.weakDBC.method);
                    fprintf('over %d boundaries: \n', BSplinePatch.weakDBC.noCnd);
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
    fprintf('\n');
    fprintf('__________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Check input
if isfield(BSplinePatch.weakDBC,'computeTangMtxResVct') && ...
        (strcmp(BSplinePatch.weakDBC.method, 'penalty') || ...
        strcmp(BSplinePatch.weakDBC.method, 'lagrangeMultipliers'))
    error('No field weakDBC.%s needs to be defined for the method %s corresponding to the application of weak Dirichlet boundary conditions', ...
        BSplinePatch.weakDBC.computeTangMtxResVct, BSplinePatch.weakDBC.method);
end

% On the application of weak Dirichlet boundary conditions
isWeakDBC = false;
if ~isempty(BSplinePatch.weakDBC)
    if isfield(BSplinePatch.weakDBC, 'noCnd')
        if BSplinePatch.weakDBC.noCnd > 0
            isWeakDBC = true;
        end
    end
end

% Save the given parameters of the B-Spline patch
parametersSaved = BSplinePatch.parameters;

% Get the number of cables
noCables = BSplinePatch.cables.No;
parametersCablesSaved = struct([]);
for iCables = 1:noCables
    parametersCablesSaved{iCables} = BSplinePatch.cables.parameters{iCables};
end

% Save the Neumann boundary conditions
NBCSaved = BSplinePatch.NBC;

% Assign zero material to the B-Spline patch
BSplinePatch = rmfield(BSplinePatch, 'parameters');
BSplinePatch.parameters.E = 0.0;
BSplinePatch.parameters.nue = 0.0;
BSplinePatch.parameters.t = parametersSaved.t;
BSplinePatch.parameters.rho = 0.0;
BSplinePatch.parameters.prestress = parametersSaved.prestress;

% Assign zero material to the cables of the B-Spline patch
for iCables = 1:noCables
    BSplinePatch.cables.parameters{iCables}.E = 0.0;
    BSplinePatch.cables.parameters{iCables}.areaCS = parametersCablesSaved{iCables}.areaCS;
    BSplinePatch.cables.parameters{iCables}.rho = 0.0;
    BSplinePatch.cables.parameters{iCables}.prestress = parametersCablesSaved{iCables}.prestress;
end

% Assign follower load
for iCnd = 1:BSplinePatch.NBC.noCnd
    BSplinePatch.NBC.isFollower(iCnd,1) = false;
end

% Assign the the flag corresponding to whether the reference geometry is
% updated
isReferenceUpdated = true;

% Initialize convergence flag to false
isConverged = false;

% Initialize the history of the Control Points through the form finding
% iterations and the residual history
CPHistory = struct([]);
resHistory = zeros(propFormFinding.maxNoIter - 1, 1);

% Initialize form finding iteration counter
counterFormFindingIterations = 1;

% Initialize the dummy arrays
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
connections = 'undefined';
dHatDot = 'undefined';
dHatDDot = 'undefined';
propCoupling = 'undefined';
KConstant = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
graph = 'undefined';
updateDirichletBCs = 'undefined';
propIDBC = 'undefined';
masterDOFs = [];
slaveDOFs = [];

% Set flag for co-simulation with EMPIRE to false
isCosimulationWithEmpire = false;

% On plotting the deformed configuration through the nonlinear iterations
plot_IGANLinear = '';

% Static analysis
propTransientAnalysis.timeDependence = 'steadyState';
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;
t = 0;

% Nonlinear analysis
propNLinearAnalysis.method = 'newtonRapshon';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 1e-3;
propNLinearAnalysis.maxIter = 10;

% Function handle to the computation of the element tangent stiffness
% matrix and load vector
if BSplinePatch.noElmnts ~= 1
    computeTangentStiffMtxesVct = @computeTangentStiffMtxResVctIGAMembraneNLinear;
else
    computeTangentStiffMtxesVct = @computeTangentStiffMtxResVctIGAMembraneNLinearOutdated;
end

% Adjust tabulation
tab = '\t';

% Re-assign the arrays
CP = BSplinePatch.CP;

% Number of Control Points in xi,eta-direction
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));

% Create an element freedom table for the patch in the array
BSplinePatch.DOFNumbering = zeros(nxi, neta, 3);
k = 1;
for cpj = 1:neta
    for cpi = 1:nxi
        BSplinePatch.DOFNumbering(cpi, cpj, 1) = k;
        BSplinePatch.DOFNumbering(cpi, cpj, 2) = k + 1;
        BSplinePatch.DOFNumbering(cpi, cpj, 3) = k + 2;

        % Update counter
        k = k + 3;
    end
end

%% 1. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
BSplinePatch.noDOFs = 3*BSplinePatch.noCPs;
noDOFsPatchLM = 0;
if isWeakDBC
    if strcmp(BSplinePatch.weakDBC.method, 'lagrangeMultipliers')
        for iCnd = 1:BSplinePatch.weakDBC.noCnd
            noDOFsPatchLMCnd = 3*length(BSplinePatch.weakDBC.lambda{iCnd}.CP(:, 1));
            noDOFsPatchLM = noDOFsPatchLM + noDOFsPatchLMCnd;
            if iCnd == 1
                index = 3*BSplinePatch.noCPs;
            else
                index = BSplinePatch.weakDBC.lambda{iCnd - 1}.EFT(length(BSplinePatch.weakDBC.lambda{iCnd - 1}.EFT));
            end
            BSplinePatch.weakDBC.lambda{iCnd}.EFT = index + 1:index + noDOFsPatchLMCnd;
        end
    end
end
BSplinePatch.noDOFs = BSplinePatch.noDOFs + noDOFsPatchLM;

%% 2. Compute the total number of DOFs for the B-Spline patch
noDOFs = BSplinePatch.noDOFs;

%% 3. Find the global DOF numbering of the DOFs where homogeneous, inhomogeneous conditions are applied

% Find the numbering of the DOFs where homogeneous Dirichlet conditions are
% prescribed
homDOFs = BSplinePatch.homDOFs;

% Get the numbering and the values of the DOFs which are prescribed
inhomDOFs = BSplinePatch.inhomDOFs;
valuesInhomDOFs = BSplinePatch.valuesInhomDOFs;

% Find the numbering of the free DOFs
freeDOFs = zeros(noDOFs, 1);
for i = 1:noDOFs
    freeDOFs(i, 1) = i;
end
freeDOFs(ismember(freeDOFs, homDOFs)) = [];
freeDOFs(ismember(freeDOFs, inhomDOFs)) = [];

%% 4. Create an element freedom table for the B-Spline patch
BSplinePatch.EFTPatches = 1:noDOFs;

%% 5. Find the total number of DOFs for the patch including the Lagrange Multipliers and make an EFT for each Lagrange Multipliers field employed
if isWeakDBC
    if strcmp(BSplinePatch.weakDBC.method, 'lagrangeMultipliers')
        for iCnd = 1:BSplinePatch.weakDBC.noCnd
            % Get the number of Control Points
            noXiLambda = length(BSplinePatch.weakDBC.lambda{iCnd}.Xi);

            % Initialize the field of the DOF numbering
            BSplinePatch.weakDBC.lambda{iCnd}.DOFNumbering = zeros(noXiLambda, 3);

            % Compute the entries of the DOF numbering array
            k = 1;
            for cpi = 1:noXiLambda
                BSplinePatch.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 1) = k;
                BSplinePatch.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 2) = k + 1;
                BSplinePatch.weakDBC.lambda{iCnd}.DOFNumbering(cpi, 3) = k + 2;

                % Update counter
                k = k + 3;
            end
        end
    end
end

%% 6. Compute the constant matrices for each patch corresponding to the application of weak boundary conditions
if isWeakDBC
    isConstMtxWeakDBC = true;
    if isfield(BSplinePatch.weakDBC, 'method')
        if strcmp(BSplinePatch.weakDBC.method, 'penalty') || ...
                (strcmp(BSplinePatch.weakDBC.method, 'nitsche') && ...
                ~BSplinePatch.weakDBC.estimationStabilPrm)
            computeWeakDBCConstantProblemMatrices = @computeWeakDBCMtxPenaltyIGAMembrane;
        elseif strcmp(BSplinePatch.weakDBC.method, 'lagrangeMultipliers')
            computeWeakDBCConstantProblemMatrices = @computeWeakDBCMtxLagrangeMultipliersIGAMembrane;
        elseif strcmp(BSplinePatch.weakDBC.method, 'nitsche') && ...
                BSplinePatch.weakDBC.estimationStabilPrm
            isConstMtxWeakDBC = false;
        else
            error('Define a valid method for the application of weak Dirichlet boundary conditions in BSplinePatches{iPatches}.weakDBC.method');
        end
        if isConstMtxWeakDBC
            BSplinePatch.KConstant = ...
                computeWeakDBCConstantProblemMatrices ...
                (BSplinePatch, connections, ...
                BSplinePatch.noDOFs, propCoupling);
        else
            BSplinePatch.KConstant = 'undefined';
        end
    else
        error('Define variable method for the application of weak Dirichlet boundary conditions in BSplinePatches{iPatches}.weakDBC');
    end
else
    BSplinePatch.KConstant = 'undefined';
end

%% 7. Place the B-Spline patch into an array
BSplinePatches = {BSplinePatch};

%% 8. Initialize Initialize the displacement field
dHatPrevious = zeros(noDOFs,1);

%% 9. Loop over all the form finding iterations
if strcmp(outMsg,'outputEnabled')
    msgPNR = sprintf(strcat(tab,'\tLooping over the form finding iterations\n', ...
        tab, '\t----------------------------------------\n\n'));
    fprintf(msgPNR);
end
dHat_null = zeros(noDOFs, 1);
while ~isConverged && counterFormFindingIterations <= propFormFinding.maxNoIter
    %% 9i. Solve the static linear problem
    [dHat, ~, rH, ~, ~, ~, ~, ~, BSplinePatches, ~, ~] = ...
        solve_IGANLinearSystem...
        (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
        connections, dHat_null, dHatDot, dHatDDot, KConstant, massMtx, ...
        dampMtx, computeTangentStiffMtxesVct, ...
        @computeUpdatedGeometryIGAThinStructureMultipatches, freeDOFs, ...
        homDOFs, inhomDOFs, valuesInhomDOFs, updateDirichletBCs, masterDOFs, ...
        slaveDOFs, solve_LinearSystem, t, propCoupling, propTransientAnalysis, ...
        propNLinearAnalysis, propIDBC, plot_IGANLinear, isReferenceUpdated, ...
        isCosimulationWithEmpire, tab, graph, '');
    if length(find(rH)) > 2
        warning('More than 1 iterations needed for convergence');
    end
    
    %% 9ii. Update the Control Point coordinates of the B-Spline patch
    BSplinePatches{1}.CP = BSplinePatches{1}.CPd;
    CPHistory{counterFormFindingIterations} = BSplinePatches{1}.CP;
    
    %% 9iii. Compute the difference between the displacements from both iteration steps
	delta_dHat = dHat - dHatPrevious;
    
    %% 9iv. Compute the residual and check convergence
    resHistory(counterFormFindingIterations,1) = norm(delta_dHat);
    if strcmp(outMsg,'outputEnabled')
        msgNR = sprintf(strcat(tab,'\t||delta_dHat|| = %d at form-finding iteration No. %d \n'), ...
            resHistory(counterFormFindingIterations, 1), counterFormFindingIterations);
        fprintf(msgNR);
    end
    if resHistory(counterFormFindingIterations, 1) < propFormFinding.tolerance
        if strcmp(outMsg, 'outputEnabled')
            fprintf(strcat(tab, ' \tForm-finding iterations converged!\n\n'));
        end
        isConverged = true;
        break;
    end
    
    %% Debugging
%     clear graph;
%     graph.index = 1;
%     graph.postprocConfig = 'current';
%     graph.resultant = 'displacement';
%     graph.component = '2norm';
%     graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
%         (BSplinePatches,dHat,graph,'');
    
    %% 9v. Save the displacement field of the previous form-finding step
    dHatPrevious = dHat;
    
    %% 9vi. Update form-finding iteration counter
    counterFormFindingIterations = counterFormFindingIterations + 1;
end

%% 10. Re-assign arrays

% Check convergence
noIter = counterFormFindingIterations - 1;
if noIter == propFormFinding.maxNoIter
    if strcmp(outMsg, 'outputEnabled')
        warning(strcat(tab, ' \tForm-finding iterations did not converge up to tolerance %f!\n\n'), propFormFinding.tolerance);
    end
end

% Write ouput B-Spline patch array
BSplinePatch = BSplinePatches{1};

% Assign back the material properties of the patch
BSplinePatch = rmfield(BSplinePatch, 'parameters');
BSplinePatch.parameters.E = parametersSaved.E;
BSplinePatch.parameters.nue = parametersSaved.nue;
BSplinePatch.parameters.t = parametersSaved.t;
BSplinePatch.parameters.rho = parametersSaved.rho;
BSplinePatch.parameters.prestress = parametersSaved.prestress;

% Assign back the material properties of the cables
for iCables = 1:noCables
    BSplinePatch.cables.parameters{iCables}.E = parametersCablesSaved{iCables}.E;
    BSplinePatch.cables.parameters{iCables}.areaCS = parametersCablesSaved{iCables}.areaCS;
    BSplinePatch.cables.parameters{iCables}.rho = parametersCablesSaved{iCables}.rho;
    BSplinePatch.cables.parameters{iCables}.prestress = parametersCablesSaved{iCables}.prestress;
end

% % Assign back the variable to whether the load is consevative or not
for iCnd = 1:BSplinePatch.NBC.noCnd
    BSplinePatch.NBC.isFollower(iCnd,1) = NBCSaved.isFollower(iCnd, 1);
end

%% 11. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Form finding analysis took %.2d seconds \n\n', computationalTime);
    fprintf('___________________Form Finding Analysis Ended____________________\n');
    fprintf('##################################################################\n\n\n');
end

end
