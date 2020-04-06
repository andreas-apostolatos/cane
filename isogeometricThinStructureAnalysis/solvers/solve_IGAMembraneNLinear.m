function [dHat, CPHistory, resHistory, isConverged, BSplinePatches, minElASize] = ...
    solve_IGAMembraneNLinear...
    (BSplinePatch, propNLinearAnalysis, solve_LinearSystem, plot_IGANLinear, ...
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
% Returns the displacement field corresponding to a geometrically
% non-linear anslysis for the isogeometric membrane.
%
%                Input : 
%         BSplinePatch : Structure containing all the information regarding 
%                        the connections between the multipatches
%  propNLinearAnalysis : Structure on the non-linear analysis settings :
%                                    #Load Steps/Stages#
%                              .noLoadSteps : Selected number of load steps 
%                                             (stages of the loading)
%                                .tolerance : Tolerance for the Newton 
%                                             iterations on the 2-norm
%                                  .maxIter : Maximum number of the Newton 
%                                             iteration
%   solve_LinearSystem : Function handle to the linear equation system
%                        solver
%      plot_IGANLinear : Function handle to plotting the current
%                        configuration together with resultants in case the
%                        variable is not an empty string ''
%                graph : On the graphics
%               outMsg : Whether or not to output message on refinement 
%                       progress
%                       'outputEnabled' : enables output information
%
%               Output :
%                 dHat : The displacement field of the patch
%            CPHistory : The deformation history of the Control Points
%           resHistory : The residual history throughout the nonlinear
%                        iterations
%          isConverged : Flag on whether the nonlinear iterations has
%                        converged or not
%       BSplinePatches : The updated array of the B-Spline patches with the
%                        estimated stabilization parameters
%           minElASize : The minimum element area size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the constant problem matrix
%
% 2. Solve the nonlinear system
%
% 3. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_________________________________________________________________________\n');
    fprintf('#########################################################################\n');
    fprintf('Static nonlinear analysis for an isogeometric membrane has been initiated \n\n');
    if isfield(BSplinePatch, 'weakDBC')
        if isfield(BSplinePatch.weakDBC, 'noCnd')
            if BSplinePatch.weakDBC.noCnd > 0
                fprintf('Weak boundary conditions using the %s method are applied \n', ...
                    BSplinePatch.weakDBC.method);
                fprintf('over %d boundaries of the B-Spline patch: \n', ...
                    BSplinePatch.weakDBC.noCnd);
                if strcmp(BSplinePatch.weakDBC.method, 'Nitsche')
                    if BSplinePatch.weakDBC.estimationStabilPrm == true
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
    fprintf('Nonlinear scheme : Newton method \n');
    fprintf('Number of load steps = %d \n', propNLinearAnalysis.noLoadSteps);
    fprintf('Residual tolerance = %d \n', propNLinearAnalysis.eps);
    fprintf('Maximum number of nonlinear iterations = %d \n', propNLinearAnalysis.maxIter);
    fprintf('__________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Analysis type
analysis.type = 'isogeometricMembraneNonlinear';

% Initialize the dummy arrays
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
connections = 'undefined';
dHatDot = 'undefined';
dHatDDot = 'undefined';
propCoupling = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
updateDirichletBCs = 'undefined';
propIDBC = 'undefined';
masterDOFs = [];
slaveDOFs = [];

% Set the flag for the co-simulation with EMPIRE to false
isCosimulationWithEmpire = false;

% Flag on whether the reference configuration is updated
isReferenceUpdated = false;

% Static analysis
propTransientAnalysis.timeDependence = 'steadyState';
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;
t = 0;

% Compute constant problem matrices in case of the application of weak
% boundary conditions
if isfield(BSplinePatch.weakDBC, 'computeConstMtx')
    computeConstantProblemMatrices = str2func(BSplinePatch.weakDBC.computeConstMtx);
else
     computeConstantProblemMatrices = 'undefined';
end

% Function handle to the computation of the element tangent stiffness
% matrix and load vector
if BSplinePatch.noElmnts ~= 1
    computeTangentStiffMtxesVct = @computeTangentStiffMtxResVctIGAMembraneNLinear;
else
    computeTangentStiffMtxesVct = @computeTangentStiffMtxResVctIGAMembraneNLinearOutdated;
end
% computeTangentStiffMtxesVct = @computeTangentStiffMtxResVctIGAMembraneNLinearOutdated;

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

% Create the element freedom table for the BSplinePatch into the array of
% the patches
BSplinePatch.EFTPatches = 1:3*BSplinePatch.noCPs;

% Place the B-Spline patch into an array
BSplinePatches = {BSplinePatch};

% Get number of DOFs
numDOFs = 3*nxi*neta;
BSplinePatches{1}.noDOFs = numDOFs;

% Find the numbering of the DOFs where homogeneous Dirichlet conditions are
% prescribed
homDOFs = BSplinePatch.homDOFs;

% Find the numbering of the free DOFs
freeDOFs = zeros(numDOFs, 1);
for i = 1:numDOFs
    freeDOFs(i, 1) = i;
end
freeDOFs(ismember(freeDOFs, homDOFs)) = [];

% Get the numbering and the values of the DOFs which are prescribed
inhomDOFs = BSplinePatch.inhomDOFs;
valuesInhomDOFs = BSplinePatch.valuesInhomDOFs;

% Initialize the displacement field
dHat = zeros(numDOFs, 1);

%% 1. Compute the constant problem matrix
if isa(computeConstantProblemMatrices, 'function_handle')
    KConstant = computeConstantProblemMatrices ...
        (BSplinePatches, connections, numDOFs, propCoupling);
else
    KConstant = 'undefined';
end

%% 2. Solve the nonlinear system
[dHat, CPHistory, resHistory, isConverged, ~, ~, ~, ~, BSplinePatches, ~, ...
    minElASize] = solve_IGANLinearSystem ...
    (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
    connections, dHat, dHatDot, dHatDDot, KConstant, massMtx, dampMtx, ...
    computeTangentStiffMtxesVct, ...
    @computeUpdatedGeometryIGAThinStructureMultipatches, freeDOFs, ...
    homDOFs, inhomDOFs, valuesInhomDOFs, updateDirichletBCs, masterDOFs, ...
    slaveDOFs, solve_LinearSystem, t, propCoupling, propTransientAnalysis, ...
    propNLinearAnalysis, propIDBC, plot_IGANLinear, isReferenceUpdated, ...
    isCosimulationWithEmpire, tab, graph, outMsg);

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Static nonlinear analysis took %.2d seconds \n\n', computationalTime);
    fprintf('______________________Static Linear Analysis Ended_______________________\n');
    fprintf('#########################################################################\n\n\n');
end

end
