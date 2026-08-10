function [dHat, FComplete, minElSize] = solve_IGAHSDTShellLinear ...
    (BSplinePatch, solve_LinearSystem, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
% 
% Returns the displacement field and the complete force vector of a
% linear HSDT shell with 5 DOF per control point: [u, v, w, θx, θy]
%
% Reference: Higher-order Shear Deformation Theory for shell analysis
%
%                Input :
%         BSplinePatch : Structure containing all the information regarding 
%                        the HSDT shell patch with 5 DOF per control point
%   solve_LinearSystem : Function handle to the linear equation system
%                        solver
%               outMsg : Whether or not to output message on analysis 
%                       progress
%                       'outputEnabled' : enables output information
%
%               Output :
%                 dHat : the displacement field [u, v, w, θx, θy] for all CPs
%            FComplete : the complete load vector
%             minElSize : The minimum element area in the mesh
%
% Function layout:
%
% 0. Read input
%
% 1. Solve the linear system
%
% 2. Appendix
%
%% Function main body
if strcmp(outMsg, 'outputEnabled')
    fprintf('_________________________________________________________\n');
    fprintf('#########################################################\n');
    fprintf('Static linear analysis for an isogeometric HSDT shell\n');
    fprintf('with 5 DOF per control point has been initiated\n');
    fprintf('DOF ordering: [u, v, w, θx, θy] per control point\n');
    fprintf('_________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Define analysis
analysis.type = 'isogeometricHSDTShellAnalysis';

% Assign dummy variables
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
connections = 'undefined';
dHatDot = 'undefined';
dHatDDot = 'undefined';
propCoupling = 'undefined';
propNLinearAnalysis = 'undefined';
KConstant = 'undefined';
computeUpdatedGeometry = 'undefined';
masterDOFs = 'undefined';
slaveDOFs = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
updateDirichletBCs = 'undefined';
propIDBC = 'undefined';
isReferenceUpdated = 'undefined';
isCosimulationWithEmpire = 'undefined';
graph = 'undefined';
plot_IGANLinear = 'undefined';

% The applied analysis is steady-state
propTransientAnalysis.timeDependence = 'steadyState';
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;
t = 0;

% Tabulation
tab = '\t';

% Re-assign the arrays
CP = BSplinePatch.CP;

% Number of Control Points in xi,eta-direction
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Create an element freedom table for the patch (5 DOF per control point)
BSplinePatch.DOFNumbering = zeros(numCPs_xi, numCPs_eta, 5);
k = 1;
for cpj = 1:numCPs_eta
    for cpi = 1:numCPs_xi
        BSplinePatch.DOFNumbering(cpi, cpj, 1) = k;     % u-displacement
        BSplinePatch.DOFNumbering(cpi, cpj, 2) = k + 1; % v-displacement
        BSplinePatch.DOFNumbering(cpi, cpj, 3) = k + 2; % w-displacement
        BSplinePatch.DOFNumbering(cpi, cpj, 4) = k + 3; % θx-rotation
        BSplinePatch.DOFNumbering(cpi, cpj, 5) = k + 4; % θy-rotation

        % Update counter (5 DOF per control point)
        k = k + 5;
    end
end

% Create the element freedom table for the BSplinePatch into the array of
% the patches (5 DOF per control point)
BSplinePatch.EFTPatches = 1:5*BSplinePatch.noCPs;

% Place the B-Spline patch into an array
BSplinePatches = {BSplinePatch};

% Get number of DOFs (5 DOF per control point)
numDOFs = 5*numCPs_xi*numCPs_eta;
BSplinePatches{1}.noDOFs = numDOFs;

% Find the numbering of the DOFs where homogeneous Dirichlet conditions are
% prescribed
homDOFs = BSplinePatch.homDOFs;

% Find the numbering of the free DOFs
freeDOFs = 1:numDOFs;
freeDOFs(ismember(freeDOFs,homDOFs)) = [];

% Get the numbering and the values of the DOFs which are prescribed
inhomDOFs = BSplinePatch.inhomDOFs;
valuesInhomDOFs = BSplinePatch.valuesInhomDOFs;

% Initialize the displacement field (5 DOF per control point)
dHat = zeros(numDOFs, 1);

%% 1. Solve the linear system
[dHat, ~, ~, ~, FComplete, ~, ~, ~, ~, ~, minElSize] = ...
    solve_IGALinearSystem ...
    (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
    connections, dHat, dHatDot, dHatDDot, KConstant, massMtx, dampMtx, ...
    @computeStiffMtxAndLoadVctIGAHSDTShellLinear, ...
    computeUpdatedGeometry, freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateDirichletBCs, masterDOFs, slaveDOFs, solve_LinearSystem, t, ...
    propCoupling, propTransientAnalysis, propNLinearAnalysis, propIDBC, ...
    plot_IGANLinear, isReferenceUpdated, isCosimulationWithEmpire, ...
    tab, graph, outMsg);

%% 2. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Static linear HSDT analysis took %.2d seconds \n\n', computationalTime);
    fprintf('Results contain 5 DOF per control point:\n');
    fprintf('  DOF 1,6,11,... = u-displacements\n');
    fprintf('  DOF 2,7,12,... = v-displacements\n');
    fprintf('  DOF 3,8,13,... = w-displacements\n');
    fprintf('  DOF 4,9,14,... = θx-rotations\n');
    fprintf('  DOF 5,10,15,... = θy-rotations\n\n');
    fprintf('______________Static HSDT Linear Analysis Ended_______________\n');
    fprintf('#############################################################\n\n\n');
end

end