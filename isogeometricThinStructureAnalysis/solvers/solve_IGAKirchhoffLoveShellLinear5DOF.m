function [dHat, FComplete, minElSize] = solve_IGAKirchhoffLoveShellLinear5DOF ...
    (BSplinePatch, solve_LinearSystem, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
% 
% Returns the displacement field and the complete force vector of a
% linear Kirchhoff-Love shell. Source Reference : 
% Modified version for 5 DOFs per control point.
%
% J. Kiendl "Isogeometric Analysis and Shape Optimal Design for Shell 
% Structures" Ph.D. Thesis, Technische Universtät München (2011)
%
%                Input :
%         BSplinePatch : Structure containing all the information regarding 
%                        the connections between the multipatches
%   solve_LinearSystem : Function handle to the linear equation system
%                        solver
%               outMsg : Whether or not to output message on refinement 
%                       progress
%                       'outputEnabled' : enables output information
%
%               Output :
%                 dHat : the displacement field
%            FComplete : the complete load vecto
%            minElArea : The minimum element area in the mesh
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
    fprintf('Static linear analysis for an isogeometric Kirchhoff-Love\n');
    fprintf('shell has been initiated (5 DOF version)\n');
    fprintf('_________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Define analysis
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

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

% Create an element freedom table for the patch in the array (5 DOFs per control point)
BSplinePatch.DOFNumbering = zeros(numCPs_xi, numCPs_eta, 5);
k = 1;
for cpj = 1:numCPs_eta
    for cpi = 1:numCPs_xi
        BSplinePatch.DOFNumbering(cpi, cpj, 1) = k;
        BSplinePatch.DOFNumbering(cpi, cpj, 2) = k + 1;
        BSplinePatch.DOFNumbering(cpi, cpj, 3) = k + 2;
        BSplinePatch.DOFNumbering(cpi, cpj, 4) = k + 3;
        BSplinePatch.DOFNumbering(cpi, cpj, 5) = k + 4;

        % Update counter (increment by 5 for 5 DOFs)
        k = k + 5;
    end
end

% Create the element freedom table for the BSplinePatch into the array of
% the patches (5 DOFs per control point)
BSplinePatch.EFTPatches = 1:5*BSplinePatch.noCPs;

% Place the B-Spline patch into an array
BSplinePatches = {BSplinePatch};

% Get number of DOFs (5 DOFs per control point)
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

% Initialize the displacement field
dHat = zeros(numDOFs, 1);

%% 1. Solve the linear system
[dHat, ~, ~, ~, FComplete, ~, ~, ~, ~, ~, minElSize] = ...
    solve_IGALinearSystem ...
    (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
    connections, dHat, dHatDot, dHatDDot, KConstant, massMtx, dampMtx, ...
    @computeStiffMtxAndLoadVctIGAKirchhoffLoveShellLinear5DOF, ...
    computeUpdatedGeometry, freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateDirichletBCs, masterDOFs, slaveDOFs, solve_LinearSystem, t, ...
    propCoupling, propTransientAnalysis, propNLinearAnalysis, propIDBC, ...
    plot_IGANLinear, isReferenceUpdated, isCosimulationWithEmpire, ...
    tab, graph, outMsg);

%% 2. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Static linear analysis took %.2d seconds \n\n', computationalTime);
    fprintf('______________Static Linear Analysis Ended_______________\n');
    fprintf('#########################################################\n\n\n');
end

end