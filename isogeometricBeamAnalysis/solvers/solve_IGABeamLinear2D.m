function [dHat, FComplete, minElSize] = ...
    solve_IGABeamLinear2D ...
    (analysis, p, Xi, CP, homDOFs, NBC, parameters, isNURBS, ...
    solve_LinearSystem, int, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement field corresponding to the 2D analysis of beam
% structures of the Bernoulli or the Timoshenko type within their
% isogeometric setting for the discretization.
%
%              Input :
%           analysis : Beam analysis type :
%                       'Bernoulli' : isogeometric Bernoulli beam analysis
%                      'Timoshenko' : isogeometric Timoshenko beam analysis
%                  p : The polynomial degree of the NURBS curve
%                 Xi : The knot vector
%                 CP : The set of Control Point coordinates and weights
%            homDOFs : The global numbering of the constrained DoFs
%                NBC : On the Neumann boundary conditions
%                               .noCnd : Number of conditions
%                      .computeLoadVct : Structure array containing the 
%                                        function name for the computation 
%                                        of the load vector
%                     .xiLoadExtension : Structure containing the extensions 
%                                        for the load application
%                       .loadAmplitude : Array containing the load 
%                                        amplitudes
%                       .loadDirection : Array containing the load 
%                                        directions
%         parameters : The technical parameters of the beam
%            isNURBS : Flag on whether the geometrical basis is a NURBS or a 
%                      B-Spline
% solve_LinearSystem : Function handle to the solution of the resulting
%                      linear equation system
%                int : The chosen intergration rule
%             outMsg : 'outputEnabled' : enables output information on the
%                      command window
%
%             Output :
%               dHat : The vector of the Control Point variables
%          FComplete : The complete force vector
%          minElSize : The minimum element size in the isogeometric mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Assign DoF numbering as CP: CP1->dof1,dof2,dof3 CP2->dof4,dof5,dof6
%
% 2. Solve the linear equation system
%
% 3. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('______________________________________________________\n');
    fprintf('######################################################\n');
    fprintf('Static Linear Analysis for the isogeometric ');
    if strcmp(analysis.type, 'Bernoulli')
        fprintf('Bernoulli ');
    elseif strcmp(analysis.type, 'Timoshenko')
        fprintf('Timoshenko ');
    end
    fprintf('\nbeam has been initiated\n');
    fprintf('______________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Create the B-Spline patch
BSplinePatch.p = p;
BSplinePatch.Xi = Xi;
BSplinePatch.CP = CP;
BSplinePatch.isNURBS = isNURBS;
BSplinePatch.int = int;
BSplinePatch.parameters = parameters;
BSplinePatch.NBC = NBC;

% On the output messaging
tab = '\t';

% Compute the number of Control Points
numCPs_xi = length(CP(:, 1));

% Compute the number of DOFs
if strcmp(analysis.type, 'Bernoulli')
    numDOFs = 2*numCPs_xi;
elseif strcmp(analysis.type, 'Timoshenko')
    numDOFs = 3*numCPs_xi;
end

% Assign a sequential numbering to the system DOFs
DOFSequentialNumbering = zeros(1, numDOFs);
for i = 1:numDOFs
    DOFSequentialNumbering(1, i) = i;
end

% Find the free DOFs of the system
freeDOFs = DOFSequentialNumbering;
freeDOFs(ismember(freeDOFs, homDOFs)) = [];

% Create the dummy variables
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
dHatDot = 'undefined';
dHatDDot = 'undefined';
propNLinearAnalysis = 'undefined';
connections = 'undefined';
propCoupling = 'undefined';
KConstant = 'undefined';
masterDOFS = 'undefined';
slaveDOFs = 'undefined';
computeUpdatedGeometry = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
updateDirichletBCs = 'undefined';
propIDBC = 'undefined';
plot_IGANLinear = 'undefined';
isReferenceUpdated = 'undefined';
isCosimulationWithEmpire = 'undefined';
propGraph = 'undefined';

% Static linear analysis
propTransientAnalysis.timeDependence = 'steadyState';
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;
t = 0;

% On the inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Get the rule for computating the master stiffness matrix according to the
% analysis type
if strcmp(analysis.type, 'Bernoulli')
    computeLinearMatricesSteadyState = @computeStiffMtxIGABernoulliBeamLinear;
elseif strcmp(analysis.type, 'Timoshenko')
    computeLinearMatricesSteadyState = @computeStiffMtxIGATimoshenkoBeamLinear;
end

% Initialize the output array
dHat = zeros(numDOFs, 1);

%% 1. Assign DoF numbering as CP: CP1->dof1,dof2,dof3 CP2->dof4,dof5,dof6

% Initialize the array of the degrees of freedom
if strcmp(analysis.type, 'Bernoulli')
    DOFNumbering = zeros(numCPs_xi, 2);
elseif strcmp(analysis.type, 'Timoshenko')
    DOFNumbering = zeros(numCPs_xi, 3);
end

% Initialize counter
k = 1;

% Loop over all the Control points
if strcmp(analysis.type, 'Bernoulli')
    for cpi = 1:numCPs_xi
        DOFNumbering(cpi, 1) = k;
        DOFNumbering(cpi, 2) = k + 1;

        % Update counter
        k = k + 2;
    end
elseif strcmp(analysis.type, 'Timoshenko')
    for cpi = 1:numCPs_xi
        DOFNumbering(cpi, 1) = k;
        DOFNumbering(cpi, 2) = k + 1;
        DOFNumbering(cpi, 3) = k + 2;

        % Update counter
        k = k + 3;
    end
end

% Attribute the DOF numbering to the B-Spline patches themselves
BSplinePatch.DOFNumbering = DOFNumbering;

% Collect all the B-Spline patches into an array
BSplinePatches = {BSplinePatch};

%% 2. Solve the linear equation system
[dHat, ~, ~, ~, FComplete, rankD, condK, minEig, ~, ~, minElSize] = ...
    solve_IGALinearSystem ...
    (analysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
    connections, dHat, dHatDot, dHatDDot, KConstant, massMtx, dampMtx, ...
    computeLinearMatricesSteadyState, computeUpdatedGeometry, freeDOFs, ...
    homDOFs, inhomDOFs, valuesInhomDOFs, updateDirichletBCs, ...
    masterDOFS, slaveDOFs, solve_LinearSystem, t, propCoupling, ...
    propTransientAnalysis, propNLinearAnalysis, propIDBC, ...
    plot_IGANLinear, isReferenceUpdated, isCosimulationWithEmpire, ...
    tab, propGraph, outMsg);

%% 3. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Static Linear Analysis took %d seconds \n\n',computationalTime);
    fprintf('_____________Static Linear Analysis Ended_____________\n');
    fprintf('######################################################\n\n\n');
end

end

