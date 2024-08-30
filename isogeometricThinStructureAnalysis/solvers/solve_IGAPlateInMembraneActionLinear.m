function [dHat, FComplete, rankD, condK, minEig, minElSize] = ...
    solve_IGAPlateInMembraneActionLinear...
    (propAnalysis, BSplinePatch, stiffMtx, solve_LinearSystem, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement field  and the complete force vector 
% corresponding to 2D isogeometric plate in membrane action
%
%             Input :
%      propAnalysis : Structure containing general information on the 
%                     analysis,
%                         .type : Analysis type
%      BSplinePatch : Polynmial orders, knot vectors and control points of
%                     the B-Spline patch as well as the load vector, its
%                     technical parameters, Dirichlet boundary conditions 
%                     and integration rule
%                           .p,q : polynomial degrees
%                        .Xi,Eta : knot vectors
%                            .CP : vector containing control point 
%                                  coordinates and weights
%                       .isNURBS : Flag on whether the geometrical basis is 
%                                  NURBS or B-Spline
%                    .parameters : Technical and geometrical parameters
%                        .FGamma : force vector
%                       .homDOFs : Array containing information on 
%                                  homogeneous Dirichlet boundary 
%                                  conditions
%                     .inhomDOFs : Array containing information on the 
%                                  inhomogeneous Dirichlet boundary 
%                                  conditions
%               .valuesInhomDOFs : Values on the inhomogeneous Dirichlet 
%                                  boundary conditions
%                           .NBC : On the Neumann boundary conditions
%           stiffMtx : The precomputed global stiffness matrix
% solve_LinearSystem : Function handle to the linear equation system solver
%             outMsg : Whether or not to output message on refinement 
%                      progress 'outputEnabled' : enables output 
%                      information
%
%            Output :
%              dHat : the displacement field
%         FComplete : the complete load vector
%             rankD : The rank deficiency of the system
%             condK : The condition number of the system
%            minEig : The minimum eigenvalue of the system matrix
%         minElSize : The minimum element area size in the isogeometric%
%                     discretization
%
% Function layout :
%
% 0. Read input
%
% 1. Solve the linear equation system
%
% 2. Appendix
%
%% Function main body
if strcmp(outMsg, 'outputEnabled')
    fprintf('_________________________________________________________\n');
    fprintf('#########################################################\n');
    fprintf('Static linear analysis for an isogeometric plate in plane\n');
    fprintf('stress action has been initiated\n');
    fprintf('_________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Check input
numKnots_xi = length(BSplinePatch.Xi);
numKnots_eta = length(BSplinePatch.Eta);
numCPs_xi = length(BSplinePatch.CP(:, 1, 1));
numCPs_eta = length(BSplinePatch.CP(1, :, 1));
checkInputForBSplineSurface ...
    (BSplinePatch.p, numKnots_xi, numCPs_xi, BSplinePatch.q, ...
    numKnots_eta, numCPs_eta);

% Get the total number of isogeometric elements
numElmnts = ...
    length(unique(BSplinePatch.Xi))*length(unique(BSplinePatch.Eta));

% Initialize dummy variables
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
dHatDot = 'undefined';
dHatDDot = 'undefined';
propNLinearAnalysis = 'undefined';
connections = 'undefined';
propCoupling = 'undefined';
propGraph = 'undefined'; 
massMtx = 'undefined';
dampMtx = 'undefined';
updateDirichletBCs = 'undefined';
computeUpdatedGeometry = 'undefined';
propIDBC = 'undefined';
isReferenceUpdated = 'undefined';
isCosimulationWithEmpire = 'undefined';
plot_IGANLinear = 'undefined';

% Function handle to the computation of the linear stiffness matrix
if numElmnts == 1
    computeStiffMtxLoadVct = ...
        @computeStiffMtxAndLoadVctIGAPlateInMembraneActionLinear;
else
    computeStiffMtxLoadVct = ...
        @computeStiffMtxAndLoadVctIGAPlateInMembraneActionLinearPagewise;
end

% Steady-state analysis
propTransientAnalysis.isTimeDependent = false;
propTransientAnalysis.computeProblemMtrcsTransient = 'undefined';
propTransientAnalysis.isStaticStep = true;

% Get number of Control Points
numCPs = length(BSplinePatch.CP(:, 1, 1))*length(BSplinePatch.CP(1, :, 1));

% Get number of DOFs
numDOFs = 2*numCPs;
BSplinePatch.noDOFs = numDOFs;

% Constant part of the stiffness matrix
if isnumeric(stiffMtx)
    numDOFs_constant = length(stiffMtx);
    if numDOFs_constant ~= numDOFs
        error('The precomputed part of the stiffness matrix has %d DOFs but the system has %d DOFs', ...
            numDOFs_constant, numDOFs);
    else
        KConstant = stiffMtx;
    end
else
    KConstant = 'undefined';
end

% Initialize output array
dHat = zeros(numDOFs, 1);

% Get master and slave DOF
masterDOFs = BSplinePatch.masterDOFs;
slaveDOFs = BSplinePatch.slaveDOFs;

% Get the Dirichlet boundary conditions
homDOFs = BSplinePatch.homDOFs;
inhomDOFs = BSplinePatch.inhomDOFs;
valuesInhomDOFs = BSplinePatch.valuesInhomDOFs;

% Get a sequencial numbering of the unconstained DOFs into a vector
freeDOFs = 1:numDOFs;
freeDOFs(ismember(freeDOFs,homDOFs)) = [];

% The analysis is steady-state
t = 0;

% Collect the BSpline patches into an array
BSplinePatches = {BSplinePatch};

%% 1. Solve the linear equation system
[dHat, ~, ~, ~, FComplete, rankD, condK, minEig, ~, ~, minElSize] = ...
    solve_IGALinearSystem ...
    (propAnalysis, dHatSaved, dHatDotSaved, dHatDDotSaved, BSplinePatches, ...
    connections, dHat, dHatDot, dHatDDot, KConstant, massMtx, dampMtx, ...
    computeStiffMtxLoadVct, computeUpdatedGeometry, freeDOFs, homDOFs, ...
    inhomDOFs, valuesInhomDOFs, updateDirichletBCs, masterDOFs, slaveDOFs, ...
    solve_LinearSystem, t, propCoupling, propTransientAnalysis, ...
    propNLinearAnalysis, propIDBC, plot_IGANLinear, isReferenceUpdated, ...
    isCosimulationWithEmpire, '\t', propGraph, outMsg);

%% 2. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Static linear analysis took %.2d seconds \n\n', computationalTime);
    fprintf('______________Static Linear Analysis Ended_______________\n');
    fprintf('#########################################################\n\n\n');
end

end
