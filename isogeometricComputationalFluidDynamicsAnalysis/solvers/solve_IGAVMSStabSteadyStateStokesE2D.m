function [up, F, minElSize] = solve_IGAVMSStabSteadyStateStokesE2D...
    (analysis, BSplinePatch, computeBodyForces, solve_LinearSystem, ...
    propIDBC, propNBC, propInt, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Solve the steady-state linear Stokes problem in 2D using the Bossak time 
% integration scheme for the temporal discretization and the NURBS 
% (isogeometric) basis for the spatial discretization as well as the
% Variational Multiscale method as a stabilization.
%
% Unconditional stability of the Bossak time integration scheme is ensured 
% if the following relations hold :
% 
% - alphaBeta <= .5 
% - beta >= gamma/2 >= .25
% - alphaBeta + gamma >= .25
%
%               Input :
%            analysis : Structure containing general information about the
%                       analysis
%                         .type : Analysis type
%        BSplinePatch : Structure containing computational information on
%                       the given B-Spline patch,
%                           p,q : Polynomial degrees along the xi-,eta-
%                                 parametric directions
%                        Xi,Eta : Knot vectors along the xi-,eta-
%                                 parametric directions
%                            CP : Control Point coordinates and weights
%                       isNURBS : Flag on whether the underlying basis is a 
%                                 B-Spline or a NURBS
%                  DOFNumbering : Numbering of the DOFs within the patch
%                       homDOFs : The global numbering of the DOFs where 
%                                 homogeneous Dirichlet boundary conditions 
%                                 are applied
%                     inhomDOFs : The global numbering of the DOFs where
%                                 inhomogeneous Dirichlet boundary 
%                                 conditions are applied
%                    parameters : Flow parameters
%                                   .nue : Dynamic viscosity
%   computeBodyForces : The body force vector [bx by]' or a function handle
%                       to its computation
%  solve_LinearSystem : Function handle to the linear equation system
%                       solver
%            propIDBC : Structure containing information on the 
%                       inhomogeneous Dirichlet boundary conditions,
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
%                                 .irb : Array .noConditions x m containing
%                                        the global numbering of the DOFs
%                                        which belong to each of the
%                                        segment
%                     .prescribedValue : Array of size .noConditions which
%                                        contains handles to functions
%                                        which determine the prescribed
%                                        values at each segment
%                          .isDominant : Flag on whether the inhomogeneous
%                                        Dirichlet boundary conditions are
%                                        dominant over the homogeneous ones
%                                        or not
%             propNBC : On the Neumann boundary conditions :
%                        .noConditions : Number of segements where those
%                                        conditions are applied
%                              .xiSpan : .noConditions x 2 array containing
%                                        the knot span extension of the
%                                        segments where those conditions
%                                        are applied in xi-direction
%                             .etaSpan : .noConditions x 2 array containing
%                                        the knot span extension of the
%                                        segments where those conditions
%                                        are applied in eta-direction
%                       .loadAmplitude : .noConditions x 2 array containing 
%                                         the rule for computing the load 
%                                         at each segment
%                       .loadDirection : .noConditions x 2 array containing
%                                        the direction of the load at each 
%                                        segment
%             propInt : On the spatial integration
%                               .type : 'default' or 'manual'
%                              .xiNGP : No. of GPs along xi-direction for 
%                                       stiffness entries
%                             .etaNGP : No. of GPs along eta-direction for 
%                                       stiffness entries
%                       .xiNGPForLoad : No. of GPs along xi-direction for 
%                                       load entries
%                      .etaNGPForLoad : No. of GPs along eta-direction for 
%                                       load entries
%                         .nGPForLoad : No. of GPs along boundary
%              outMsg : On printing information during analysis in the
%                       command window
%
%              Output :
%                  up : The solution vector of velocity and pressure
%                   F : The applied load vector of both boundary and body
%                       fluxes
%           minElSize : The minimum element area size over the isogeometric
%                       mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Find the prescribed and the free DOFs of the system
%
% 2. Create a structure representing a multipatch geometry with only one patch
%
% 3. Solve the steady-state problem
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_______________________________________________________________\n');
    fprintf('###############################################################\n');
    fprintf('Linear staedy-state analysis for isogeometric 2D incompressible \n');
    fprintf('Stokes flow problem has been initiated. \n\n');
    fprintf('_______________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Assign dummy variables
uSaved = 'undefined';
uDotSaved = 'undefined';
uDDotSaved = 'undefined';
uDDot = 'undefined';
connections = 'undefined';
uDot = 'undefined';
KConstant = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
computeUpdatedGeometry = 'undefined';
masterDOFs = 'undefined';
slaveDOFs = 'undefined';
propCoupling = 'undefined';
propNLinearAnalysis = 'undefined';
plot_IGANLinear = 'undefined';
isReferenceUpdated = 'undefined';
isCosimulationWithEmpire = 'undefined';
propGraph = 'undefined';

% Transient analysis properties
propFldDynamics.computeProblemMtrcsTransient = 'undefined';

% Function handle to the computation of the body forces
BSplinePatch.computeBodyForces = computeBodyForces;

% Function handle to the update of the inhomogeneous Dirichlet boundary
% conditions
updateDirichletBCs = @updateInhomDOFsIGAIncompressibleFlow;

% Global numbering of the DOFs where homogeneous and inhomogeneous
% Dirichlet boundary conditions are applied
homDOFs = BSplinePatch.homDOFs;
inhomDOFs = BSplinePatch.inhomDOFs;
valuesInhomDOFs = BSplinePatch.valuesInhomDOFs;

% Compute number of Control Points
numxi = length(BSplinePatch.CP(:, 1, 1));
numeta = length(BSplinePatch.CP(1, :, 1));
numCPs = numxi*numeta;

% Compute the number of degrees of freedom
numDOFs = 3*numCPs;
BSplinePatch.noDOFs = numDOFs;

% Assign a sequential numbering to the system DOFs
DOFNumbering = zeros(1, numDOFs);
for i = 1:numDOFs
    DOFNumbering(1, i) = i;
end

% Initialize force vector
F = zeros(numDOFs, 1);

% Initialize output array
up = zeros(numDOFs, 1);

% Steady-state analysis
t = 0;

% Tabulation for writting information in the command window
tab = '\t';

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs, inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs, prescribedDoFs)) = [];

%% 2. Create a structure representing a multipatch geometry with only one patch
BSplinePatches = {BSplinePatch};

%% 3. Solve the steady-state problem
[up, ~, ~, ~, ~, ~, ~, ~, ~, ~, minElSize] = ...
    solve_IGALinearSystem ...
    (analysis, uSaved, uDotSaved, uDDotSaved, BSplinePatches, ...
    connections, up, uDot, uDDot, KConstant, massMtx, dampMtx, ...
    @computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4StokesE2D, ...
    computeUpdatedGeometry, freeDOFs, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, updateDirichletBCs, masterDOFs, slaveDOFs, ...
    solve_LinearSystem, t, propCoupling, propFldDynamics, ...
    propNLinearAnalysis, propIDBC, plot_IGANLinear, isReferenceUpdated, ...
    isCosimulationWithEmpire, tab, propGraph, outMsg);

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Transient linear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('_______________Transient nonLinear Analysis Ended______________\n');
    fprintf('###############################################################\n\n\n');
end

end
