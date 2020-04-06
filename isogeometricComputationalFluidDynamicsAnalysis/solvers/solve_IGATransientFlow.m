function [upHistory, resHistory, minElSize] = ...
    solve_IGATransientFlow ...
    (analysis, BSplinePatch, computeBodyForces, computeInitCnds, ...
    computesSteadyStateMatrices, solve_LinearSystem, ...
    solve_IGAEquationSystem, propIDBC, propFldDynamics, ...
    propNLinearAnalysis, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Solve the transient nonlinear Navier-Stokes problem in 2D using the 
% Bossak time integration scheme for the temporal discretization and the 
% NURBS (isogeometric) basis for the spatial discretization.
%
% Unconditional stability of the Bossak time integration scheme is ensured 
% if the following relations hold :
% 
% - alphaBeta <= .5 
% - beta >= gamma/2 >= .25
% - alphaBeta + gamma >= .25
%
%                       Input :
%                    analysis : Structure containing general information 
%                               about the analysis,
%                                  .type : Analysis type
%                BSplinePatch : Structure containing information about the
%                               computational B-Spline patch,
%                                  .p,.q : Polynomial degrees along xi- and
%                                          eta- parametric directions
%                               .Xi,.Eta : Knots vectors along xi- and
%                                          eta- parametric directions
%                                    .CP : Control Point coordinates and
%                                          weights
%                               .isNURBS : Flag on whether the basis is a
%                                          NURBS or a B-Spline
%                               .homDOFs : Global numbering of the DOFs 
%                                          where homogeneous Dirichlet 
%                                          boundary conditions are applied
%                             .inhomDOFs : Global numbering of the DOFs 
%                                          where inhomogeneous Dirichlet 
%                                          boundary conditions are applied
%                       .valuesInhomDOFs : Prescribed values of the DOFs 
%                                          where inhomogeneous Dirichlet 
%                                          boundary conditions are applied. 
%                                          This array can be empty and get 
%                                          dynamically allocated for the 
%                                          transient inhomogeneous 
%                                          Dirichlet boundary conditions
%           computeBodyForces : The body force vector [bx by]' or function 
%                               handle to its computation
%             computeInitCnds : Function handle to the computation of the
%                               initial conditions
% computesSteadyStateMatrices : Function handle to the computation of the 
%                               problem matrices corresponding to the 
%                               steady-state problem
%          solve_LinearSystem : Function handle to the linear equation 
%                               system solver
%     solve_IGAEquationSystem : Function handle to the solution of the
%                               equation system resulting out of the
%                               isogeometric discretization
%                    propIDBC : On the inhomogeneous Dirichlet boundary 
%                               conditions,
%                                  .noCnd : Number of segements where those
%                                           conditions are applied
%                                 .xiSpan : .noConditions x 2 array 
%                                           containing the knot span 
%                                           extension of the segments where 
%                                           those conditions are applied in 
%                                           xi-direction
%                                .etaSpan : .noConditions x 2 array 
%                                           containing the knot span 
%                                           extension of the segments where 
%                                           those conditions are applied in 
%                                           eta-direction
%                    .prescribedDirection : .noConditions x 1 array 
%                                           containing the direction of the 
%                                           load application for each 
%                                           condition
%                     .isUniqueOnBoundary : .noConditions x 1 array 
%                                           containing flags on whether 
%                                           each of those conditions are 
%                                           unique over their application 
%                                           boundary
%                        .prescribedValue : Array of size .noConditions 
%                                           which contains handles to 
%                                           functions which determine the 
%                                           prescribed values at each 
%                                           segment
%                             .isDominant : Flag on whether the 
%                                           inhomogeneous Dirichlet 
%                                           boundary conditions are
%                                           dominant over the homogeneous 
%                                           ones or not
%                     propNBC : On the Neumann boundary conditions,
%                                      .numCnd : Number of segements where 
%                                                those conditions are 
%                                                applied
%                                      .xiSpan : .numCnd x 2 array 
%                                                containing the knot span 
%                                                extension of the segments 
%                                                where those conditions are 
%                                                applied in xi-direction
%                                     .etaSpan : .numCnd x 2 array 
%                                                containing the knot span 
%                                                extension of the segments 
%                                                where those conditions are 
%                                                applied in eta-direction
%                               .loadAmplitude : .numCnd x 2 array 
%                                                containing the rule for 
%                                                computing the load at each 
%                                                segment
%                               .loadDirection : .noConditions x 2 array 
%                                                containing the direction 
%                                                of the load at each 
%                                                segment
%             computeInitCnds : Function handle to the computation of the 
%                               initial boundary conditions computation
%             propFldDynamics : Structure containing information on the
%                               transient analysis,
%                                 .method : The time integration method
%                              .alphaBeta : (parameter for Bossak scheme)
%                                  .gamma : (parameter for Bossak scheme)
%                                 .TStart : Start time of the simulation
%                                   .TEnd : End time of the simulation
%                            .noTimeSteps : Number of time steps
%                                     .dt : Time step
%         propNLinearAnalysis : Structure containing information on the non
%                               -linear analysis,
%                                  .method : The nonlinear solution scheme
%                                     .eps : The residual tolerance
%                                 .maxIter : The maximum number of 
%                                            nonlinear iterations
%                      outMsg : On printing information during analysis in 
%                               the command window
%
%                      Output :
%                   upHistory : The history of the velocity and pressure 
%                               field throughout the transient analysis
%                  resHistory : The history of the residual vector
%                   minElSize : The minimum element area size over the 
%                               isogeometric discretization
%
% Function layout :
%
% 0. Read input
%
% 1. Remove DOFs which belong to two different boundary conditions
%
% 2. Find the prescribed and the free DOFs of the system
%
% 3. Solve the transient system applying the Bossak time integration scheme
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________\n');
    fprintf('###################################################################\n');
    fprintf('Nonlinear transient analysis for isogeometric 2D incompressible \n');
    fprintf('Navier-Stokes flow problem using the Bossak time integration scheme \n');
    fprintf('has been initiated. \n\n');
    if strcmp(propNLinearAnalysis.method, 'Newton')
        fprintf('Nonlinear scheme : Newton method \n');
    end
    fprintf('Residual tolerance: %d \n', propNLinearAnalysis.eps);
    fprintf('Maximum number of iterations: %d \n', propNLinearAnalysis.maxIter);
    fprintf('Maximum number of nonlinear iterations = %d \n', propNLinearAnalysis.maxIter);
    if propFldDynamics.isAdaptive
        fprintf('Adaptive time step chosen \n');
    else
        fprintf('Constant time step: %d (seconds) \n', propFldDynamics.dt);
    end
    fprintf('___________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Assign dummy variables
connections = 'undefined';
masterDOFs = 'undefined';
slaveDOFs = 'undefined';
nodesALE = 'undefined';
computeConstantMatrices = 'undefined';
computeUpdatedMesh = 'undefined';
computeUpdatedGeometry = 'undefined';
propCoupling = 'undefined';
propPostproc = 'undefined';
caseName = 'undefined';
pathToOutput = 'undefined';
propOutput = 'undefined';
isReferenceUpdated = 'undefined';

% Function handle to the computation of the mass matrix
computeMassMtx = @computeMassMtx4IGAVMSStabNSE2D;

% Function handle to the update of the prescribed values for the
% inhomogeneous Dirichlet boundary conditions
updateInhomDOFs = @updateInhomDOFsIGAIncompressibleFlow;

% Initialize the array of the values imposed on the DOFs with inhomogeneous
% Dirichlet boundary conditions
valuesInhomDOFs = [];

% Function handle to the computation of the body forces
BSplinePatch.computeBodyForces = computeBodyForces;

% Global numbering of the DOFs where homogeneous and inhomogeneous
% Dirichlet boundary conditions are applied
homDOFs = BSplinePatch.homDOFs;
inhomDOFs = BSplinePatch.inhomDOFs;

% Empire co-simulation
propEmpireCoSimulation.isCoSimulation = false;

% Title
title = 'Transient isogeometric incompressible flow analysis';

% Compute number of Control Points
numCPs_xi = length(BSplinePatch.CP(:, 1, 1));
numCPs_eta = length(BSplinePatch.CP(1, :, 1));
numCPs = numCPs_xi*numCPs_eta;

% Compute the number of degrees of freedom
numDOFs = 3*numCPs;
BSplinePatch.noDOFs = numDOFs;

% Tabulation for writting information in the command window
tab = '\t';

% Insert the B-Spline patch into an array of patches
BSplinePatches = {BSplinePatch};

%% 1. Remove DOFs which belong to two different boundary conditions
DOFsInBothHomAndInhomDBC = [];
if propIDBC.isDominant
    for iIDBC = 1:length(inhomDOFs)
        index = find(homDOFs == inhomDOFs(iIDBC));
        if norm(index) ~= 0
            DOFsInBothHomAndInhomDBC = ...
                mergesorted(DOFsInBothHomAndInhomDBC, index);
        end
    end
    homDOFs(DOFsInBothHomAndInhomDBC) = [];
else
    for iDBC = 1:length(homDOFs)
        index = find(inhomDOFs == homDOFs(iDBC));
        if norm(index) ~= 0
            DOFsInBothHomAndInhomDBC = ...
                mergesorted(DOFsInBothHomAndInhomDBC, index);
        end
    end
    inhomDOFs(DOFsInBothHomAndInhomDBC) = [];
end

%% 2. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs, inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = 1:numDOFs;
freeDOFs(ismember(freeDOFs,prescribedDoFs)) = [];

%% 3. Solve the transient system applying the Bossak time integration scheme
[upHistory, resHistory, ~, ~, minElSize] = ...
    solve_IGATransientAnalysis ...
    (analysis, BSplinePatches, connections, freeDOFs, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, updateInhomDOFs, masterDOFs, slaveDOFs, nodesALE, ...
    computeInitCnds, solve_IGAEquationSystem, computeConstantMatrices, ...
    computeMassMtx, computesSteadyStateMatrices, computeUpdatedMesh, ...
    computeUpdatedGeometry, solve_LinearSystem, ...
    propCoupling, propFldDynamics, propNLinearAnalysis, propPostproc, ...
    propIDBC,caseName, pathToOutput, title, propOutput, ...
    isReferenceUpdated, propEmpireCoSimulation, tab, outMsg);

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Transient nonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('________________Transient nonLinear Analysis Ended_________________\n');
    fprintf('###################################################################\n\n\n');
end

end