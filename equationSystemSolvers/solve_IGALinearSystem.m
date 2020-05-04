function [u, CPHistory, residual, isConverged, FComplete, rankD, ...
    condK, minEig, BSplinePatches, propCoupling, minElAreaSize] = ...
    solve_IGALinearSystem ...
    (propAnalysis, uSaved, uDotSaved, uDDotSaved, BSplinePatches, connections, ...
    u, uDot, uDDot, constMtx, massMtx, dampMtx, computeLinearMtrcsSteadyState, ...
    computeUpdatedGeometry, freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateDirichletBCs, masterDOFs, slaveDOFs, solve_LinearSystem, t, ...
    propCoupling, propTransientAnalysis, propNLinearAnalysis, propIDBC, ...
    plot_IGANLinear, isReferenceUpdated, isCosimulationWithEmpire, tab, ...
    propGraph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the solution to a linear system which correspond to the
% isogeometric discretization of the underlying field.
%
%                          Input :
%                   propAnalysis : Structure containing general information 
%                                  on the analysis,
%                                       .type : Analysis type
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
%  computeLinearMtrcsSteadyState : Function handle for the computation of 
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
%                                       .method : Time integration method
%                                    .alphaBeta : Bossak parameter
%                                        .gamma : Bossak parameter
%                                       .TStart : Start time of the 
%                                                 simulation
%                                         .TEnd : End time of the 
%                                                 simulation
%                                           .nT : Number of time steps
%                                           .dt : Time step
%                 .computeProblemMtrcsTransient : Function handle to the
%                                                 computation of the
%                                                 problem matrices
%                                                 according to the defined
%                                                 time integration scheme
%                                                 given the system matrix 
%                                                 of the steady-state
%                                                 system as well as the 
%                                                 mass matrix of the
%                                                 problem
%                            .computeUpdatedVct : Function handle to the 
%                                                 computation of the 
%                                                 updated values of the 
%                                                 discretized
%                                 .isStaticStep : Flag on whether the
%                                                 current step corresponds
%                                                 to the static one, namely
%                                                 before the time loop
%            propNLinearAnalysis : Dummy variable for this function
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
%                plot_IGANLinear : Dummy variable for this function
%             isReferenceUpdated : Flag on whether the reference geometry 
%                                  is updated
%       isCosimulationWithEmpire : Flag on whether co-simulation through
%                                  Empire is assumed
%                            tab : Tabulation for the output messages
%                      propGraph : Dummy variable for this function
%                         outMsg : On printing information during analysis 
%                                  on the command window
%
%                         Output :
%                              u : The converged discrete solution vector
%                      CPHistory : Dummy output for this function
%                       residual : Zero since it is a linear analysis
%                    isConverged : Dummy output for this function
%                      FComplete : The complete force vector of the system
%                          rankD : The rank deficiency of the linear system
%                          condK : The condition number of the linear 
%                                  system
%                         minEig : The minimum eigenvalue to the linear 
%                                  system
%                   propCoupling : Updated structure of the coupling
%                                  properties for multipatch system
%                  minElAreaSize : The minimum element size over the 
%                                  isogeometric mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Update the time-dependent inhomogeneous Dirichlet boundary conditions
%
% 2. Compute the linear matrices of the system
%
% 3. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
%
% 4. Solve the linear equation system
%
% 5. Re-assemble to the complete vector of unknowns
%
% 6. Compute the complete force vector
%
% 7. Appendix
%
%% Function main body

%% 0. Read input

% Dummy variables
CPHistory = 'undefined';
loadFactor = 'undefined';
tanMtxLoad = 'undefined';
noPatch = 'undefined';
noTimeStep = 'undefined';
iNLinearIter = 'undefined';
noWeakDBCCnd = 'undefined';
isReferenceUpdated = 'undefined';
isConverged = true;

% The residual for a linear analysis is zero
residual = 0;

% Check if the linear solver is called from a transient analysis
isTransient = false;
if ~ischar(propTransientAnalysis)
    if isfield(propTransientAnalysis, 'timeDependence')
        if strcmp(propTransientAnalysis.timeDependence, 'transient')
            isTransient = true;
        end
    end
end

%% 1. Update the time-dependent inhomogeneous Dirichlet boundary conditions
if ~ischar(propIDBC)
    if isa(updateDirichletBCs, 'function_handle')
        [homDOFs, inhomDOFs, valuesInhomDOFs] = ...
            updateDirichletBCs(BSplinePatches, homDOFs, propIDBC, t);
    else
        error('Variable updateInhomDOFs should define a function handle');
    end
end

%% 2. Compute the linear matrices of the system
if strcmp(outMsg, 'outputEnabled') && ~isTransient
    fprintf(strcat(tab, '>> Computing the system matrix and right hand side vector\n'));
end
[stiffMtx, RHS, BSplinePatches, propCoupling, minElAreaSize] = ...
    computeLinearMtrcsSteadyState ...
    (constMtx, tanMtxLoad, u, uSaved, uDot, uDotSaved, BSplinePatches, ...
    connections, propCoupling, loadFactor, noPatch, noTimeStep, ...
    iNLinearIter, noWeakDBCCnd, t, propTransientAnalysis, ...
    isReferenceUpdated, strcat(tab, '\t'), outMsg);
if isa(propTransientAnalysis.computeProblemMtrcsTransient, 'function_handle') && ...
        ~propTransientAnalysis.isStaticStep
    [stiffMtx, RHS] = propTransientAnalysis.computeProblemMtrcsTransient ...
        (u, uSaved, uDot, uDotSaved, uDDot, uDDotSaved, massMtx, dampMtx, ...
        stiffMtx, RHS, propTransientAnalysis);
end
if strcmp(outMsg, 'outputEnabled') && ~isTransient
    rankD = length(RHS(freeDOFs)) - rank(stiffMtx(freeDOFs, freeDOFs));
    if rankD ~= 0
        fprintf(strcat(tab,'>> The system has rank deficiency equal to %d\n'), rankD);
    end
    condK = cond(stiffMtx(freeDOFs, freeDOFs));
    fprintf(strcat(tab, '>> The condition number of the system is %d\n'), condK); 
else
    rankD = 'undefined';
    condK = 'undefined';
end
isMinEig = false;
if isstruct(propCoupling)
    if isfield(propCoupling, 'method')
        if ischar(propCoupling.method)
            if strcmp(propCoupling.method, 'nitsche')
                isMinEig = true;
            end
        end
    end
end
if isMinEig
    [~, eigVal] = eig(stiffMtx(freeDOFs, freeDOFs));
    eigVal(~eigVal) = Inf;
    minEig = min(min(eigVal));
    if strcmp(outMsg, 'outputEnabled')
        fprintf(strcat(tab, '>> The minimum eigenvalue of the system is %d\n'), minEig);
    end
else
    minEig = 'undefined';
end

%% 3. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
if norm(valuesInhomDOFs) ~= 0
    RHS = RHS - stiffMtx(:, inhomDOFs)*valuesInhomDOFs';
end

%% 4. Solve the linear equation system
if strcmp(outMsg, 'outputEnabled') && ~isTransient
    fprintf(strcat(tab, '>> Solving the linear system of %d equations\n'), length(RHS(freeDOFs)));
end
[uRed, hasLinearSystemConverged] = solve_LinearSystem ...
    (stiffMtx(freeDOFs, freeDOFs), RHS(freeDOFs), u(freeDOFs));
if ~hasLinearSystemConverged
    error('Linear equation solver has not converged');
end

%% 5. Re-assemble to the complete vector of unknowns
u(freeDOFs) = uRed;
u(homDOFs) = 0;
u(inhomDOFs) = valuesInhomDOFs;

%% 6. Compute the complete force vector
FComplete = stiffMtx*u;

%% 7. Appendix
if strcmp(outMsg,'outputEnabled') && ~isTransient
    fprintf('\n');
end

end
