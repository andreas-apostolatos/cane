function [u, FComplete, hasConverged, minElSize] = solve_FEMLinearSystem...
    (propAnalysis, uSaved, uDotSaved, uDDotSaved, mesh, F, computeBodyForces, ...
    parameters, u, uDot, uDDot, massMtx, dampMtx, precompStiffMtx, ...
    precomResVct, computeProblemMatricesSteadyState, DOFNumbering, ...
    freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, uMeshALE, ...
    solve_LinearSystem, propTransientAnalysis, propNLinearAnalysis, int, ...
    tab, outMsg)
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
% classical finite element discretization of the underlying field.
%
%                             Input :
%                      propAnalysis : Structure containing general 
%                                     information about the analysis,
%                                       .type : The analysis type
%                            uSaved : The discrete solution field of the 
%                                     previous time step
%                         uDotSaved : The time derivative of the discrete 
%                                     solution field of the previous time 
%                                     step
%                        uDDotSaved : The second order time derivative of 
%                                     the discrete solution field of the 
%                                     previous time (dummy variable for 
%                                     this function)
%                              mesh : Nodes and elements of the underlying 
%                                     mesh
%                                 F : The boundary force vector
%                 computeBodyForces : The body force vector
%                        parameters : The parameters of the physical field
%                                 u : Initial guess for the primary field 
%                                     (just an emtpy array must be given to 
%                                      this function)
%                              uDot : Initial guess for the time derivative 
%                                     of the primary field (dummy input for 
%                                     this function)
%                           massMtx : Mass matrix
%                           dampMtx : Damping matrix
%                   precompStiffMtx : Precomputed part of the stiffness 
%                                     matrix
%                      precomResVct : Precomputed part of the residual 
%                                     vector
% computeProblemMatricesSteadyState : Function handle for the computation 
%                                     of the matrices and vectors nessecary 
%                                     for computing the solution vector
%                                     corresponding to the steady-state
%                      DOFNumbering : The global numbering of the DOFs in a
%                                     3-dimensional array
%                          freeDOFs : The global numbering of the 
%                                     uncostrained DOFs
%                           homDOFs : The global numbering of the DOFs 
%                                     where homogeneous Dirichlet boundary 
%                                     conditions are applied
%                         inhomDOFs : The global numbering of the DOFs 
%                                     where inhomogeneous Dirichlet 
%                                     boundary conditions are applied
%                   valuesInhomDOFs : The vector containing the magnitude 
%                                     of the prescribed values on the DOFs 
%                                     with inhomogeneous Dirichlet boundary 
%                                     conditions
%                          uMeshALE : ALE mesh velocity (dummy variable for 
%                                     this function)
%                solve_LinearSystem : Function handle to the computation of 
%                                     the solution of a linear equation 
%                                     system
%             propTransientAnalysis : Transient analysis parameters:
%                                         .method : Time integration method
%                                      .alphaBeta : Bossak parameter
%                                          .gamma : Bossak parameter
%                                             .T0 : Start time of the 
%                                                   simulation
%                                           .TEnd : End time of the 
%                                                   simulation
%                                             .nT : Number of time steps
%                                             .dt : Time step (numeric or 
%                                                   adaptive)
%               propNLinearAnalysis : Nonlinear analysis parameters (dummy 
%                                     input variable for this function)
%                                       .scheme : The nonlinear solution 
%                                                 scheme
%                                    .tolerance : The residual tolerance
%                                      .maxIter : The maximum number of 
%                                                 nonlinear iterations
%                               int : On the spatial integration
%                                           .type : 'default' or 'manual'
%                                  .parentElement : 'tri', 'quad', usw.
%                                         .noXiGP : if .parentElement == 
%                                                   'quad'
%                                        .noEtaGP : if .parentElement == 
%                                                   'quad'
%                                      .noGP : if .parentElement == 'tri'
%                               tab : Tabulation for the output in the 
%                                     command window
%                            outMsg : On printing information during 
%                                     analysis in the command window
%
%                            Output :
%                                 u : The converged discrete solution 
%                                     vector
%                         FComplete : The complete flux vector of the 
%                                     system
%                      hasConverged : Flag on whether the nonlinear 
%                                     iterations have converged (for linear 
%                                     computations it is returned always 
%                                     true)
%                         minElSize : The minimum element size over the 
%                                     isogeometric mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the linear matrices of the steady-state problem
%
% 2. Compute the stiffness matrix and the residual vector for the transient problem
%
% 3. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
%
% 4. Solve the linear equation system
%
% 5. Re-assemble to the complete vector of unknowns
%
% 6. Compute the complete force vector
%
%% Function main body

%% 0. Read input

% Assign back dummy output variables
loadFactor = 'undefined';
hasConverged = 'undefined';

%% 1. Compute the linear matrices of the steady-state problem
if isa(computeProblemMatricesSteadyState, 'function_handle')
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(tab, '>> Computing the stiffness matrix of the system\n'));
    end
    [stiffMtx, resVct, minElSize] = computeProblemMatricesSteadyState...
        (propAnalysis, u, uSaved, uDot, uDotSaved, precompStiffMtx, ...
        precomResVct, DOFNumbering, mesh, F, loadFactor, computeBodyForces, ...
        propTransientAnalysis, parameters, int);
    if ~ischar(precompStiffMtx)
        stiffMtx = stiffMtx + precompStiffMtx;
    end
    if ~ischar(precomResVct)
        resVct = resVct + precomResVct;
    end
else
    minElSize = 'undefined';
    if ~ischar(precompStiffMtx)
        stiffMtx = precompStiffMtx;
    end
    if ~ischar(precomResVct)
        resVct = precomResVct;
    end
end

%% 2. Compute the stiffness matrix and the residual vector for the transient problem
if isfield(propTransientAnalysis, 'timeDependence')
    if ~strcmp(propTransientAnalysis.timeDependence, 'steadyState')
        if isa(propTransientAnalysis.computeProblemMtrcsTransient, 'function_handle') && ...
            ~propTransientAnalysis.isStaticStep
            [stiffMtx, resVct] = ...
                propTransientAnalysis.computeProblemMtrcsTransient...
                (u, uSaved, uDot, uDotSaved, uDDot, uDDotSaved, massMtx, ...
                dampMtx, stiffMtx, resVct, propTransientAnalysis);
        elseif ~isa(propTransientAnalysis.computeProblemMtrcsTransient, 'function_handle') && ...
                ~(strcmp(propTransientAnalysis.timeDependence, 'steadyState') || strcmp(propTransientAnalysis.timeDependence, 'pseudotransient'))
            error('Variable propTransientAnalysis.computeProblemMtrcsTransient is undefined but the simulation is transient')
        end
    end
end

%% 3. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
if norm(valuesInhomDOFs) ~= 0
    resVct = resVct - stiffMtx(:, inhomDOFs)*valuesInhomDOFs';
end

%% 4. Solve the linear equation system
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab, '>> Solving the linear system of %d equations\n'), length(freeDOFs));
end
[uRed, hasLinearSystemConverged] = ...
    solve_LinearSystem(stiffMtx(freeDOFs, freeDOFs),resVct(freeDOFs), u(freeDOFs));
if ~hasLinearSystemConverged
    error('Linear equation solver has not converged');
end

%% 5. Re-assemble to the complete vector of unknowns
u(freeDOFs) = uRed;
u(homDOFs) = 0;
u(inhomDOFs) = valuesInhomDOFs;

%% 6. Compute the complete flux vector
if strcmp(outMsg, 'outputEnabled')
    fprintf(strcat(tab, '>> Computing the complete force vector\n'));
end
FComplete = stiffMtx*u;

end
