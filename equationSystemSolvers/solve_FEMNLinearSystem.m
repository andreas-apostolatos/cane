function [u, FComplete, isConverged,minElASize] = ...
    solve_FEMNLinearSystem ...
    (propAnalysis, uSaved, uDotSaved, uDDotSaved, mesh, F, computeBodyForces, ...
    propParameters, u, uDot, uDDot, massMtx, dampMtx, precompStiffMtx, ...
    precomResVct, computeNLinearMtrcsSteadyState, DOFNumbering, ...
    freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, uMeshALE, ...
    solve_LinearSystem, propTransientAnalysis , t, propNLinearAnalysis, ...
    propGaussInt, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the solution to a nonlinear system which correspond to the
% classical finite element discretization of the underlying field.
%
%                          Input :
%                   propAnalysis : Structure containing information on the
%                                  analysis,
%                                       analysis.type : Analysis type
%                         uSaved : The discrete solution field of the 
%                                  previous time step
%                      uDotSaved : The time derivative of the discrete 
%                                  solution 
%                                  field of the previous time step
%                     uDDotSaved : The second order time derivative of the 
%                                  discrete solution field of the previous 
%                                  time step (dummy variable for this 
%                                  function)
%                           mesh : The nodes and the elements of the 
%                                  underlying mesh
%                              F : The boundary force vector
%              computeBodyForces : Function handle to the computation of 
%                                  the body force vector
%                 propParameters : The parameters of the physical field
%                              u : Initial guess for the primary field
%                           uDot : Initial guess for the time derivative of 
%                                  the primary field
%                          uDDot : Initial guess for the second time 
%                                  derivative of the primary field
%                        massMtx : The mass matrix of the system
%                        dampMtx : The damping matrix of the system
%                precompStiffMtx : Precomputed part of the stiffness matrix
%                   precomResVct : Precomputed part of the residual vector
% computeNLinearMtrcsSteadyState : Function handle for the computation of 
%                                  the matrices and vectors nessecary for 
%                                  the update of the discrete solution 
%                                  field corresponding to steady-state
%                   DOFNumbering : The global numbering of the DOFs
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
%                       uMeshALE : Velocity of the mesh in case of an ALE
%                                  formulation
%             solve_LinearSystem : Function handle to the solver for the 
%                                  linear equation system
%          propTransientAnalysis : Transient analysis parameters:
%                                       .method : Time integration method
%                                    .alphaBeta : Bossak parameter
%                                        .gamma : Bossak parameter
%                                           .T0 : Start time of the 
%                                                 simulation
%                                         .TEnd : End time of the 
%                                                 simulation
%                                           .nT : Number of time steps
%                                           .dt : Time step (numeric or 
%                                                 adaptive)
%                              t : The current time of the transient 
%                                  simulation
%            propNLinearAnalysis : Nonlinear analysis parameters
%                                       .method : The nonlinear solution 
%                                                 scheme
%                                          .eps : The residual tolerance
%                                      .maxIter : The maximum number of 
%                                                 nonlinear iterations
%                   propGaussInt : On the spatial integration
%                                           .type : 'default', 'user'
%                                     .domainNoGP : Number of Gauss Points 
%                                                   for the domain 
%                                                   integration
%                                   .boundaryNoGP : Number of Gauss Points 
%                                                   for the boundary 
%                                                   integration
%                            tab : Tabulation for the output in the command
%                                  window
%                         outMsg : On printing information during analysis 
%                                  in the command window
%
%                         Output :
%                              u : The converged discrete solution vector
%                      FComplete : The complete flux vector of the system
%                    isConverged : Flag on whether the nonlinear iterations 
%                                  have converged
%                     minElASize : The minimum element edge size in the mesh
%
% Function layout :
%
% 1. Loop over all load steps
% ->
%    1i. Compute the load factor
%
%   1ii. Loop over all the Newton-Rapson iterations
%   ->
%        1ii.1. Compute the tangent stiffness matrix and the residual vector for the steady-state problem
%
%        1ii.2. Compute the tangent stiffness matrix and the residual vector for the transient problem
%
%        1ii.3. Compute the right-hand side (RHS) residual vector in equation
%
%        1ii.4. Check condition for convergence on the residual vector
%
%        1ii.5. Solve the linearized equation system
%
%        1ii.6. Re-assemble to the complete vector of unknowns
%   <-
%
%  1iii. If the nonlinear iterations have not converged break the loop
% <-
%
%% Function main body

% Flag on whether there exist more than one load steps in the simulation
isLoadSteps = true;
tabLoadSteps = '\t';
if propNLinearAnalysis.noLoadSteps == 1
    isLoadSteps = false;
    tabLoadSteps = '';
end

%% 1. Loop over all load steps
for iLS = 1:propNLinearAnalysis.noLoadSteps
    if strcmp(outMsg,'outputEnabled') && isLoadSteps
        fprintf(strcat(tab, 'Load step %d/%d \n'), iLS, propNLinearAnalysis.noLoadSteps);
        fprintf(strcat(tab, '----------------\n'));
        fprintf('\n');
    end
    %% 1i. Compute the load factor
    loadFactor = iLS/propNLinearAnalysis.noLoadSteps;
    
    %% 1ii. Loop over all the Newton-Rapson iterations
    if strcmp(outMsg, 'outputEnabled')
        msgPNR = sprintf(strcat(tab, tabLoadSteps, 'Looping over all the Newton iterations\n', ...
            tab, tabLoadSteps, '-------------------------------------- \n \n'));
        fprintf(msgPNR);
    end
    for iNR = 1:propNLinearAnalysis.maxIter
        %% 1ii.1. Compute the tangent stiffness matrix and the residual vector for the steady-state problem
        [tanMtx, resVct, minElASize] = computeNLinearMtrcsSteadyState ...
            (propAnalysis, u, uSaved, uDot, uDotSaved, uMeshALE, ...
            precompStiffMtx, precomResVct, DOFNumbering, mesh, F, ...
            loadFactor, computeBodyForces, propTransientAnalysis, t, ...
            propParameters, propGaussInt);
        
        %% 1ii.2. Compute the tangent stiffness matrix and the residual vector for the transient problem
        if isfield(propTransientAnalysis, 'timeDependence')
            if ~strcmp(propTransientAnalysis.timeDependence, 'STEADY_STATE')
                if isa(propTransientAnalysis.computeProblemMtrcsTransient, 'function_handle')
                    [tanMtx, resVct] = ...
                        propTransientAnalysis.computeProblemMtrcsTransient...
                        (u, uSaved, uDot, uDotSaved, uDDot, uDDotSaved, ...
                        massMtx, dampMtx, tanMtx, resVct, ...
                        propTransientAnalysis);
                elseif ~isa(propTransientAnalysis.computeProblemMtrcsTransient, 'function_handle') && ...
                        ~(strcmp(propTransientAnalysis.timeDependence, 'steadyState') || strcmp(propTransientAnalysis.timeDependence, 'pseudotransient'))
                    error('Variable propTransientAnalysis.computeProblemMtrcsTransient is undefined but the simulation is transient')
                end
            end
        end

        %% 1ii.3. Compute the right-hand side (RHS) residual vector in equation
        RHS = - resVct;
        if norm(valuesInhomDOFs) ~= 0 && iNR == 1
            RHS = RHS + tanMtx(:, inhomDOFs)*loadFactor*valuesInhomDOFs';
        end

        %% 1ii.4. Check condition for convergence on the residual vector

        % Compute the norm of the residual vector over the free DOFs
        residualNorm = norm(resVct(freeDOFs));

        % Issue a message on the evolution of the residual vector
        if strcmp(outMsg, 'outputEnabled')
            msgNR = sprintf(strcat(tab, tabLoadSteps, '\t||resVct|| = %d at nonlinear iteration No. = %d \n'), residualNorm,iNR);
            fprintf(msgNR);
        end

        % Check the convergence of the Newton iterations
        if residualNorm <= propNLinearAnalysis.eps && iNR ~= 1
            if strcmp(outMsg, 'outputEnabled')
                msgANR = sprintf(strcat('\n', tab, tabLoadSteps, '\tNonlinear iterations converged! \n \n \n'));
                fprintf(msgANR);
            end
            isConverged = true;
            break;
        end

        % If the Newton iterations do not converge after specified limit:
        if iNR == propNLinearAnalysis.maxIter
            if strcmp(outMsg, 'outputEnabled')
                if ~ischar(propTransientAnalysis)
                    warning('Nonlinear iterations did not converge for the fixed time step of %d', propTransientAnalysis.dt);
                else
                    warning('Nonlinear iterations did not converge');
                end
            end

            % Flag on the convergence of the Newton iterations
            isConverged = false;
        end

        %% 1ii.5. Solve the linearized equation system
        [deltauRed, isLinearSysConverged] = solve_LinearSystem ...
            (tanMtx(freeDOFs, freeDOFs), RHS(freeDOFs), u(freeDOFs));
        if ~isLinearSysConverged
            error('Linear equation system solver has not converged');
        end

        %% 1ii.6. Re-assemble to the complete vector of unknowns
        u(freeDOFs) = u(freeDOFs) + deltauRed;
        if iNR == 1
            u(homDOFs) = 0;
            u(inhomDOFs) = valuesInhomDOFs;
        end
    end
    
    %% 1iii. If the nonlinear iterations have not converged break the loop
    if ~isConverged
        break;
    end
end

%% 2. Compute the complete flux vector after convergence
FComplete = tanMtx*u;

end
