function [u,CPHistory,resHistory,hasConverged,FComplete,rankD,condK,...
    minEig,propCoupling,minElAreaSize] = ...
    solve_IGALinearSystem...
    (analysis,uSaved,uDotSaved,uDDotSaved,BSplinePatches,connections,u,uDot,...
    uDDot,KConstant,massMtx,computeLinearMtrcsSteadyState,computeUpdatedGeometry,...
    freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,masterDOFs,slaveDOFs,...
    solve_LinearSystem,t,propCoupling,propTransientAnalysis,...
    propNLinearAnalysis,plot_IGANLinear,tab,graph,outMsg)
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
%                       analysis : .type : The analysis type
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
%                             .homDOFs : The global numbering of the
%                                        DOFs where homogeneous Dirichlet
%                                        boundary conditions are applied
%                           .inhomDOFs : The global numbering of the
%                                        DOFs where homogeneous Dirichlet
%                                        boundary conditions are applied
%                     .valuesInhomDOFs : Prescribed values to the DOFs 
%                                        where homogeneous Dirichlet
%                                        boundary conditions are applied
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
%                      KConstant : The matrices which stay constant 
%                                  throughout the nonlinear computations
%                        massMtx : The mass matrix of the problem
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
%                                         .TEnd : End time of the simulation
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
%                plot_IGANLinear : Dummy variable for this function
%                            tab : Tabulation for the output messages
%                          graph : Dummy variable for this function
%                         outMsg : On printing information during analysis 
%                                  on the command window
%
%                         Output :
%                              u : The converged discrete solution vector
%                      CPHistory : Dummy output for this function
%                     resHistory : Dummy output for this function
%                   hasConverged : Dummy output for this function
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
% 1. Compute the linear matrices of the system
%
% 2. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
%
% 3. Solve the linear equation system
%
% 4. Re-assemble to the complete vector of unknowns
%
% 5. Compute the complete force vector
%
% 6. Appendix
%
%% Function main body

%% 0. Read input

% Assign back dummy output variables
CPHistory = 'undefined';
loadFactor = 'undefined';
resHistory = 'undefined';
hasConverged = true;

% Damping is still not implemented
damMtx = 'undefined';

%% 1. Compute the linear matrices of the system
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'>> Computing the system matrix and right hand side vector\n'));
end
[stiffMtx,RHS,minElAreaSize] = computeLinearMtrcsSteadyState(KConstant,u,...
            uSaved,uDot,uDotSaved,BSplinePatches,connections,...
            propCoupling,propTransientAnalysis,t,strcat(tab,'\t'),...
            loadFactor,outMsg);
if isa(propTransientAnalysis.computeProblemMtrcsTransient,'function_handle') && ...
        ~propTransientAnalysis.isStaticStep
    [stiffMtx,RHS] = propTransientAnalysis.computeProblemMtrcsTransient...
        (u,uSaved,uDot,uDotSaved,uDDot,uDDotSaved,massMtx,damMtx,...
        stiffMtx,RHS,propTransientAnalysis);
end
if strcmp(outMsg,'outputEnabled')
    rankD = length(RHS(freeDOFs)) - rank(stiffMtx(freeDOFs,freeDOFs));
    if rankD ~= 0
        fprintf(strcat(tab,'>> The system has rank deficiency equal to %d\n'),rankD);
    end
    condK = cond(stiffMtx(freeDOFs,freeDOFs));
    fprintf(strcat(tab,'>> The condition number of the system is %d\n'),condK); 
    [~,eigVal] = eig(stiffMtx(freeDOFs,freeDOFs));
    eigVal(~eigVal) = Inf;
    minEig = min(min(eigVal));
    fprintf(strcat(tab,'>> The minimum eigenvalue of the system is %d\n'),minEig); 
else
    rankD = 'undefined';
    condK = 'undefined';
    minEig = 'undefined';
end

%% 2. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
if norm(valuesInhomDOFs(inhomDOFs)) ~= 0
    RHS = RHS - stiffMtx(:,inhomDOFs)*valuesInhomDOFs(inhomDOFs);
end

%% 3. Solve the linear equation system
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'>> Solving the linear system of %d equations\n'),length(RHS(freeDOFs)));
end
[uRed,hasLinearSystemConverged] = solve_LinearSystem...
    (stiffMtx(freeDOFs,freeDOFs),RHS(freeDOFs),u(freeDOFs));
if ~hasLinearSystemConverged
    error('Linear equation solver has not converged');
end

%% 4. Re-assemble to the complete vector of unknowns
u(freeDOFs) = uRed;
u(homDOFs) = 0;
u(inhomDOFs) = valuesInhomDOFs(inhomDOFs);

%% 5. Compute the complete force vector
FComplete = stiffMtx*u;

%% 6. Appendix
if strcmp(outMsg,'outputEnabled')
    fprintf('\n');
end

end
