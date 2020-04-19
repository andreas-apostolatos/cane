function [tanMtx,resVct] = computeProblemMtrcsImplicitEuler...
    (u,uSaved,uDot,uDotSaved,uDDot,uDDotSaved,massMtx,damMtx,tanMtx,...
    resVct,propTransientAnalysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
%% Function documentation
%
% Returns the system matrix and right hand side vector corresponding to the
% explicit Euler time integration scheme.
%
%                 input :
%                     u : Solution of the primary field from the previous 
%                         iteration step
%                uSaved : Solution of the primary field from the previous 
%                         time step
%                  uDot : Solution of the rate of the primary field from 
%                         the previous iteration step
%             uDotSaved : Solution of the rate of the primary field from 
%                         the previous time step
%            uDDotSaved : Solution of the second time derivative of the 
%                         primary field from the previous time step
%               massMtx : The mass matrix of the problem
%                damMtx : System damping matrix
%                tanMtx : System matrix corresponding to the steady-state 
%                         problem
%             resVctRHS : Right-hand side (RHS)/residual vector 
%                         corresponding to the steady-state problem
% propTransientAnalysis : Transient analysis parameters:
%                                   .method : Time integration method
%                                .alphaBeta : Bossak parameter
%                                    .gamma : Bossak parameter
%                                   .TStart : Start time of the simulation
%                                     .TEnd : End time of the simulation
%                                       .nT : Number of time steps
%                                       .dt : Time step
%
%                output :
%                tanMtx : The updated system matrix corresponding to the
%                         implicit Euler time integration scheme
%                resVct : Right-hand side (RHS)/residual vector 
%                         corresponding to the implicit Euler time 
%                         integration scheme
%
% Function main body :
%
% 1. Compute the problem matrix considering the inertia forces for the implicit Euler time integration scheme
%
% 2. Compute the right hand-side (RHS)/residual vector considering the inertia forces for the implicit Euler time integration scheme
%
%% Function main body

%% 0. Read input
dt = propTransientAnalysis.dt;

%% 1. Compute the problem matrix considering the inertia forces for the implicit Euler time integration scheme
if ischar(damMtx)
    tanMtx = tanMtx + ...  % steady-state
        massMtx*(1/dt);    % transient (implicit Euler)
else
    error('Damping not yet supported for this time integration type');
end

%% 2. Compute the right hand-side (RHS)/residual vector considering the inertia forces for the implicit Euler time integration scheme
if ischar(damMtx)
    resVct = massMtx*u/dt;  % transient (implicit Euler)
else
    error('Damping not yet supported for this time integration type');
end

end