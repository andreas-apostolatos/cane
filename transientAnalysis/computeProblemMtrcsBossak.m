function [tanMtx, resVct] = computeProblemMtrcsBossak ...
    (u, uSaved, uDot, uDotSaved, uDDot, uDDotSaved, massMtx, damMtx, ...
    tanMtx, resVct, propTransientAnalysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the system matrix and right hand side vector corresponding to the
% Bossak time integration scheme.
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
%                 uDDot : Solution of the second time derivative of the 
%                         primary field
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
%                                   .alphaB : Bossak parameter
%                                    .betaB : Bossak parameter
%                                   .gammaB : Bossak parameter
%                                   .TStart : Start time of the simulation
%                                     .TEnd : End time of the simulation
%                              .noTimeSteps : Number of time steps
%                                       .dt : Time step (numeric or 
%                                             adaptive)
%
%                output :
%                tanMtx : The updated system matrix corresponding to the
%                         Bossak time integration scheme
%                resVct : Right-hand side (RHS)/residual vector 
%                         corresponding to the Bossak time integration 
%                         scheme
%
% Function main body :
%
% 0. Read input
%
% 1. Compute the problem matrix considering the inertia and the damping forces for the Bossak time integration scheme
%
% 2. Compute the right hand-side (RHS)/residual vector considering the inertia and the damping forces for the Bossak time integration scheme
%
%% Function main body

%% 0. Read input
dt = propTransientAnalysis.dt;
alphaB = propTransientAnalysis.alphaB;
betaB = propTransientAnalysis.betaB;
gammaB = propTransientAnalysis.gammaB;

%% 1. Compute the problem matrix considering the inertia and the damping forces for the Bossak time integration scheme
if ~ischar(damMtx)
    tanMtx = tanMtx + ...  % steady-state
        (1 - alphaB)/betaB/dt^2*massMtx + gammaB/betaB/dt*damMtx; % transient (Bossak)
else
    tanMtx = tanMtx + ...  % steady-state
        (1 - alphaB)/betaB/dt^2*massMtx; % transient (Bossak)
end

%% 2. Compute the right hand-side (RHS)/residual vector considering the inertia and the damping forces for the Bossak time integration scheme
if ~ischar(damMtx)
    resVct = resVct + ...  % steady-state
        ((1 - alphaB)/betaB/dt^2*massMtx + gammaB/betaB/dt*damMtx)*u - ...
        ((1 - alphaB)/betaB/dt^2*massMtx + gammaB/betaB/dt*damMtx)*uSaved - ...
        ((1 - alphaB)/betaB/dt*massMtx - (betaB - gammaB)/betaB*damMtx)*uDotSaved - ...
        ((1 - 2*betaB - alphaB)/2/betaB*massMtx - dt*(2*betaB - gammaB)/2/betaB*damMtx)*uDDotSaved;  % transient (BETI)
else
    resVct = resVct + ...  % steady-state
        ((1 - alphaB)/betaB/dt^2*massMtx)*u - ...
        ((1 - alphaB)/betaB/dt^2*massMtx)*uSaved - ...
        ((1 - alphaB)/betaB/dt*massMtx)*uDotSaved - ...
        ((1 - 2*betaB - alphaB)/2/betaB*massMtx)*uDDotSaved;  % transient (BETI)
end

end
