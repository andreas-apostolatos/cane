function [tanMtx, resVct] = computeProblemMtrcsBossakIGA4NSE ...
    (u, uSaved, uDot, uDotSaved, uDDot, uDDotSaved, massMtx, damMtx, tanMtx, ...
    resVct, propFldDynamics)
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
% Bossak time integration scheme for the isogeometric discretization of the 
% Navier-Stokes equations.
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
%                damMtx : Dummy array for this function
%                tanMtx : System matrix corresponding to the steady-state 
%                         problem
%             resVctRHS : Right-hand side (RHS)/residual vector 
%                         corresponding to the steady-state problem
%       propFldDynamics : Structure containing information on the time
%                         integration for the transient fluid dynamics
%                         analysis,
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
%% Function main body

% This function needs to be investigated because the mass matrix is
% dependent on the stabilization terms

end