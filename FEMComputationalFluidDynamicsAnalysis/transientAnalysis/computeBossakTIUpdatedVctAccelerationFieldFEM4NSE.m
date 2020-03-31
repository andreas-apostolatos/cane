function [upDot, upDDot] = ...
    computeBossakTIUpdatedVctAccelerationFieldFEM4NSE ...
    (up, upSaved, upDotSaved, upDDotSaved, propFldDynamics)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the updated acceleration field corresponding to the application
% of the Bossak scheme in the classical finite element discretization of 
% the incompressible Navier-Stokes equations in 2D. The function has also a 
% dummy argument upDDot which would be used in a possible computational 
% structural dynamics problem.
%
%             Input :
%                up : The converged discrete solution of the current time 
%                     step
%           upSaved : The converged discrete solution of the previous time 
%                     step
%        upDotSaved : The updated discrete time derivative of the solution 
%                     from the previous iteration step
%       upDDotSaved : The updated discrete second order time derivative of 
%                     the solution from the previous iteration step (dummny 
%                     variable for this function)
%   propFldDynamics : Structure containing information on the time
%                     integration of the fluid dynamics analysis,
%                               .method : Time integration method
%                            .alphaBeta : Bossak parameter
%                                .gamma : Bossak parameter
%                               .TStart : Start time of the simulation
%                                 .TEnd : End time of the simulation
%                                   .nT : Number of time steps
%                                   .dt : Time step (numeric or adaptive)
%
%           Output :
%            upDot : Discrete acceleration field
%           upDDot : Second order time derivative of the velocity field 
%                    (dummny output variable for this function)
%
%% Function main body

% Discrete acceleration field
upDot = 1/propFldDynamics.gamma/propFldDynamics.dt*(up - upSaved) - ...
            (1 - propFldDynamics.gamma)/propFldDynamics.gamma*upDotSaved;
        
% Second order time derivative of the velocity field
upDDot = 'undefined';

end
