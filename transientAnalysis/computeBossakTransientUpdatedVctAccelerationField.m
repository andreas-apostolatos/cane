function [uDot, uDDot] = computeBossakTransientUpdatedVctAccelerationField ...
    (u, uSaved, uDotSaved, uDDotSaved, propTransientAnalysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the updated values of the discretized field and its time 
% derivatives nessecary for the backward Euler time integration scheme
%
%                 Input :
%                     u : the discrete primary field of the current time
%                         step
%                uSaved : the discrete primary field of the previous time
%                         step
%             uDotSaved : Discrete first time derivative of the primary
%                         field from the previous time step
%            uDDotSaved : Dicrete second time derivative of the primary
%                         field from the previous time step
% propTransientAnalysis : On the transient analysis :
%                          .method : The time integration method
%                          .TStart : Start time of the simulation
%                            .TEnd : End time of the simulation
%                              .nT : Number of time steps
%                              .dt : Time step  
%                Output : 
%                  uDot : updated first time derivative at the current
%                         time step
%                 uDDot : updated second time derivative at the current
%                         time step
%
% Function layout :
%
% 0. Read input
%
% 1. Update the first time derivative of the displacement field
%
% 2. Update the second time derivative of the displacement field
%
%% Function main body

%% 0. Read input
dt = propTransientAnalysis.dt;
betaB = propTransientAnalysis.betaB;
gammaB = propTransientAnalysis.gammaB;

%% 1. Update the first time derivative of the displacement field
uDot = (1 - gammaB/betaB)*uDotSaved - gammaB/betaB/dt*uSaved + ...
    dt*(2*betaB - gammaB)/2/betaB*uDDotSaved + gammaB/betaB/dt*u;

%% 2. Update the second time derivative of the displacement field
uDDot = 1/betaB/dt^2*(u - uSaved) - 1/betaB/dt*uDotSaved - ...
    (1 - 2*betaB)/2/betaB*uDDotSaved;

end
