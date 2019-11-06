function [uDot,uDDot] = computeBETITransientUpdatedVctAccelerationField...
    (u,uSaved,uDotSaved,uDDotSaved,propTransientAnalysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentationx
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
%                         field from the previous time step (dummy variable
%                         for this function)
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
% 1.Update the first time derivative of the displacement field
%
% 2. Update the second time derivative of the displacement field
%
%% Function main body

%% 0. Read input
dt = propTransientAnalysis.dt;

%% 1.Update the first time derivative of the displacement field
uDot = (u - uSaved)/dt;

%% 2. Update the second time derivative of the displacement field
uDDot = (u - uSaved)/dt^2 - uDotSaved/dt;

end
