function [uDot,uDDot] = computeBETIUpdatedVctAccelerationFieldFEMPlateInMembraneAction...
    (u,uSaved,uDotSaved,uDDotSaved,propTransientAnalysis)
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
%                     u : the discrete primary field of timestep n
%                uSaved : the discrete primary field saved
%             uDotSaved : first derivative saved
%            uDDotSaved : second derivative saved
% propTransientAnalysis : On the transient analysis :
%                          .method : The time integration method
%                          .TStart : Start time of the simulation
%                            .TEnd : End time of the simulation
%                              .nT : Number of time steps
%                              .dt : Time step  
%                Output : 
%                  uDot : updated first time derivative at time step n + 1
%                 uDDot : updated second time derivative at time step n + 1
%
%% Function main body

% Update the first time derivative of the displacement field
uDot = (u - uSaved)/propTransientAnalysis.dt;

% Update the second time derivative of the displacement field
uDDot = (u - uSaved -uDotSaved*propTransientAnalysis.dt)/(propTransientAnalysis.dt^2);

end
