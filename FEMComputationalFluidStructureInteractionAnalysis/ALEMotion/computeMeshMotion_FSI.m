function [dx, dy, dz] = ...
    computeMeshMotion_FSI(x, y, z, t, propUser)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the user-defined mesh motion corresponding to 
% Fluid-Structure-Interaction 
%
%       Input :
%    dx,dy,dz : The resulting mesh motion in terms of displacement
%     x,y,z,t : Chorochronical location
%    propUser : User-defined displacement,
%                 .u,.v,.w : Displacement along the X-, Y- and Z-coorindate
%
%% Function main body

dx = propUser.u;
dy = propUser.v;
dz = propUser.w;

end