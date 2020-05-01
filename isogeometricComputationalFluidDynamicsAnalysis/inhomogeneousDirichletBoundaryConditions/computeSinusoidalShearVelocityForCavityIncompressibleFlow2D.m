function velocity = ...
    computeSinusoidalShearVelocityForCavityIncompressibleFlow2D ...
    (x, y, z, t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the velocity at the given location y on a vertical (x-directed) 
% line which is defined by a lower and an upper value w.r.t. its x-coordinate.
% Especially suited for the cavity benchmark.
%
%    Input :
%    x,y,z : Cartesian coordinates of the node where the prescribed 
%            (inhomogeneous) Dirichlet boundary conditions are applied
%        t : The time instance where to compute the prescribed Dirichlet 
%            boundary conditions
%
%   Output :
% velocity : The concetration of the scalar quantity at the given location    
%
%% Function main body

% Maximum concetration of the scalar in the center of the application line
Smax = 6.0;

% Height of the channel
Length = 1.0;

% X-component of the left edge
xlo = 0;

% X-component of the right edge
yup = Length;

% Start time of the simulation
Tstart = 0;

% end time of the simulation
Tend = 10;

% Spatial frequency of the prescribed motion
omega = 1e1;

% Compute the scalar value at the given location
velocity = 4*Smax*Length*x*(x-Length)*cos(omega*t);

end
