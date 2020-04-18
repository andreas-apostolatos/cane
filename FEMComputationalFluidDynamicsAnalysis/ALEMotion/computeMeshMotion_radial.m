function [dx, dy, dz] = ...
    computeMeshMotion_radial(x, y, z, t, propUser)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Matthew Keller
%
%% Function documentation
%
% Returns the radial displacement of a point which is assumed to be on a
% circle with center (propUser.x_Mid, propUser.y_Mid) in 2D.
%
%     Input :
%     x,y,z : x-,y- and z-coordinates of the point
%         t : Time instance
%  propUser : User-defined properties,
%                .x0, .y0 : The center of the circle
%                     .dr : Radial displacement
%
%    Output :
%  dx,dy,dz : Displacement in x-,y- and z-direction of the point
%
%% Function main body

%% 0. Check input

% Get the center of the circle
if ~isfield(propUser, 'x0')
    error('propUser must define x0');
end
x0 = propUser.x0;
if ~isfield(propUser, 'y0')
    error('propUser must define y0');
end
y0 = propUser.y0;

% Get the radial displacement
if ~isfield(propUser, 'dr')
    error('propUser must define dr');
end
dr = propUser.dr;

%% 1. Compute the angle by which the point is turned on the circle
rcostheta = abs(x - x0);
rsintheta = abs(y - y0);
angle = atand(rsintheta/rcostheta);

%% 2. Compute the displacement along x- and y-direction given the incremental radial displacement
if x < x0
   if y < y0
        dx = -dr*cosd(angle);
        dy = -dr*sind(angle);
   elseif y > y0
        dx = -dr*cosd(angle);
        dy = dr*sind(angle);
   elseif y == y0
        dx = -dr;
        dy = 0;
   end
elseif x > x0
   if y < y0
        dx = dr*cosd(angle);
        dy = -dr*sind(angle);
   elseif y > y0
        dx = dr*cosd(angle);
        dy = dr*sind(angle);
   elseif y == y0
        dx = dr;
        dy = 0;
   end
elseif x == x0
   if y < y0
        dy = -dr;
        dx = 0; 
   elseif y > y0
        dy = dr;
        dx = 0; 
   elseif y == y0
        dy = 0;
        dx = 0; 
   end
end
dz = 0;

end