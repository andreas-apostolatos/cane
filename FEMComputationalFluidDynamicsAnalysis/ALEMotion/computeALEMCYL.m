function [dx,dy,dz] = ...
    computeALEMCYL(x,y,z,t,propUser)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos, Matthew Keller
%
%% Function documentation
%
% Returns the movement of a node according to a user-defined function for
% the movement of the boundary
%
%                    Input :
%                        x : x coordinate of the current node
%                        y : y coordinate of the current node
%                        z : z coordinate of the current node
%                        t : Current time step value
%                 propUser : Extra user-defined properties
%
%                   Output :
%                       dx : movement in x coordinate of the current node
%                       dy : movement in y coordinate of the current node
%                       dz : movement in z coordinate of the current node

%% Function main body

delta_p1 = propUser.delta_p1;
iterate_p1 = propUser.iterate_p1;
Perturb_Flag = propUser.Perturb_Flag;
x_Mid = propUser.x_Mid;
y_Mid = propUser.y_Mid;
partition_x = abs(x-x_Mid);
partition_y = abs(y-y_Mid);
angle = atand(partition_y/partition_x);

if strcmp(Perturb_Flag, 'perturb')
    perturb_val = -delta_p1;
elseif strcmp(Perturb_Flag, 'set_final')
    perturb_val = iterate_p1;
end

if ~ischar(t)
   if x < x_Mid
       if y < y_Mid
            dx = perturb_val*cosd(angle);
            dy = perturb_val*sind(angle);
       elseif y > y_Mid
            dx = perturb_val*cosd(angle);
            dy = -perturb_val*sind(angle);
       elseif y == y_Mid
            dx = perturb_val;
            dy = 0;
       end
   elseif x > x_Mid
       if y < y_Mid
            dx = -perturb_val*cosd(angle);
            dy = perturb_val*sind(angle);
       elseif y > y_Mid
            dx = -perturb_val*cosd(angle);
            dy = -perturb_val*sind(angle);
       elseif y == y_Mid
            dx = -perturb_val;
            dy = 0;
       end
   elseif x == x_Mid
       if y < y_Mid
            dy = perturb_val;
            dx = 0; 
       elseif y > y_Mid
            dy = -perturb_val;
            dx = 0; 
       elseif y == y_Mid
            dy = 0;
            dx = 0; 
       end
   end
else
  dy = 0;
  dx = 0;
end
dz = 0;
end
