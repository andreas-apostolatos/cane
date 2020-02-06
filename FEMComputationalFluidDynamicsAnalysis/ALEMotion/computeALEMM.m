function [dx,dy,dz] = ...
    computeALEMM(x,y,z,t,propUser)
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

p1 = propUser.p1;
p2 = propUser.p2;
delta_p1 = propUser.delta_p1;
delta_p2 = propUser.delta_p2;
iterate_p1 = propUser.iterate_p1;
iterate_p2 = propUser.iterate_p2;
Perturb_Flag = propUser.Perturb_Flag;
x_Mid = propUser.x_Mid;

if ~ischar(t)
   if strcmp(Perturb_Flag, 'dx')
       if x < x_Mid
           dx = 0.5*delta_p2;
       elseif x > x_Mid
           dx = -0.5*delta_p2;           
       elseif x == x_Mid
           dx = 0;           
       end         
       dy = 0;       
   elseif strcmp(Perturb_Flag, 'dy')
       dx = 0;
       dy = (y/p1)*delta_p1;
   elseif strcmp(Perturb_Flag, 'dxdy')
       if x < x_Mid
           dx = -0.5*iterate_p2;
       elseif x > x_Mid
           dx = 0.5*iterate_p2;           
       elseif x == x_Mid
           dx = 0;       
       end
       dy = -(y1/p1)*iterate_p1;
   end
else
  dy = 0;
  dx = 0;
end
dz = 0;

end
