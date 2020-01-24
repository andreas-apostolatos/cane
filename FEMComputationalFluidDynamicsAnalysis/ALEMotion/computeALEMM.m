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

if ~ischar(t)
   dy = (y/p1)*t;
else
  dy = 0;
end
dx = 0;
dz = 0;

end
