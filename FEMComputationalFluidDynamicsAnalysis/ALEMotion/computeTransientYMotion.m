function [dx, dy, dz] = ...
    computeTransientYMotion(x, y, z, t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the movement of a node according to a user-define function for
% the movement of the boundary
%
%                    Input :
%                        x : x coordinate of the current node
%                        y : y coordinate of the current node
%                        z : z coordinate of the current node
%                        t : Current time step value
%
%                   Output :
%                       dx : movement in x coordinate of the current node
%                       dy : movement in y coordinate of the current node
%                       dz : movement in z coordinate of the current node
%
%% Function main body
amplification = .5;
frequency = 10;
dy = amplification*sin(frequency*pi*t);
dx = 0;
dz = 0;

end
