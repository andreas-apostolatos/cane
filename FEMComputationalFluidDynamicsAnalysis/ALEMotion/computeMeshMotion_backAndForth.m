function [dx, dy, dz] = ...
    computeMeshMotion_backAndForth(x, y, z, t, varargin)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the updated nodal coordinates for a point which is moving
% horizontally in a sinusoidal manner. This function is used for the unit
% test case created for unit testing the mesh motion algorithm.
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
amplification = .3;
frequency = 1/2;
if ~ischar(t)
    dx = amplification*sin(frequency*pi*t);
else
    dx = 0;
end
dy = 0;
dz = 0;

end
