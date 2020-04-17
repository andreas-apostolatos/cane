function [dx,dy,dz] = ...
    computeALEMF(x, y, z, t, varargin)
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
% amplification = 1e-4;
% if (-4.0<=x) && (x<=4.0)
%     dy = 0.5*0.30*(1-cos(2*pi*t/26.3158/amplification));
% elseif (4.0<x) && (x<6.5)
%     dy = 0.5*(0.5*0.30*(1-cos(2*pi*t/26.3158/amplification)))*(1-tanh(10*(x-5.25)));
% elseif (-6.5<x) && (x<-4.0)
%     dy = 0.5*(0.5*0.30*(1-cos(2*pi*t/26.3158/amplification)))*(1-tanh(10*(-x-5.25)));
% else
%     dy = 0;
% end

amplification = .5;
frequency = 10;
if ~ischar(t)
  dy = amplification*sin(frequency*pi*t);
else
  dy = 0;
end
dx = 0;
dz = 0;

end
