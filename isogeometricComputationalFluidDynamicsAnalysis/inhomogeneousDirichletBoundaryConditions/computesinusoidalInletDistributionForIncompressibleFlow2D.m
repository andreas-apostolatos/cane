%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function velocity = sinusoidalInletDistributionForIncompressibleFlow2D(x,y,z,t)
%% Function documentation
%
% Returns the velocity at the given location y on a vertical (y-directed) 
% line which is defined by a lower and an upper value w.r.t. its y-coordinate.
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
Height = 4.0;

% Y-component of the lower edge
ylo = -Height/2;

% Y-component of the upper edge
yup = Height/2;

% Start time of the simulation
Tstart = 10;

% end time of the simulation
Tend = 0;

% Spatial frequency of the prescribed motion
omega = 4;

% Compute the scalar value at the given location
if t<=5.13*10^-1
    if y>=ylo && y<=yup
        velocity = Smax*abs(sin(omega*pi*y/(yup-ylo)))*cos(5*pi*t/(Tend-Tstart));
    else
        velocity = 0.0;
    end
else
    if y>=ylo && y<=yup
        velocity = Smax*abs(sin(omega*pi*y/(yup-ylo)))*cos(5*pi*t/(Tend-Tstart));
    else
        velocity = 0.0;
    end
end

end

