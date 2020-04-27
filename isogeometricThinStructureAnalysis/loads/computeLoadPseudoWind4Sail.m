function FAmp = computeLoadPseudoWind4Sail(x, y, z, t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the amplitude of a force resembling a wind load onto a sail.
%
%   Input :
%   x,y,z : The physical coordinates of point onto the wind sail
%       t : The time instance
%
%  Output :
%    FAmp : The amplitude of the load at the given point
%
%% Function main body

% For the sail
FAmp = 1e3; % 1e3

% End time for the load application
TStarLoad = -6;
TEndLoad = 100;

% % For the four point sail
% FAmp = 1e6;

% % End time for the load application
% TStarLoad = 19;
% TEndLoad = 20;

if t <= TEndLoad && t > TStarLoad
%     Length = 5;
%     Height = 9;
%     FAmp = 1e3*sin(2*pi()*x/Length/1e0)*sin(2*pi()*z/Height/1e0)*sin(2*pi()*t/10);
% FAmp = FAmp*sin(pi()*t/TEndLoad/2);
FAmp = FAmp*(t - TStarLoad)/(TEndLoad - TStarLoad);
% FAmp = FAmp*t/TEndLoad;
else
    FAmp = 0;
    return;
end

end
