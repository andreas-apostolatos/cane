function FAmp = computeLoadPseudoWind4GenoaSail(x, y, z, t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the amplitude of a force resembling a wind load onto a genoa sail
% which is not in stall position.
%
%   Input :
%   x,y,z : The physical coordinates of point onto the wind sail
%       t : The time instance
%
%  Output :
%    FAmp : The amplitude of the load at the given point
%
%% Function main body
if t < 20
    Length = 5;
    Height = 9;
    FAmp = - abs(1e3*sin(2*pi()*x/Length/1e0)*sin(2*pi()*z/Height/1e0)*sin(2*pi()*t/10));
else
    FAmp = 0;
end

end
