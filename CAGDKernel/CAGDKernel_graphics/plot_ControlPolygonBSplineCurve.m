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
function plot_ControlPolygonBSplineCurve(CP)
%% Function documenation
% 
% Plots Control Points and Control Polygon corresponding to a B-Spline 
% curve
%
%   Input :
%      CP : The control points for a B-Spline curve
%
%  Output :
%           graphics
%
%% Function main body

% Plot the intermediate Control Points and their connections
for l=1:length(CP(:,1))-1
    plot3(CP(l,1),CP(l,2),CP(l,3),'--or');
    plot3(CP(l:l+1,1),CP(l:l+1,2),CP(l:l+1,3),'--or');
end

end