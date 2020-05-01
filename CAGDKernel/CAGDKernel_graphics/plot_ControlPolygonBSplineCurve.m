function plot_ControlPolygonBSplineCurve(CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
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
for l=1:length(CP(:,1))-1
    plot3(CP(l,1),CP(l,2),CP(l,3),'--or');
    plot3(CP(l:l+1,1),CP(l:l+1,2),CP(l:l+1,3),'--or');
end

end
