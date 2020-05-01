function plot_ControlPolygonBSplineSurface(CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots control points and control polygon for the given set of Control
% Points
%
%   Input :
%      CP : Set of Control Point coordinates
%
%  Output : Graphics
%
%% Function main body
for iCP = 1:length(CP(1, :, 1)) - 1
    for jCP = 1:length(CP(:, 1, 1)) - 1
        plot3(CP(jCP, iCP:iCP + 1, 1), CP(jCP, iCP:iCP + 1, 2), ...
            CP(jCP, iCP:iCP + 1, 3), '--or');
        plot3(CP(jCP:jCP + 1, iCP, 1), CP(jCP:jCP + 1, iCP, 2), ...
            CP(jCP:jCP + 1, iCP, 3), '--or');
    end
    
    % Update counter
    jCP = jCP + 1;
    plot3(CP(jCP, iCP:iCP + 1, 1), CP(jCP, iCP:iCP + 1, 2), ...
        CP(jCP, iCP:iCP + 1, 3), '--or');
end
iCP = iCP + 1;
for jCP =1:length(CP(:, 1, 1)) - 1
    plot3(CP(jCP:jCP + 1, iCP, 1), CP(jCP:jCP + 1, iCP, 2), ...
        CP(jCP:jCP + 1, iCP, 3), '--or');
end

end
