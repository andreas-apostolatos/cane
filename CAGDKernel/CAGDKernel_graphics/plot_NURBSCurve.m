function plot_NURBSCurve(p, q, Xi, Eta, CP, grid, par, xiEta, ...
    xiEta_start, xiEta_end)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the parametric curve on the NURBS surface
%
%    Input :
%      p,q : Polynomial degrees
%      U,V : Knot vectors in u,v-direction
%       CP : Set of control points and weights
%     grid : How many points to use
%      par : Integer that takes values 1 or 2,
%               1-fix parameter xi on surface
%               2-fix parameter eta on surface
%       xiEta : coordinate of the parametric line
% xiEta_start : Start coordinate of the curve
%   xiEta_end : End coordinate of the curve
%
%   Output : graphics
%
% Function layout :
%
%
% 0. Read input
%
% 1. Create the points on the curve
%
% 2. Plot the line
%
% 3. Create the element edges
%
% 4. Plot the control polygon
%
% 5. Define graph properties
%
%% Function main body

%% 0. Read input
if par ~= 1 && par ~= 2
    error('Variable par can be either 1 or 2');
end

%% 1. Create the points on the curve
[Xp, Yp, Zp] = create_curve ...
    (p, q, Xi, Eta, CP, grid, par, xiEta, xiEta_start, xiEta_end);

%% 2. Plot the line
% surf(Xp,Yp,'FaceColor','green','EdgeColor','none');
line(Xp,Yp,Zp);
hold on;
 
%% 3. Create the element edges
% create_el_edges(p,U,CP)

%% 4. Plot the control polygon
% create_conpolygon_curve(CP) 

%% 5. Define graph properties
axis equal;
camlight left;
lighting phong;
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
zlabel('z', 'FontSize', 14);
hold off;

end
