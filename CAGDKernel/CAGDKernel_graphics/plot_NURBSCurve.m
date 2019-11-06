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
function plot_NURBSCurve(p,q,U,V,CP,grid,par,uv,uv_start,uv_end)
%% Function documentation
%
% Plots the parametric curve on the NURBS surface
%
%    Input :
%      p,q : Polynomial degrees
%      U,V : Knot vectors in u,v-direction
%       CP : Set of control points and weights
%     grid : How many points to use
%      par : 1-fix parameter v on surface
%            2-fix parameter u on surface
%       uv : coordinate of the parametric line
% uv_start : Start coordinate of the curve
%   uv_end : End coordinate of the curve
%
%   Output : graphics
%           
%% Function main body

[Xp,Yp,Zp] = create_curve(p,q,U,V,CP,grid,par,uv,uv_start,uv_end);

% geometry
%surf(Xp,Yp,'FaceColor','green','EdgeColor','none');
line(Xp,Yp,Zp);
hold on;
 
% element edges
%create_el_edges(p,U,CP)

% control points and polygon
%create_conpolygon_curve(CP) 

axis equal;
camlight left; lighting phong;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
 
hold off;