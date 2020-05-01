function index = plot_BSplineCurveOnBSplineSurface ...
    (p_curve, Xi_curve, CP_curve, isNURBS_curve, ...
    p_surface, q_surface, Xi_surface, Eta_surface, CP_surface, ...
    isNURBS_surface, numEval, color_curve, color_surface, lineWidth, graph)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots a B-Spline curve which is embedded in the parametric space of a
% B-Spline surface.
%
%           Input :
%         p_curve : The polynomial order of the curve
%        Xi_curve : The knot vector of the curve in its parametric space
%        CP_curve : The set of Control Point coordinates and weights which
%                   are living in the parametric space of the surface 
%                   Xi_surface times Eta_surface
%   isNURBS_curve : Flag on whether the curve basis is a B-Spline or a 
%                   NURBS
%    p*,q_surface : The polynomial degrees of the B-Spline surface
% Xi*,Eta_surface : The knot vectors of the B-Spline surface
%      CP_surface : The Control Point coordinates and the weights of the
%                   B-Spline surface which are living in the physical space
% isNURBS_surface : Flag on whether the surface basis is a B-Spline or a 
%                   NURBS
%            grid : Number of grid points for the creation of the surface
%     color_curve : Color of the embedded curve
%   color_surface : Color of the surface
%       lineWidth : Width of the curve
%           graph : Structure containing,
%                       .index : Index of the current graph
%
% Function layout :
%
%% Function main body

%% 1. Plot the B-Spline surface in the physical space
figure(graph.index)
plot_BSplineSurfacePatch(p_surface,q_surface,Xi_surface,Eta_surface,...
    CP_surface,isNURBS_surface,color_surface);
hold on;

%% 2. Plot the embedded B-Spline curve in the surface's parametric space
[Xp,Yp,Zp] = createBSplineSurfaceEmbeddedBSplineCurveOnCartesianSpace....
    (p_surface,q_surface,Xi_surface,Eta_surface,CP_surface,isNURBS_surface,...
    p_curve,Xi_curve,CP_curve,isNURBS_curve,numEval);
line(Xp,Yp,Zp,'Linewidth',lineWidth,'color',color_curve);
% plot_ControlPolygonBSplineCurveOnBSplineSurface...
%     (CP_curve,p_surface,q_surface,Xi_surface,Eta_surface,CP_surface,isNURBS_surface);
axis equal;
hold off;

%% 3. Update graphs index
index = graph.index + 1;

end
