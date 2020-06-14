function index = plot_BSplineCurveOnParametricSpaceOfBSplineSurface...
    (p_curve, Xi_curve, CP_curve, isNURBS_curve, Xi_surface, Eta_surface, ...
    numEval, color_curve, color_surface, lineWidth, graph)
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
% isNURBS_surface : Flag on whether the surface basis is a B-Spline or a 
%                   NURBS
%         numEval : Number of grid points for the creation of the surface
%     color_curve : Color of the embedded curve
%   color_surface : Color of the surface
%       lineWidth : Width of the curve
%           graph : Structure containing,
%                       .index : Index of the current graph
%
% Function layout :
%
% 0. Read input
%
% 1. Plot the surface patch in its parametric space
%
% 2. Plot the B-Spline surface in its parametric space
%
% 3. Plot the embedded B-Spline curve in the surface's parametric space
%
% 4. Update graphs index
%
%% Function main body

%% 0. Read input
xiGrid = 10;
etaGrid = 10;

%% 1. Plot the surface patch in its parametric space
surf([Xi_surface(1) Xi_surface(end); Xi_surface(1) Xi_surface(end)],...
    [Eta_surface(1) Eta_surface(1); Eta_surface(end) Eta_surface(end)],...
    [0 0;0 0],'FaceColor',color_surface,'EdgeColor','none');
view(2)
hold on;

%% 2. Plot the B-Spline surface in its parametric space
figure(graph.index)
plot_knotsForBSplineSurfaceOnPArametricSpace...
    (Xi_surface,Eta_surface,xiGrid,etaGrid);

%% 3. Plot the embedded B-Spline curve in the surface's parametric space
[Xp,Yp,Zp] = createBSplineCurveOnCartesianSpace(p_curve,Xi_curve,CP_curve,isNURBS_curve,numEval);
plot3(Xp,Yp,Zp,'Color',color_curve,'LineWidth',lineWidth);
plot_ControlPolygonBSplineCurve(CP_curve);
hold off;

%% 4. Update graphs index
index = graph.index + 1;

end
