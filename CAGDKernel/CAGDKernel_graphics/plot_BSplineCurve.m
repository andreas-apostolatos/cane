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
function plot_BSplineCurve(p,Xi,CP,isNURBS,noEval,color,lineWidth)
%% Function documentation
%
% Plots a B-Spline curve.
%
%     Input :
%          p: Polynomial order
%        Xi : The knot vector
%        CP : Set of control points and weights
%    noEval : Number of sampling points to be used
%     color : The color of the line to be plotted
% lineWidth : The width of the line
%
%   Output : graphics
%
% Function layout :
%
% 1. Get the coordinates of the sampling points on the curve
%
% 2. Create the geometry
%
% 3. Plot the Control Polygon
%
% 4. Define graph properties
%           
%% Function main body

%% 1. Get the coordinates of the sampling points on the curve
[Xp,Yp,Zp] = createBSplineCurveOnCartesianSpace...
    (p,Xi,CP,isNURBS,noEval);

%% 2. Create the geometry
line(Xp,Yp,Zp,'Linewidth',lineWidth,'color',color);
hold on;

%% 3. Plot the Control Polygon
plot_ControlPolygonBSplineCurve(CP);

%% 4. Define graph properties
axis equal;

end