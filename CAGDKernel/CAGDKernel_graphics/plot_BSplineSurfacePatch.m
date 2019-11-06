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
function plot_BSplineSurfacePatch(p,q,Xi,Eta,CP,isNURBS,color)
% Function documentation
%
% Plots a B-Spline surface in 3D.
%
%   Input :
%     p,q : The polynomial degrees in xi-,eta-directions
%  Xi,Eta : The knot vectors in xi-,eta-directions
%      CP : The Control Point coordinate and weights in xi-,eta-directions
% isNURBS : Flag on whether the basis is a NURBS or a B-Spline
%   color : The color of the surface to be plotted
% 
%  Output :
%           graphics
%
% Function layout :
%
% 1. Create the B-Spline surface
%
% 2. Plot the surface
%
% 3. Plot the element edges
%
% 4. Plot the Control polygon for the surface
% 
%% Function main body

%% 0. Read input

% On the grid of the graphs
xiGrid = 49;
etaGrid = 49;

% element edges
isDeformed = 0;

% Dummy variables
prestress = 'undefined';
compPrestress = 'undefined';

%% 1. Create the B-Spline surface
[Xp1,Yp1,Zp1] = createBSplineSurfaceOnCartesianSpace...
    (p,q,Xi,Eta,CP,isNURBS,prestress,compPrestress,xiGrid,etaGrid);

%% 2. Plot the surface
surf(Xp1,Yp1,Zp1,'FaceColor',color,'EdgeColor','none');
hold on;

%% 3. Plot the element edges
plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isDeformed,xiGrid,etaGrid);

%% 4. Plot the Control polygon for the surface
plot_ControlPolygonBSplineSurface(CP);

end