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
function plot_BSplineSurfaceMultipatches(BSplinePatches,color)
% Function documentation
%
% Plots a multipatch B-Spline surface in 3D.
%
%   	   Input :
% BSplinePatches : A structure containing an array of B-Spline patches
%     		   BSplinePatches{i}.p,.q : The polynomial degrees in 
%	  			 	    xi-,eta-directions
%  		BSplinePatches{i}.Xi,.Eta : The knot vectors in xi-,
%				 	    eta-directions
%      		     BSplinePatches{i}.CP : The Control Point coordinates
%					    and weights in xi-,eta-
%					    directions
% 		BSplinePatches{i}.isNURBS : Flag on whether the basis is a 
%					    NURBS or a B-Spline
%   	   color : The color of the surface to be plotted
% 
%  Output :
%           graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over the patches
% ->
%    1i. Get the patch properties
%
%   1ii. Create the B-Spline surface
%
%  1iii. Plot the surface
%
%   1iv. Plot the element edges
%
%    1v. Plot the Control polygon for the surface
% <-
% 
%% Function main body

%% 0. Read input

% On the grid of the graphs
xiGrid = 49;
etaGrid = 49;

% Number of patches
noPatches = length(BSplinePatches);

% element edges
isDeformed = 0;

%% 1. Loop over the patches
for counterPatches = 1:noPatches
  %% 1i. Get the patch properties
  p = BSplinePatches{counterPatches}.p;
  q = BSplinePatches{counterPatches}.q;
  Xi = BSplinePatches{counterPatches}.Xi;
  Eta = BSplinePatches{counterPatches}.Eta;
  CP = BSplinePatches{counterPatches}.CP;
  isNURBS = BSplinePatches{counterPatches}.isNURBS;
  
  %% 1ii. Create the B-Spline surface
  [Xp1,Yp1,Zp1] = createBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,xiGrid,etaGrid);

  %% 1iii. Plot the surface
  surf(Xp1,Yp1,Zp1,'FaceColor',color,'EdgeColor','none');
  hold on;

  %% 1iv. Plot the element edges
  plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isDeformed,xiGrid,etaGrid);

  %% 1v. Plot the Control polygon for the surface
  plot_ControlPolygonBSplineSurface(CP);
end

end