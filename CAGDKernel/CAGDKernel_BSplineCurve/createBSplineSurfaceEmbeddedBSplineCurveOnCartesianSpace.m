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
function [X,Y,Z] = createBSplineSurfaceEmbeddedBSplineCurveOnCartesianSpace....
    (p_surface,q_surface,Xi_surface,Eta_surface,CP_surface,isNURBS_surface,...
    p_curve,Xi_curve,CP_curve,isNURBS_curve,nEval)
%% Function documentation
%
% Returns three arrays containing the spatial positions of the points on 3D
% space belonging to the curve
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
%
%          Output :
%           X,Y,Z : The Cartesian coordinates of the embedded B-Spline
%                   curve in the physical space
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the parametric locations
%
%    1i. Find the span index where xi is located
%
%   1ii. Compute the basis functions at xi
%
%  1iii. Compute the Cartesian coordinates of xi in the parametric space of the surface
%
%   1iv. Compute the Cartesian coordinates of xi in the physical space of the surface
%
%    1v. Update the parametric location of the point
%
% 2. Write the output
%
%% Function main body

%% 0. Read input

% The number of Knots
mxi = length(Xi_curve);

% The number of base vectors
nxi = length(CP_curve(:,1));

% Check input
checkInputForBSplineCurve(p_curve,mxi,nxi);

% Compute an increment
dxi = (Xi_curve(length(Xi_curve)) - Xi_curve(1))/(nEval-1);

% Initialize output array
XYZ = zeros(nEval,3);

% points in u-direction
xiTilde = Xi_curve(p_curve+1);

%% 1. Loop over all the parametric locations
for j=1:nEval
    %% 1i. Find the span index where xi is located
    spandIndex = findKnotSpan(xiTilde,Xi_curve,nxi);
    
    %% 1ii. Compute the basis functions at xi
    nDeriv = 0;
    R_curve = computeIGABasisFunctionsAndDerivativesForCurve...
        (spandIndex,p_curve,xiTilde,Xi_curve,CP_curve,isNURBS_curve,nDeriv);
    
    %% 1iii. Compute the Cartesian coordinates of xi in the parametric space of the surface
    xiEta = computeCartesianCoordinatesOfAPointOnBSplineCurve...
        (p_curve,spandIndex,xiTilde,Xi_curve,CP_curve,R_curve);
    
    %% 1iv. Compute the Cartesian coordinates of xi in the physical space of the surface
    xiSpan = findKnotSpan(xiEta(1),Xi_surface,length(CP_surface(:,1,1)));
    etaSpan = findKnotSpan(xiEta(2),Eta_surface,length(CP_surface(1,:,1)));
    R_surface = computeIGABasisFunctionsAndDerivativesForSurface...
        (xiSpan,p_surface,xiEta(1),Xi_surface,etaSpan,q_surface,xiEta(2),Eta_surface,CP_surface,isNURBS_surface,0);
    XYZ(j,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface...
        (xiSpan,p_surface,xiEta(1),Xi_surface,etaSpan,q_surface,xiEta(2),Eta_surface,CP_surface,R_surface);
    
    %% 1v. Update the parametric location of the point
    xiTilde = xiTilde + dxi;
end

%% 2. Write the output
X = XYZ(:,1);
Y = XYZ(:,2);
Z = XYZ(:,3);

end