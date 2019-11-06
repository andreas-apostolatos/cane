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
function [X,Y,Z] = createBSplineCurveOnBSplineSurface...
    (p,q,Xi,Eta,CP,isNURBS,grid,xiStart,etaStart,xiEnd,etaEnd)
%% Function documentation
%
% Returns three arrays, containing the coordinates of the points on the
% NURBS surface in a grid of gridu x gridv lines
%
%    Input :
%      p,q : Polynomial degrees
%      U,V : Knot vectors in u,v-direction
%       CP : Set of control points and weights
%     grid : Number of sampling points to be used
%    u1,v1 : Starting point parameters of the curve
%    u2,v2 : Ending point parameters of the curve
%
%   Output :
%        X : Array containing the x-coordinates of the points on the surface
%        Y : Array containing the y-coordinates of the points on the surface
%        Z : Array containing the z-coordinates of the points on the surface
%
% Function Layout :
%
% 0. Read input 
%
% 1. Loop over all the sampling points on the curve
%
% 2. Write the coordinates into the individual arrays
%
%% Function main body

%% 0. Read input 

% Number of control points in u,v-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Initialize output array
XYZ = zeros(grid,3);

% Number of knots in u-direction
mxi = length(Xi);

% Number of knots in v-direction
meta = length(Eta);

% Compute a step size for lambda
step = 1/grid;

%% 1. Loop over all the sampling points on the curve
for i = 1:grid + 1
    % Get the linear combination factor
    lambda = (i-1)*step;
    
    % Update the parametric location
    xi = (1-lambda)*xiStart + lambda*xiEnd;
    eta = (1-lambda)*etaStart + lambda*etaEnd;
    
    % Check the input parameters
    checkInputForBSplineCurve(p,mxi,nxi);
    checkInputForBSplineCurve(q,meta,neta);
    
    % Find the correct knot span
    xiSpan = findKnotSpan(xi,Xi,nxi);
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    % Compute the IGA basis functions
    R = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
    
    % Get the point on the surface
    XYZ(i,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R);
    
end

%% 2. Write the coordinates into the individual arrays
X = XYZ(:,1);
Y = XYZ(:,2);
Z = XYZ(:,3);

end