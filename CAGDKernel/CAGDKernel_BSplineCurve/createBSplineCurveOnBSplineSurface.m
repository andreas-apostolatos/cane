function [X, Y, Z] = createBSplineCurveOnBSplineSurface ...
    (p, q, Xi, Eta, CP, isNURBS, grid, xiStart, etaStart, xiEnd, etaEnd)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns three arrays, containing the coordinates of the points on the
% NURBS surface in a grid of gridu x gridv lines
%
%            Input :
%              p,q : Polynomial degrees
%           Xi,Eta : Knot vectors in xi,eta-direction
%               CP : Set of control points and weights
%             grid : Number of sampling points to be used
% xiStart,etaStart : Starting point parameters of the curve
%     xiEnd,etaEnd : Ending point parameters of the curve
%
%           Output :
%            X,Y,Z : Arrays containing the x,y and z-coordinates of the 
%                    points on the surface
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
