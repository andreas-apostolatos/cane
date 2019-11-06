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
function [X,Y,Z] = createBSplineCurveOnCartesianSpace(p,Xi,CP,isNURBS,nEval)
%% Function documentation
%
% Returns three arrays containing the spatial positions of the points on 3D
% space belonging to the curve
%
%   Input : 
%       p : The polynomial degree of the B-Spline curve
%      Xi : The knot vector
%      CP : The set of Control Point and coordinates
% isNURBS : On whether the basis is a NURBS or a B-Spline
%   nEval : Number of evaluation points
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
%  1iii. Compute the Cartesian coordinates of xi
%
%   1iv. Update the parametric location of the point
%
% 2. Write the output
%
%% Function main body

%% 0. Read input

% The number of Knots
mxi = length(Xi);

% The number of base vectors
nxi = length(CP(:,1));

% Check input
checkInputForBSplineCurve(p,mxi,nxi);

% Compute an increment
dxi = (Xi(length(Xi)) - Xi(1))/(nEval-1);

% Initialize output array
XYZ = zeros(nEval,3);

% points in u-direction
xi = Xi(p+1);

%% 1. Loop over all the parametric locations
for j=1:nEval
    %% 1i. Find the span index where xi is located
    spandIndex = findKnotSpan(xi,Xi,nxi);
    
    %% 1ii. Compute the basis functions at xi
    nDeriv = 0;
    R = computeIGABasisFunctionsAndDerivativesForCurve...
        (spandIndex,p,xi,Xi,CP,isNURBS,nDeriv);
    
    %% 1iii. Compute the Cartesian coordinates of xi 
    XYZ(j,1:3) = computeCartesianCoordinatesOfAPointOnBSplineCurve...
        (p,spandIndex,xi,Xi,CP,R);
    
    %% 1iv. Update the parametric location of the point
    xi = xi + dxi;
end

%% 2. Write the output
X = XYZ(:,1);
Y = XYZ(:,2);
Z = XYZ(:,3);

end