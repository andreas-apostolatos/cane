function plot_knotsForBSplineCurveOnCartesianSpace(p, Xi, CP, isNURBS)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Draws the element edges for the B-Spline curve, i.e the knots on the
% geometry
%
%      Input :
%          p : Polynomial degrees
%         Xi : Knot vector
%         CP : Set of Control Point coordinates and weights
%    isNURBS : Flag on whether the basis is a NURBS or a B-Spline
%
%     Output : 
%              Graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the element edges (knots)
%
%    1i. Compute the curve parameter
%
%   1ii. Find the current knot span
%
%  1iii. Compute the NURBS basis functions at xi
%
%   1iv. Compute the Cartesian coordinates of the point on space
%
% 2. Plot all the element edges
%         
%% Function main body

%% 0. Read input

% Number of Control Points in xi,eta-direction
nxi = length(CP(:, 1, 1));

% Initialize array to be plotted
Point = zeros(length(Xi),3);

%% 1. Loop over all the element edges (knots)
for j = 1:length(Xi)
    %% 1i. Compute the curve parameter
    xi = Xi(j);
    
    %% 1ii. Find the current knot span
    knotSpan = findKnotSpan(xi, Xi, nxi);
    
    %% 1iii. Compute the NURBS basis functions at xi
    numDeriv = 0;
    R = computeIGABasisFunctionsAndDerivativesForCurve ...
        (knotSpan, p, xi, Xi, CP, isNURBS, numDeriv);
    
    %% 1iv. Compute the Cartesian coordinates of the point on space
    Point(j, 1:3) = computeCartesianCoordinatesOfAPointOnBSplineCurve ...
        (p, knotSpan, xi, Xi, CP, R);
end

%% 2. Plot all the element edges
plot3(Point(:, 1), Point(:, 2), Point(:, 3), '.', 'Markersize', 10, 'Color', 'green');
axis equal;
grid on;
xlabel('x', 'FontSize', 18);
ylabel('y', 'FontSize', 18);
zlabel('z', 'FontSize', 18);

end
