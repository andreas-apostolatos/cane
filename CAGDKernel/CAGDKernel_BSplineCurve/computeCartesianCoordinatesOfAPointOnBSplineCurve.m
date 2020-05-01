function X = computeCartesianCoordinatesOfAPointOnBSplineCurve ...
    (p, knotSpanIndex, xi, Xi, CP, R)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the Cartesian coordinates of a point on a BSpline Curve
%
%         Input :
%             p : The polynomial degree of the curve
% knotSpanIndex : The span where u is contained
%            xi : The coordinate in the unit interval where base vectors 
%                 are to be evaluate
%            Xi : The knot vector of the NURBS curve
%            CP : The Control Points of the NURBS curve
%             R : The NURBS basis functions R(i,1), i=1,...,nBasisFncs
%
%        Output :
%             X : Vector containing the Cartesian coordinates of the NURBS 
%                 curve at xi
%
%% Function main body

% Find the correct knot span
if knotSpanIndex == 0
    knotSpanIndex = findspan(xi,Xi,length(CP(:,1)));  
end

% Initialize the point coordinates
X(1) = 0;
X(2) = 0;
X(3) = 0;

% Compute the point coordinates iteratively
for b = 0:p
    X(1) = R(b+1,1)*CP(knotSpanIndex-p+b,1) + X(1);
    X(2) = R(b+1,1)*CP(knotSpanIndex-p+b,2) + X(2);
    X(3) = R(b+1,1)*CP(knotSpanIndex-p+b,3) + X(3); 
end

end
