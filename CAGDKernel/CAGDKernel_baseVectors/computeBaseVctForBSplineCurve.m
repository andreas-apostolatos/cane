function [G, H] = computeBaseVctForBSplineCurve ...
    (knotSpanIndex, p, CP, dR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the base vector and its first derivative for a B-Spline curve.
%
%         Input :
% knotSpanIndex : The span where u is contained
%             p : The polynomial degree of the curve
%            CP : The Control points of the NURBS curve
%            dR : The basis functions and their derivatives at a given 
%                 parametric location
%
%        Output :
%             G : The base vector to the B-Spline curve at the given
%                 parametric location
%             H : The curvature vector at the parametric location
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the base vector
%
%% Function main body

%% 0. Read input

% Check input
[~,m] = size(dR);

% Tangent vector
if m < 2
    error('For the computation of the base vector the first derivatives are needed');
end
G = zeros(3,1);
if m > 2
    H = zeros(3,1);
else
    H = 'undefined';
end

%% 1. Compute the base vector
for b = 0:p
    for a = 1:3
        % Compute the index
        index = knotSpanIndex - p + b;
        
        % Compute the tangent base vector to the curve
        G(a,1) = G(a,1) + dR(b+1,2)*CP(index,a);
        
        % Compute the curvature vector
        if m > 2
            H(a,1) = H(a,1) + dR(b+1,3)*CP(index,a);
        end
    end
end

end
