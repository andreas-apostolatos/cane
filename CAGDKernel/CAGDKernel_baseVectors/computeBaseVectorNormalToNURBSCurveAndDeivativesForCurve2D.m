function [G,dG] = computeBaseVectorNormalToNURBSCurveAndDeivativesForCurve2D(knotSpanIndex,p,CP,dR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns two arrays; G is the array containing the base vector of the
% NURBS curve at the given curve parameter and the not normalized
% perpendicular to the curve vector. Array dG contains the accelaration
% vector as well as perpendicular to the acceleration vector.
%
%         Input :
% knotSpanIndex : The span where u is contained
%             p : The polynomial degree of the curve
%            CP : The Control points of the NURBS curve
%            dR : The basis functions and their derivatives at a given 
%                 parametric location
%
%        Output :
%             G : 3x2 array containing the coordinates of the tangent to 
%                 the curve base vector and the normal to the curve vector
%            dG : 3x2 array containing the accelaration vector and the 
%                 normal to the acceleration vector 
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the base vector and the accelration vector
%
% 2. Compute the perpedicular to the base and accelration vector
%
%% Function main body

%% 0. Read input

% Tangent vector
G = zeros(3,2);

% Acceleration vector
dG = zeros(3,2);

%% 1. Compute the base vector and the accelration vector
for b = 0:p
    for a = 1:3
        % Compute the index
        index = knotSpanIndex - p + b;
        
        % Compute the tangent base vector to the curve
        G(a,1) = G(a,1) + dR(b+1,2)*CP(index,a);
        
        % Compute the accelaration vector to the curve
        dG(a,1) = dG(a,1) + dR(b+1,3)*CP(index,a);
    end
end
 
%% 2. Compute the perpedicular to the base and accelration vector

% Normal vector with respect to xy-axis only
G(1,2) = - G(2,1) ;
G(2,2) = G(1,2) ;
G(3,2) = G(3,1) ;

% Normal acceleration vector
dG(1,2) = - dG(2,1) ;
dG(2,2) = dG(1,2) ;
dG(3,2) = dG(3,1) ;

end


  
  
