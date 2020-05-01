function [G1, G2] = computeBaseVectorAndNormalToNURBSCurve2D ...
    (knotSpan, p, CP, dR)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the base vectors of a NURBS curve at the parametric position u
% together with the normal to that position vector which is not normalized.
%
%    Input :
% knotSpan : The span where xi is contained
%        p : The polynomial degree of the curve
%       CP : The Control points of the NURBS curve
%       dR : The basis functions and their derivatives. To get the first
%            derivatives use R(i,2), i=1,...,nBasisFcns
%
%   Output :
%       G1 : The tangent to the curve base vector G1
%       G2 : The normal to the curve vector such that G2(x) = - G1(y) and
%            G2(y) = G1(x)
%
% Function layout
%
% 1. Compute the base vector G1
%
% 2. Compute the not normalized normal to the 2D curve vector G2
%
%% Function main body

%% 1. Compute the base vector G1
G1 = zeros(2,1);
for b = 0:p
	for a = 1:2
        G1(a,1) = dR(b+1,2)*CP(knotSpan-p+b,a) + G1(a,1);
	end
end
 
%% 2. Compute the not normalized normal to the 2D curve vector G2
G2 = zeros(2,1) ;
G2(1) = - G1(2) ;
G2(2) = G1(1) ;

end


  
  
