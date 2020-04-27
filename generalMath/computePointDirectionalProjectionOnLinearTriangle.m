function zetaLambda = computePointDirectionalProjectionOnLinearTriangle ...
    (P, X1, X2, X3, n)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the natural coordinates over the unit triangle corresponding to
% the directional projection of point P in the direction of n onto the
% triangle.
%
%       Input :
%           P : The point to be directionally projected onto the triangle
%    X1,X2,X3 : The vertices of the triangle in a counterclockwise fashion
%           n : The normal vector pointing to the direction of the projection
%
%      Output :
% zeta1,zeta2 : The coordinates of the projected point in the triangular
%               coordinates
%      lambda : The distance between the point to be projected and the
%               projected one
% isProjected : Flag on whether the projection has been performed
%               successfully
%
% Function layout :
%
%% Function main body
zetaLambda = [X1(1,1)-X3(1,1) X2(1,1)-X3(1,1) -n(1,1)
              X1(2,1)-X3(2,1) X2(2,1)-X3(2,1) -n(2,1)
              X1(3,1)-X3(3,1) X2(3,1)-X3(3,1) -n(3,1)]\[P(1,1)-X3(1,1)
                                                        P(2,1)-X3(2,1)
                                                        P(3,1)-X3(3,1)];

end

