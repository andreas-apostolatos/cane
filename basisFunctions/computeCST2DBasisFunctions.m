function [N, Area, isInside] = computeCST2DBasisFunctions ...
    (vertexI, vertexJ, vertexK, x, y)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the triangular shape functions given the three vertices of the
% triangle and the location on the x-y plane. The vertices must be provided
% in a counterclock-wise fashion:
%
%               k
%              / \
%             /   \
%            /     \ 
%           /       \
%          /         \
%         /           \
%        i-------------j
%
%                      Input :
%    vertexI,vertexJ,vertexK : The coordinates of the triangular vertices 
%                              counterclockwise
%                        x,y : Physical location on the x,y plane to 
%                              compute the shape functions
%
%                     Output :
%                          N : The evaluated basis functions at x,y: 
%                              N = [Ni Nj Nk]
%                       Area : The area of the triangle
%                   isInside : Flag on whether the point 
%
% Function layout :
%
% 0. Read input
% 
% 1. Compute the area of the triangle
%
% 2. Compute the permutations
%
% 3. Compute the basis functions for the linear triangle at (x,y)
%
% 4. Return a flag on whether the point was in or out the triangle
%
%% Function main body

%% 0. Read input

% Initialize output flag
isInside = true;

%% 1. Compute the area of the triangle

Area = .5*det([1 vertexI(1,1) vertexI(1,2)
               1 vertexJ(1,1) vertexJ(1,2)
               1 vertexK(1,1) vertexK(1,2)]);

%% 2. Compute the permutations

% For basis function Ni:
% zi:
zi = vertexJ(1,1)*vertexK(1,2)-vertexK(1,1)*vertexJ(1,2);
% yjk:
yjk = vertexJ(1,2)-vertexK(1,2);
% xkj:
xkj = vertexK(1,1) - vertexJ(1,1);

% For basis function Nj:
% zj:
zj = (vertexK(1,1)*vertexI(1,2)-vertexI(1,1)*vertexK(1,2));
% yik:
yik = -(vertexI(1,2)-vertexK(1,2));
% xki:
xki = -(vertexK(1,1) - vertexI(1,1));

% For basis function Nj:
% zk:
zk = vertexI(1,1)*vertexJ(1,2)-vertexJ(1,1)*vertexI(1,2);
% yij:
yij = vertexI(1,2)-vertexJ(1,2);
% xji:
xji = vertexJ(1,1) - vertexI(1,1);

%% 3. Compute the basis functions for the linear triangle at (x,y)

% Ni:
Ni = (zi+yjk*x+xkj*y)/2/Area;

% Nj:
Nj = (zj+yik*x+xki*y)/2/Area;

% Nk:
Nk = (zk+yij*x+xji*y)/2/Area;

% vector containing all the basis functions
N = [Ni 
     Nj 
     Nk];

%% 4. Return a flag on whether the point was in or out the triangle
if Ni < 0 || Nj < 0 || Nk < 0
    isInside = false;
end

end
