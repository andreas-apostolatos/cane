function [dN, areaTri, isInside] = ...
    computeCST2DBasisFunctionsAndFirstDerivatives ...
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
% Returns the triangular shape functions and their derivatives with respect 
% to x and y coordinates given the three vertices of the triangle and the 
% location on the x-y plane. In addition it returns the element area of the 
% given triangle allows to compute for multiple elements, so the input has 
% one additional non-singleton dimension. The vertices must be provided in 
% a counterclock-wise fashion:
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
%    vertexI,vertexJ,vertexK : The coordinates of the vertices in a 
%                              counterclockwise fashion
%                        x,y : The physical location on the x,y plane where
%                              to compute the shape functions
%
%                     Output :
%                         dN : The evaluated basis functions and their 
%                              first derivatives  at 
%                              x,y: dN = [Ni dNi/dx dNi/dy
%                                         Nj dNj/dx dNj/dy
%                                         Nk dNk/dx dNk/dy]
%                    areaTri : The area of the triangular element
%                   isInside : Flag indicating whether the point on where
%                              the basis functions are to be evaluated is
%                              inside or outside the triangle
%
% Function layout :
%
% 1. Compute the area of the triangle
%
% 2. Compute the permutations
%
% 3. Compute the basis functions for the linear triangle at (x,y)
%
% 4. Compute the derivatives of the basis functions for the linear triangle w.r.t. to x at (x,y)
%
% 5. Compute the derivatives of the basis functions for the linear triangle w.r.t. to y at (x,y)
%
% 6. Assemble to the vector containing all the basis functions and their derivatives at (x,y)
% 
% 7. Re-shape the output for singleton array
%
%% Function main body

%% 0. Read input

% Number of elemements
noElmnts = size( vertexI, 1 );

% Initialize output flag
isInside = true(noElmnts, 1, 1);

%% 1. Compute the area of the triangle
% area = .5*p3Ddet([1 vertexI(1,1) vertexI(1,2);
%                   1 vertexJ(1,1) vertexJ(1,2);
%                   1 vertexK(1,1) vertexK(1,2)]);
areaTri = 0.5*p3Ddet(phorzcat(ones(noElmnts, 3, 1), pvertcat(ptranspose(vertexI(:, 1:2)), ...
                                                             ptranspose(vertexJ(:, 1:2)), ...
                                                             ptranspose(vertexK(:, 1:2)))));

%% 2. Compute the permutations

% For basis function Ni:
% zi:
zi = vertexJ(:,1,1).*vertexK(:,2,1)-vertexK(:,1,1).*vertexJ(:,2,1);
% yjk:
yjk = vertexJ(:,2,1)-vertexK(:,2,1);
% xkj:
xkj = vertexK(:,1,1)-vertexJ(:,1,1);

% For basis function Nj:
% zj:
zj = (vertexK(:,1,1).*vertexI(:,2,1)-vertexI(:,1,1).*vertexK(:,2,1));
% yik:
yik = -(vertexI(:,2,1)-vertexK(:,2,1));
% xki:
xki = -(vertexK(:,1,1) - vertexI(:,1,1));

% For basis function Nj:
% zk:
zk = vertexI(:,1,1).*vertexJ(:,2,1)-vertexJ(:,1,1).*vertexI(:,2,1);
% yij:
yij = vertexI(:,2,1)-vertexJ(:,2,1);
% xji:
xji = vertexJ(:,1,1) - vertexI(:,1,1);

%% 3. Compute the basis functions for the linear triangle at (x,y)

% Ni:
Ni = (zi+yjk.*x+xkj.*y)./2./areaTri;

% Nj:
Nj = (zj+yik.*x+xki.*y)./2./areaTri;

% Nk:
Nk = (zk+yij.*x+xji.*y)./2./areaTri;

% Update the boolean flag
isInside( Ni < 0 | Nj < 0 | Nk < 0 ) = false;

%% 4. Compute the derivatives of the basis functions for the linear triangle w.r.t. to x at (x,y)

% Ni:
dNidx = yjk./2./areaTri;

% Nj:
dNjdx = yik./2./areaTri;

% Nk:
dNkdx = yij./2./areaTri;

%% 5. Compute the derivatives of the basis functions for the linear triangle w.r.t. to y at (x,y)

% Ni:
dNidy = xkj./2./areaTri;

% Nj:
dNjdy = xki./2./areaTri;

% Nk:
dNkdy = xji./2./areaTri;

%% 6. Assemble to the vector containing all the basis functions and their derivatives at (x,y)
dN = [cat(3, Ni, dNidx, dNidy), cat(3, Nj, dNjdx, dNjdy), cat(3, Nk, dNkdx, dNkdy)];

%% 7. Re-shape the output for singleton array
if noElmnts == 1
    dN = permute( dN, [2, 3, 1] );
end

end
