function [dN, Volume, isInside] = ...
    computeCST3DBasisFunctionsAndFirstDerivatives ...
    (vertexI, vertexJ, vertexK, vertexL, x, y, z)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tetrahedral shape functions and their derivatives with respect 
% to x, y and z coordinates given the four vertices of the tetrahedron and the 
% location in the x-y-z space. In addition it returns the element volume of the 
% given tetrahedron. The vertices must be provided in a counterclock-wise 
% fashion:
%
%    top view:                    side view:   
%
%               k                             l
%              /|\                           / \\
%             / | \                         /   \ \
%            /  |  \                       /     \  \
%           /   l   \                     /       \   \
%          /  /   \  \                   /         \   k
%         / /       \ \                 /           \ /
%        i-------------j               i-------------j
%
%                      Input :
%            vertexI,vertexJ,
%            vertexK,vertexZ : The coordinates of the vertices in a 
%                              counterclockwise fashion
%                      x,y,z : The physical location in the x,y,z space 
%                              where to compute the shape functions
%
%                     Output :
%                         dN : The evaluated basis functions and their 
%                              first derivatives  at 
%                              x,y,z: dN = [Ni dNi/dx dNi/dy dNi/dz
%                                           Nj dNj/dx dNj/dy dNj/dz
%                                           Nk dNk/dx dNk/dy dNk/dk
%                                           Nl dNl/dx dNl/dy dNl/dk]
%                     Volume : The volume of the tetrahedral element
%                   isInside : Flag indicating whether the point on where
%                              the basis functions are to be evaluated is
%                              inside or outside the tetrahedron
%
% Function layout :
%
% 1. Compute the volume of the tetrahedron
%
% 2. Compute the permutations
%
% 3. Compute the basis functions for the linear tetrahedron at (x,y,z)
%
% 4. Compute the derivatives of the basis functions for the linear tetrahedron w.r.t. to x at (x,y,z)
%
% 5. Compute the derivatives of the basis functions for the linear tetrahedron w.r.t. to y at (x,y,z)
%
% 5. Compute the derivatives of the basis functions for the linear tetrahedron w.r.t. to z at (x,y,z)
%
% 7. Assemble to the vector containing all the basis functions and their derivatives at (x,y,z)
%
%% Function main body

%% 0. Read input

nElem = size( vertexI, 1 );

% Initialize output flag
isInside = true(nElem, 1, 1);

vertexI = ptranspose( vertexI );
vertexJ = ptranspose( vertexJ );
vertexK = ptranspose( vertexK );
vertexL = ptranspose( vertexL );

%% 1. Compute the volume of the tetrahedron
%Volume = 1/6 * det([vector_a; vector_b; vector_c]) 
vector_a = ptranspose([(vertexJ(:,1)-vertexI(:,1)) (vertexJ(:,2)-vertexI(:,2)) (vertexJ(:,3)-vertexI(:,3))]);
vector_b = ptranspose([(vertexK(:,1)-vertexI(:,1)) (vertexK(:,2)-vertexI(:,2)) (vertexK(:,3)-vertexI(:,3))]);
vector_c = ptranspose([(vertexL(:,1)-vertexI(:,1)) (vertexL(:,2)-vertexI(:,2)) (vertexL(:,3)-vertexI(:,3))]);
Volume = 1/6 * p3Ddet( pvertcat( vector_a, vector_b, vector_c ) );

%% 2. Compute the permutations

% dummy for the left column of the determinant computations
one = ones( nElem, 1 );

% constant factor for the following multiplications
oneOverSixVolume = ( 1 ./ (6*Volume) );

% For basis function Ni: 
% ai: 
ai = + pstimes( p3Ddet( pvertcat( vertexJ, vertexK, vertexL ) ), oneOverSixVolume );
% bi: 
bi = - pstimes( p3Ddet( pvertcat( phorzcat( one, vertexJ(:,2), vertexJ(:,3) ), ...
                                  phorzcat( one, vertexK(:,2), vertexK(:,3) ), ...
                                  phorzcat( one, vertexL(:,2), vertexL(:,3) ) ) ), ...
                oneOverSixVolume );
% ci:
ci = + pstimes( p3Ddet( pvertcat( phorzcat( one, vertexJ(:,1), vertexJ(:,3) ), ...
                                  phorzcat( one, vertexK(:,1), vertexK(:,3) ), ...
                                  phorzcat( one, vertexL(:,1), vertexL(:,3) ) ) ), ...
                oneOverSixVolume );
% di:
di = - pstimes( p3Ddet( pvertcat( phorzcat( one, vertexJ(:,1), vertexJ(:,2) ), ...
                                  phorzcat( one, vertexK(:,1), vertexK(:,2) ), ...
                                  phorzcat( one, vertexL(:,1), vertexL(:,2) ) ) ), ...
                oneOverSixVolume );
                                    
% For basis function Nj:
% aj: 
aj = - pstimes( p3Ddet( pvertcat( vertexI, vertexK, vertexL ) ), oneOverSixVolume );
% bj: 
bj = + pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,2), vertexI(:,3) ), ...
                                  phorzcat( one, vertexK(:,2), vertexK(:,3) ), ...
                                  phorzcat( one, vertexL(:,2), vertexL(:,3) ) ) ), ...
                oneOverSixVolume );
% cj:
cj = - pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,1), vertexI(:,3) ), ...
                                  phorzcat( one, vertexK(:,1), vertexK(:,3) ), ...
                                  phorzcat( one, vertexL(:,1), vertexL(:,3) ) ) ), ...
                oneOverSixVolume );
                                    
% dj:
dj = + pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,1), vertexI(:,2) ), ...
                                  phorzcat( one, vertexK(:,1), vertexK(:,2) ), ...
                                  phorzcat( one, vertexL(:,1), vertexL(:,2) ) ) ), ...
                oneOverSixVolume );

% For basis function Nk:
% ak: 
ak = + pstimes( p3Ddet( pvertcat( vertexI, vertexJ, vertexL ) ), oneOverSixVolume );
% bk: 
bk = - pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,2), vertexI(:,3) ), ...
                                  phorzcat( one, vertexJ(:,2), vertexJ(:,3) ), ...
                                  phorzcat( one, vertexL(:,2), vertexL(:,3) ) ) ), ...
                oneOverSixVolume );
% ck:
ck = + pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,1), vertexI(:,3) ), ...
                                  phorzcat( one, vertexJ(:,1), vertexJ(:,3) ), ...
                                  phorzcat( one, vertexL(:,1), vertexL(:,3) ) ) ), ...
                oneOverSixVolume );
% dk:
dk = - pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,1), vertexI(:,2) ), ...
                                  phorzcat( one, vertexJ(:,1), vertexJ(:,2) ), ...
                                  phorzcat( one, vertexL(:,1), vertexL(:,2) ) ) ), ...
                oneOverSixVolume );
                                    
% For basis function Nl: 
% al: 
al = - pstimes( p3Ddet( pvertcat( vertexI, vertexJ, vertexK ) ), oneOverSixVolume );

% bl: 
bl = + pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,2), vertexI(:,3) ), ...
                                  phorzcat( one, vertexJ(:,2), vertexJ(:,3) ), ...
                                  phorzcat( one, vertexK(:,2), vertexK(:,3) ) ) ), ...
                oneOverSixVolume );
% cl:
cl = - pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,1), vertexI(:,3) ), ...
                                  phorzcat( one, vertexJ(:,1), vertexJ(:,3) ), ...
                                  phorzcat( one, vertexK(:,1), vertexK(:,3) ) ) ), ...
                oneOverSixVolume );
% dl:
dl = + pstimes( p3Ddet( pvertcat( phorzcat( one, vertexI(:,1), vertexI(:,2) ), ...
                                  phorzcat( one, vertexJ(:,1), vertexJ(:,2) ), ...
                                  phorzcat( one, vertexK(:,1), vertexK(:,2) ) ) ), ...
                oneOverSixVolume );

            
%% 3. Compute the basis functions for the linear tetrahedron at (x,y,z)

% Ni:
Ni = ( ai + bi .* x + ci .* y + di .* z );

% Nj:
Nj = ( aj + bj .* x + cj .* y + dj .* z );

% Nk:
Nk = ( ak + bk .* x + ck .* y + dk .* z );

% Nl:
Nl = ( al + bl .* x + cl .* y + dl .* z );

isInside( Ni < 0 | Nj < 0 | Nk < 0 | Nl < 0 ) = false;

%% 7. Assemble to the vector containing all the basis functions and their derivatives at (x,y,z)
dN = pvertcat( phorzcat( Ni, bi, ci, di ), ...
               phorzcat( Nj, bj, cj, dj ), ...
               phorzcat( Nk, bk, ck, dk ), ...
               phorzcat( Nl, bl, cl, dl ) );

end
