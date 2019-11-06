function shear_y = load_shear_y(p,i,u,U,q,j,v,V,CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns Neumann boundary condition that has been recovered from 
% analytical solution for the recangular plate fixed at one edge and
% subject to beding moment at the other end
%
%  Input : 
%    p,q : Polynomial degrees
%    i,j : Knot span indeces
%    u,v : Parametric coordinates on the surface
%     CP : The set of control points and weights
%
% Output :
%  shear_y : The value of the analytical load at the surface location (u,v) 

%% Function main body

% The load maximum amplitude
Q = 1e4;

% The height of the plate
h = 2;

% Initialize the coordinates
% x = 0;
y = 0;

% Compute the basis functions affecting the knot span
Rb = nurbs_basis_functions2D(i,p,u,U,j,q,v,V,CP);

% initialize counter
k = 0;

for c = 0:q 
    for b = 0:p
        % update counter
        k = k + 1;
        
        % compute the location on x-y plane
        % x = Rb(k)*CP(i-p+b,j-q+c,1) + x;
        y = Rb(k)*CP(i-p+b,j-q+c,2) + y;  
    end
end

% Compute the load at the point x-y
shear_y = 2*Q*y/h;

end
