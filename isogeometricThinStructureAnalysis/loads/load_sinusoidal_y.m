function sinusoidal_y = load_sinusoidal_y ...
    (p, i, xi, Xi, q, j, eta, Eta, CP)
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
% subject to a sinusoidal loading at the opposite edge
%
%   Input : 
%     p,q : Polynomial degrees
%     i,j : Knot span indeces
%  xi,eta : Parametric coordinates on the surface
%  Xi,Eta : Knot vectors of the surface
%      CP : The set of control points and weights
%
%  Output :
% shear_y : The value of the analytical load at the surface location 
%           (xi,eta) 
%
%% Function main body

% The load maximum amplitude
Q = 1e3;

% The height of the plate
h = 2;

% Initialize the coordinates
% x = 0;
y = 0;

% Compute the basis functions affecting the knot span
Rb = nurbs_basis_functions2D(i,p,xi,Xi,j,q,eta,Eta,CP);

% initialize counter
k = 0;

for c = 0:q 
    for b = 0:p
        % update counter
        k = k + 1;
        
        % compute the location on x-y plane
        y = Rb(k)*CP(i-p+b,j-q+c,2) + y;  
    end
end

% Compute the load at the point x-y
% sinusoidal_y = Q*y/h - (3e-1)*Q*sin(pi*y/h); %Q*sin(pi*y/h)*cos(pi*y/h);
pot=1;
sinusoidal_y = 1e3*Q*(y^pot)/(h^pot);%cos(pi*y/h);

end

