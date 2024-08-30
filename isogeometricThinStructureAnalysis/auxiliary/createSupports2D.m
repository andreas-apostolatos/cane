function [xs, ys, zs] = createSupports2D(CP, homDOFs)
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% returns coordinates for the triangles at the support locations
%
%    Input : 
%       CP : Control Point coordinates
%  homDOFs : Array containing the global numbering of the DOFs where
%            homogeneous Dirichlet boundary conditions are applied
%
%   Output :
% xs,ys,zs : x-,y- and z-coordinates of the support triangle vertices
%
%% Function main body

% Total number of Control Points
nu = length(CP(:, 1, 1));

% scaling factors for the support triangles
up = max(max(max(max(CP))));
lo = min(min(min(min(CP))));

% Average the factor with respect to the maximum and minimum values 
fac = (up - lo)/5;

% Initialize the output arrays
xs = zeros(length(homDOFs), 4);
ys = zeros(length(homDOFs), 4);
zs = zeros(length(homDOFs), 4);

for iHom = 1:length(homDOFs)
    % Get the corresponding Control Point number p and indices CP(i,j)
    h = homDOFs(iHom)/2;
    p = ceil(h);
    j = ceil(p/nu);
    i = p - (j - 1)*nu;
    %(rb is odd -> horizontal support)
    if (p ~= h)
        xs(iHom, 1) = CP(i, j, 1);
        xs(iHom, 2) = CP(i, j, 1) - 0.1732*fac;
        xs(iHom, 3) = CP(i, j, 1) - 0.1732*fac;
        xs(iHom, 4) = xs(iHom, 1);
        ys(iHom, 1) = CP(i, j, 2);
        ys(iHom, 2) = CP(i, j, 2) + 0.1*fac;
        ys(iHom, 3) = CP(i, j, 2) - 0.1*fac;
        ys(iHom, 4) = ys(iHom, 1);
        zs(iHom, 1:4) = CP(i, j, 3);
    %(rb is even -> vertical support)
    else        
        xs(iHom, 1) = CP(i, j, 1);
        xs(iHom, 2) = CP(i, j, 1) - 0.1*fac;
        xs(iHom, 3) = CP(i, j, 1) + 0.1*fac;
        xs(iHom, 4) = xs(iHom, 1);
        ys(iHom, 1) = CP(i, j, 2);
        ys(iHom, 2) = CP(i, j, 2) - 0.1732*fac;
        ys(iHom, 3) = CP(i, j, 2) - 0.1732*fac;
        ys(iHom, 4) = ys(iHom, 1);
        zs(iHom, 1:4) = CP(i, j, 3);
    end
end

end