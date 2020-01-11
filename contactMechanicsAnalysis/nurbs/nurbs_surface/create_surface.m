function [X,Y,Z] = create_surface(p,q,U,V,CP,gridu,gridv)
%% Function documentation
%
% Returns three arrays, containing the coordinates of the points on the
% NURBS surface in a grid of gridu x gridv lines
%
%   Input :
%     p,q : Polynomial degrees
%     U,V : Knot vectors in u,v-direction
%      CP : Control point coordinates and weights
%   gridu : No. of lines in u-direction
%   gridv : No. of lines in v-direction
%
%  Output :
%       X : Array containing the x-coordinates of the points on the surface
%       Y : Array containing the y-coordinates of the points on the surface
%       Z : Array containing the z-coordinates of the points on the surface
%
%% Function main body

% Number of knots in u,v-direction
mu = length(U);
mv = length(V);

% Number of control points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% Check the compatibility of the NURBS parameters
check_input(p,mu,nu,q,mv,nv)

% Compute a tolerance value
eps=10e-10;

% Initialize counter for the lines in v-direction
l=1;  

% Incremental step for the lines in v-direction
r=(V(mv)-V(1))/gridu;    
v=V(1);

% Initialize output array
XYZ = zeros(gridu,gridv,3);

% Loop over all the parametric coordinate locations
while v <= V(mv)+eps
    j = findspan(v,V,nv);
    % Incremental step for the lines in u-direction
    s=(U(mu)-U(1))/gridv;  
    u=U(1);
    % Initialize counter for the lines in u-direction
    k=1;
    while u <= U(mu)+eps
    	i = findspan(u,U,nu);
    	XYZ(k,l,1:3) = point_on_surface(p,i,u,U,q,j,v,V,CP);
    	% Update counter for the lines in u-direction
    	k=k+1;
        % Update the u-parametric coordinate
        u=u+s;
    end
  % Update counter for the lines in v-direction
  l=l+1;
  % Update the v-parametric coordinate
  v=v+r;
end

% Write the coordinates into the individual arrays
X = XYZ(:,:,1);
Y = XYZ(:,:,2);
Z = XYZ(:,:,3);

end