function plot_NURBS_surface_3D(p,q,U,V,CP,color,u,v,is_on_u)
%% Function documentation
%
% plots the surface, elements and control points as well as the normal and
% tangent to the specified parametric locations vectors
%
%   Input :
%     p,q : The polynomial degrees in u,v-directions
%     U,V : The knot vectors in u,v-directions
%      CP : The set of Control Point coordinates and weights
%     u,v : The parametric locations on where to compute the normal and the
%           tangent to the surface boundary vectors
% is_on_u : The parametrization of the surface boundary on which the normal 
%           and the tangent vectors will be computed:
%           0 : Surface boundary is a v-parametric line
%           1 : Surface boundary is a u-parametric line
%
%  Output :
%           Graphics
%
% Function layout :
%
% 1. On the visualization of the NURBS surface
%
% 2. On the visualization of the normal and tangent to the surface boundary vectors
%
%% 1. On the visualization of the NURBS surface

gridu = 49;
gridv = 49;

[X,Y,Z] = create_surface(p,q,U,V,CP,gridu,gridv);

% geometry
surf(X,Y,Z,'FaceColor',color,'EdgeColor','none');
hold on;

% element edges
isDeformed = 0;

create_edges(p,q,U,V,CP,isDeformed,gridu,gridv);

% control points and polygon
% create_control_polygon(CP)

% camlight left; lighting phong;
colormap('default');
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
hold off;

%% 2. On the visualization of the normal and tangent to the surface boundary vectors

% Number of Control Points at u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% Find the span indices
spanu = findspan(u,U,nu);
spanv = findspan(v,V,nv);

% Compute the the normal and the tangent to the surface boundary vectors
[n,t] = normal_and_tangent_vectors_to_the_surface_boundary(spanu,p,u,U,spanv,q,v,V,CP,is_on_u);

% Get the start and the end point of the vectors
% On the normal to the boundary vector:
N_start = point_on_surface(p,spanu,u,U,q,spanv,v,V,CP);
N_end = N_start + n';
% On the tangent to the boundary vector:
T_start = point_on_surface(p,spanu,u,U,q,spanv,v,V,CP);
T_end = N_start + t';

hold on;
% Plot the boundary normal vector:
plot_vector(N_start,N_end);
hold on;
% Plot the boundary tangent vector:
plot_vector(T_start,T_end);

hold off;

end