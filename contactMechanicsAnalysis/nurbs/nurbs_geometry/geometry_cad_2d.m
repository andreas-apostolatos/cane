function [x,y] = geometry_cad_2d(p,U,q,V,CP,u_division,v_division)
%% Returns the B-rep of a geometry created using NURBS
%
%
%% Function main body

%% 0. Read input

% Number of Control Points in u-direction
nu = length(CP(:,1,1));

% Number of Control Points in v-direction
nv = length(CP(1,:,1));

% differential increment in u-direction
du = (U(length(U))-U(1))/u_division;

% differential increment in u-direction
dv = (V(length(V))-V(1))/v_division;

% Number of boundary nodes
no_nodes = 2*(u_division+v_division);

% Initilize node array
nodes = zeros(no_nodes,3);

% Initialize node counter
nodes_counter = 1;

%% Get the boundary nodes in (u,v)=(:,0) parametric line

% initialize coordinates
u = U(1);
v = V(1);

% Loop over all u-divisions
for i=1:u_division+1
    % Find the span index in u-direction
    spanu = findspan(u,U,nu);

    % Find the span index in u-direction
    spanv = findspan(v,V,nv);

    % Assign the coordinates of the current node on the mesh
    nodes(nodes_counter,:) = point_on_surface(p,spanu,u,U,q,spanv,v,V,CP);
    
    % Update the location on the u-parametric line
    u = u + du;
    
    % Update node counter
    nodes_counter = nodes_counter + 1;
end

%% Get the boundary nodes in (u,v)=(1,:) parametric line

% initialize coordinates
u = U(length(U));
v = V(1)+dv;

% Loop over all u-divisions
for i=1:v_division
    % Find the span index in u-direction
    spanu = findspan(u,U,nu);

    % Find the span index in u-direction
    spanv = findspan(v,V,nv);

    % Assign the coordinates of the current node on the mesh
    nodes(nodes_counter,:) = point_on_surface(p,spanu,u,U,q,spanv,v,V,CP);
    
    % Update the location on the u-parametric line
    v = v + dv;
    
    % Update node counter
    nodes_counter = nodes_counter + 1;
end

%% Get the boundary nodes in (u,v)=(:,1) parametric line

% Assign a tolerance value
tolerance = 1e-8;

% initialize coordinates
u = U(length(U))-du;
v = V(length(V));

% Loop over all u-divisions
for i=1:u_division
    % Find the span index in u-direction
    spanu = findspan(u,U,nu);

    % Find the span index in u-direction
    spanv = findspan(v,V,nv);

    % Assign the coordinates of the current node on the mesh
    nodes(nodes_counter,:) = point_on_surface(p,spanu,u,U,q,spanv,v,V,CP);
    
    % Update the location on the u-parametric line
    u = u - du;
    if abs(u)<tolerance
        u = 0;
    end
    
    % Update node counter
    nodes_counter = nodes_counter + 1;
end

%% Get the boundary nodes in (u,v)=(0,:) parametric line

% initialize coordinates
u = U(1);
v = V(length(V))-dv;

% Loop over all u-divisions
for i=1:v_division-1
    % Find the span index in u-direction
    spanu = findspan(u,U,nu);

    % Find the span index in u-direction
    spanv = findspan(v,V,nv);

    % Assign the coordinates of the current node on the mesh
    nodes(nodes_counter,:) = point_on_surface(p,spanu,u,U,q,spanv,v,V,CP);
    
    % Update the location on the u-parametric line
    v = v - dv;
    if abs(v)<tolerance
        v = 0;
    end
    
    % Update node counter
    nodes_counter = nodes_counter + 1;
end

%% Assign the  output arrays
x = nodes(:,1);
y = nodes(:,2);

end

