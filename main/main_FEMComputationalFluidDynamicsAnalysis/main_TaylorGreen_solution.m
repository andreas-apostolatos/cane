% define square boundary
discretization = 99;

% define coordinates
x_discretization = 0:discretization;
y_discretization = 0:discretization;

% number of values
noValues = length(x_discretization);

% problem parameters
size = 2*pi;
t = 0.1;
ro = 1;
nu = 1;

% Discrete coordinates
x_coordinates = size * (x_discretization/discretization);
y_coordinates = size * (y_discretization/discretization);

% Initialize mesh nodes
nodes = zeros(noValues*noValues,2);

% assign mesh nodes from discrete coordinates
k = 1;
for n = 1:noValues
    for m = 1:noValues
        % Assign coordinates to nodes
        nodes(k,1) = x_coordinates(m);
        nodes(k,2) = y_coordinates(n);
        % Updae counter
        k=k+1;
    end
end

% create reference solution
noNodes = length(nodes);

% initialize solutions
u = zeros(noNodes,1);
v = zeros(noNodes,1);
p = zeros(noNodes,1);

% compute solution for each node
for g = 1:noNodes
    % get coodinates
    x = nodes(g,1);
    y = nodes(g,2);
    
    % compute solution
    u(g) = -exp(-2*nu*t)*cos(x)*sin(y);
    v(g) = exp(-2*nu*t)*sin(x)*cos(y);
    p(g) = -ro*0.25*( (cos(2*x)+cos(2*y))*exp(-4*nu*t) );
end

% reshape variables
u = reshape(u,noValues,noValues);
v = reshape(v,noValues,noValues);
p = reshape(p,noValues,noValues);

% plot pressure solution
figure('Name','X-velocity');
pcolor(x_coordinates,x_coordinates,u);
title(['Analytical solution of U at t = ',num2str(t)]);
shading interp;
axis equal;
axis tight;
xlabel('----X---->');
ylabel('----Y---->');
colorbar('Location','EastOutside');

% plot pressure solution
figure('Name','Y-velocity');
pcolor(x_coordinates,x_coordinates,v);
title(['Analytical solution of V at t = ',num2str(t)]);
shading interp;
axis equal;
axis tight;
xlabel('----X---->');
ylabel('----Y---->');
colorbar('Location','EastOutside');

% plot pressure solution
figure('Name','Pressure');
pcolor(x_coordinates,x_coordinates,p);
title(['Analytical solution of Pressure at t = ',num2str(t)]);
shading interp;
axis equal;
axis tight;
xlabel('----X---->');
ylabel('----Y---->');
colorbar('Location','EastOutside');




