function index = plot_two_patches_reference(p1,q1,U1,V1,CP1,p2,q2,U2,V2,CP2,rb1,rb2,fl1,fl2,ucoup1,vcoup1,ucoup2,vcoup2,graph)
%% Function documentation
% Plots two NURBS patches in the reference configuration
%
%         Input : 
%   p1,q1,p2,q2 : polynomial degrees
%   U1,V1,U2,V2 : knot vectors in U,V-directions
%       CP1,CP2 : set of control points and weights
%       rb1,rb2 : Structures containing information on the supports
%       fl1,fl2 : Force vectors corresponding to each membrane structure
%   ucoup,vcoup : Coupling region
%         graph : structure containing all information on the plots
%
%        Output :
%         index : the index of the graph
%
%% Function main body

% Booleans
IsUndeformed = 0;

% Number the current plot
figure(graph.index)

% Plot the fisrt NURBS patch
[Xp1,Yp1,Zp1] = create_surface(p1,q1,U1,V1,CP1,50,50);
surf(Xp1,Yp1,Zp1,'FaceColor','green','EdgeColor','none');
hold;
create_edges(p1,q1,U1,V1,CP1,IsUndeformed,50,50)
% control points and control polygon
create_control_polygon(CP1)
% Create the supports
[xs1,ys1,zs1] = create_supports(CP1,rb1);
% Plot the supports
for k =1:length(xs1(:,1))
    plot3(xs1(k,:),ys1(k,:),zs1(k,:),'Linewidth',2,'Color','black');
end
% Create force arrows
[xf1,yf1,zf1] = create_force_arrows(CP1,fl1);
    
% Plot the load arrows
for k =1:length(xf1(:,1))
    plot3(xf1(k,:),yf1(k,:),zf1(k,:),'Linewidth',5);
    plot3(xf1(k,1),yf1(k,1),zf1(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

% On the visualization of the coupling interface
% Read input
nu1 = length(CP1(:,1,1));
nv1 = length(CP1(1,:,1));
nu2 = length(CP2(:,1,1));
nv2 = length(CP2(1,:,1));

% Number of sampling points
npoints = 49;

% Check in which parametric line we are on for both knot vectors:
% for patch 1
if vcoup1(1)==vcoup1(2)
    % The coupling line is on a u-curve
    uv1 = vcoup1(1);
    
    % The span in v-direction is fixed
    spanv1 = findspan(uv1,U1,nv1);
    
    % Initialize parametric coordinate
    u1 = ucoup1(1);
    
    % Compute the incremental step
    du1 = (ucoup1(2)-ucoup1(1))/npoints;
else
    % The coupling line is on a v-curve
    uv1 = ucoup1(1);
    
    % The span in u-direction is fixed
    spanu1 = findspan(uv1,U1,nu1);
    
    % Initialize parametric coordinate
    v1 = vcoup1(1);
    
    % Compute the incremental step
    dv1 = (vcoup1(2)-vcoup1(1))/npoints;
end

% for the patch 2
if vcoup2(1)==vcoup2(2)
    % The coupling line is on a u-curve
    uv2 = vcoup2(1);
    
    % The span in v-direction is fixed
    spanv2 = findspan(uv2,U2,nv2);
    
    % Initialize parametric coordinate
    u2 = ucoup2(1);
    
    % Compute the incremental step
    du2 = (ucoup2(2)-ucoup2(1))/npoints;
else
    % The coupling line is on a v-curve
    uv2 = ucoup2(1);
    
    % The span in u-direction is fixed
    spanu2 = findspan(uv2,U2,nu2);
    
    % Initialize parametric coordinate
    v2 = vcoup2(1);
    
    % Compute the incremental step
    dv2 = (vcoup2(2)-vcoup2(1))/npoints;
end

% Initialize arrays
curve1 = zeros(npoints+1,3);
curve2 = zeros(npoints+1,3);

% Loop over all the sampling points of the coupling curve for patch 1
for i=1:npoints+1
    if vcoup1(1)==vcoup1(2)
        % Find the knot span where we are inside
        spanu1 = findspan(u1,U1,nu1);
        
        % Compute the point on the curve
        curve1(i,1:3) = point_on_surface(p1,spanu1,u1,U1,q1,spanv1,uv1,V1,CP1);
              
        % Update coordinate
        u1 = u1 + du1;
    else
        % Find the knot span where we are inside
        spanv1 = findspan(v1,V1,nv1);
        
        % Compute the point on the curve
        curve1(i,1:3) = point_on_surface(p1,spanu1,uv1,U1,q1,spanv1,v1,V1,CP1);
        
        % Update coordinate
        v1 = v1 + dv1;
    end
end


% Plot the coupling interface for patch 1
plot3(curve1(:,1),curve1(:,2),curve1(:,3),'LineWidth',2,'Color','magenta');

% Plot the second NURBS patch
[Xp2,Yp2,Zp2] = create_surface(p2,q2,U2,V2,CP2,50,50);
surf(Xp2,Yp2,Zp2,'FaceColor','green','EdgeColor','none');
create_edges(p2,q2,U2,V2,CP2,IsUndeformed,50,50)
% control points and control polygon
create_control_polygon(CP2)
% Create the supports
[xs2,ys2,zs2] = create_supports(CP2,rb2);
% Plot the supports
for k =1:length(xs2(:,1))
    plot3(xs2(k,:),ys2(k,:),zs2(k,:),'Linewidth',2,'Color','black');
end
%Create force arrows
[xf2,yf2,zf2] = create_force_arrows(CP2,fl2);
    
% Plot the load arrows
for k =1:length(xf2(:,1))
    plot3(xf2(k,:),yf2(k,:),zf2(k,:),'Linewidth',5);
    plot3(xf2(k,1),yf2(k,1),zf2(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

% Loop over all the sampling points of the coupling curve for patch 2
for i=1:npoints+1
    if vcoup2(1)==vcoup2(2)
        % Find the knot span where we are inside
        spanu2 = findspan(u2,U2,nu2);
        
        % Compute the point on the curve
        curve2(i,1:3) = point_on_surface(p2,spanu2,u2,U2,q2,spanv2,uv2,V2,CP2);
              
        % Update coordinate
        u2 = u2 + du2;
    else
        % Find the knot span where we are inside
        spanv2 = findspan(v2,V2,nv2);
        
        % Compute the point on the curve
        curve2(i,1:3) = point_on_surface(p2,spanu2,uv2,U2,q2,spanv2,v2,V2,CP2);
        
        % Update coordinate
        v2 = v2 + dv2;
    end
end

% Plot the coupling interface for patch 2
plot3(curve2(:,1),curve2(:,2),curve2(:,3),'LineWidth',2,'Color','magenta');

% Adjust graph settings
camlight left; 
lighting phong;
view(2);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title ('original geometry');
hold off;

% Update plot index by 1
index = graph.index+1;

end

