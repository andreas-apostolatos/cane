function create_edges(p,q,U,V,CP,isDeformed,gridu,gridv)
%% Function documentation
%
% Draws the element edges for the NURBS surface, i.e the knots on the
% geometry
%
%      Input :
%        p,q : Polynomial degrees
%        U,V : Knot vectors in u,v-direction
%         CP : Control Point coordinates and weights
% isDeformed :
%      gridu : Points to use in u-direction
%      gridv : Points to use in v-direction
%
%     Output : Graphics
%         
%% Function main body

% Number of knots in u,v-direction
mu = length(U);
mv = length(V);

% Number of Control Points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% Assign a tolerance value
eps=10e-10;

% Edges in u-direction
% Initialize counter for the edges in u-direction
l=1;  
for j2 = q+1:mv-q
    v = V(j2);
    j = findspan(v,V,nv);
    s = (U(mu)-U(1))/gridu;  
    u = U(1);
    k=1;
    while u <= U(mu)+eps
        i = findspan(u,U,nu);
        Point(k,l,1:3) = point_on_surface(p,i,u,U,q,j,v,V,CP);
        k=k+1;
        u=u+s;
    end
    l=l+1;
end

% Edges in v-direction
for i2 = p+1:mu-p
    u = U(i2);
    i = findspan(u,U,nu);
    s = (V(mv)-V(1))/gridv; 
    v = V(1);
    k = 1;
    while v <= V(mv)+eps
        j = findspan(v,V,nv);
        Point(k,l,1:3) = point_on_surface(p,i,u,U,q,j,v,V,CP);
        k = k+1;
        v = v+s;
    end
    l = l+1;
end

% Plot all the element edges
if isDeformed == 0
    plot3(Point(:,:,1),Point(:,:,2),Point(:,:,3),'Color','black');
    axis equal;
    grid on;
    xlabel('x','FontSize',18);
    ylabel('y','FontSize',18);
    zlabel('z','FontSize',18);
else
    plot3(Point(:,:,1),Point(:,:,2),Point(:,:,3),'Color','black','LineStyle','-.');
    axis equal;
    grid on;
    xlabel('x','FontSize',18);
    ylabel('y','FontSize',18);
    zlabel('z','FontSize',18);
end

end