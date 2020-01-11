function plot_reference_configuration(p,q,U,V,CP,rb,f)
% Plots the geometry together with the boundary conditions and loads arrows
%
%   Input : 
%     p,q : The polynomial degrees in u-,v-directions
%     U,V : The knot vectors in u-,v- directions
%      CP : The set of the Control points and weights
%      rb : The set of boundary conditions
%       f : The force vector
%
%  Output : 
%           graphics
%
% Function layout :
%
%% Function main body

% Create arrays needed for the plotting of the surface, the supports as
% well as the force arrows
gridu = 49;
gridv = 49;
[Xp,Yp,Zp] = create_surface(p,q,U,V,CP,gridu,gridv);
[xs,ys,zs] = create_supports(CP,rb);
[xf,yf,zf] = create_force_arrows(CP,f);

% Assign an index to the figure
% figure(graph.index)

% geometry
surf(Xp,Yp,Zp,'FaceColor','green','EdgeColor','none');
hold on;
 
% element edges
isDeformed = 0;
create_edges(p,q,U,V,CP,isDeformed,gridu,gridv)
  
%supports
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end
  
% load arrows  
for k =1:length(xf(:,1))
    plot3(xf(k,:),yf(k,:),zf(k,:),'Linewidth',5);
    plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end
  
% control points and polygon
create_control_polygon(CP)

axis equal;
% view(0,240);
view(30,15)
camlight left; lighting phong;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
% xlim([-4 12]);
% ylim([0 1]);
% zlim([-8 1]);
view(50,20);
%grid off;
% set(gca,'xtick',[-100 100]);
% set(gca,'ytick',[-100 100]);
% set(gca,'ztick',[-100 100]);
% axis off;
 
% hold off;

end