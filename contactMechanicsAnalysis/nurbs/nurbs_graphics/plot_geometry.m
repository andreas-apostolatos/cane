%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          =============                                %
%                          || AUTHORS ||                                %
%                          =============                                %
%                                                                       %
%            Efstathios E. Theotokoglou, Professor NTUA                 %
%                   (stathis@central.ntua.gr)                           %
%                                                                       %
%            Andreas Apostolatos, Research Associate TUM                %
%                 (andreas.apostolatos@tum.de)                          %
%                                                                       %
%                       ===================                             %
%                       || plot_geometry ||                             %
%                       ===================                             %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = plot_geometry(p,q,U,V,CP,graph)
% Plots unmeshed geometry
%
%   Input : 
%     p,q : The polynomial degrees in u-,v-directions
%     U,V : The knot vectors in u-,v- directions
%      CP : The set of the Control points and weights
%   graph : Information on the graphics
%
%  Output : 
%   index : The index of the current graph
%
% Function layout :
%
%% Function main body

% Create arrays needed for the plotting of the surface, the supports as
% well as the force arrows
gridu = 99;
gridv = 99;
[Xp,Yp,Zp] = create_surface(p,q,U,V,CP,gridu,gridv);

% Assign an index to the figure
% figure(graph.index)

% geometry
surf(Xp,Yp,Zp,'FaceColor','green','EdgeColor','none');
hold on;
 
% element edges
isDeformed = 0;
create_edges(p,q,U,V,CP,isDeformed,gridu,gridv)
  
% control points and polygon
create_control_polygon(CP)
axis equal off;
axis on;
view(2);

camlight left; lighting phong;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
% xlim([-4 12]);
% ylim([0 1]);
% zlim([-8 1]);
% view(50,20);
%grid off;
% set(gca,'xtick',[-100 100]);
% set(gca,'ytick',[-100 100]);
% set(gca,'ztick',[-100 100]);
% axis off;
 
% hold off;

index = graph.index + 1;

end