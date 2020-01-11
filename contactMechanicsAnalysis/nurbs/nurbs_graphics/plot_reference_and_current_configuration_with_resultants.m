function plot_reference_and_current_configuration_with_resultants(p,q,U,V,CP1,CP2,rb,f,E,t,nue,ss,comp)
% plots two windows. The first one shows the undeformed geometry with
% loads and supports. The second one shows the strains or stresses on the
% deformed geometry.
% Parameters:
% CP1:     undeformed control points
% CP2:     deformed control points
% ss:      1=plot strain;   2=plot stress
% comp:    component of strain/stress:
%          1=xx, 2=yy, 3=xy, 4=principal_1, 5=principal_2

mu = length(U);
mv = length(V);
nu = length(CP2(:,1,1));
nv = length(CP2(1,:,1));
Ds = E/(1-nue^2)*[1 nue 0; nue 1 0; 0 0 (1-nue)/2];

%===================================================================
% 1. original geometry, supports and load arrows
[Xp1,Yp1,Zp1] = create_surf(p,q,U,V,CP1);
[xs1,ys1,zs1] = create_supports(CP1,rb);
[xf1,yf1,zf1] = create_arrows(CP1,f);

% FIRST WINDOW: UNDEFORMED
subplot(2,1,1);
hold on;
% geometry
surf(Xp1,Yp1,Zp1,'FaceColor','green','EdgeColor','none');


% element edges
create_el_edges(p,q,U,V,CP1)

% supports
for k =1:length(xs1(:,1))
  plot3(xs1(k,:),ys1(k,:),zs1(k,:),'Linewidth',2,'Color','black');
end
axis equal;

% load arrows
for k =1:length(xf1(:,1))
  plot3(xf1(k,:),yf1(k,:),zf1(k,:),'Linewidth',5);
  plot3(xf1(k,1),yf1(k,1),zf1(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

% control points and polygon
  create_conpolygon(CP1)

camlight left; lighting phong;
view(3);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title ('original geometry');
% hold off;
%=================================================================
% SECOND WINDOW: STRAIN/STRESS ON DEFORMED STRUCTURE

tol=10e-10;
l=1;  % counting index of lines
r=(V(mv)-V(1))/99;    %incremental step for v
v=V(1)+0.0001;
while v <= V(mv)+tol
  j = findspan(v,V,nv);
  s=(U(mu)-U(1))/99;  %incremental step for u
  u=U(1);
  k=1;
  while u <= U(mu)+tol
    i = findspan(u,U,nu);
    P(k,l,1:3) = get_point_surf(p,i,u,U,q,j,v,V,CP2);
    eps_actual = strain(i,p,u,U,j,q,v,V,CP1,CP2);
    eps(k,l,1:3) = eps_actual(1:3);
    eps(k,l,4) = 0.5*(eps(k,l,1)+eps(k,l,2) + sqrt((eps(k,l,1)-eps(k,l,2))^2 + eps(k,l,3)^2));
    eps(k,l,5) = 0.5*(eps(k,l,1)+eps(k,l,2) - sqrt((eps(k,l,1)-eps(k,l,2))^2 + eps(k,l,3)^2));
    n(k,l,1:3) = t*Ds*eps_actual(1:3);
    n(k,l,4) = 0.5*(n(k,l,1)+n(k,l,2) + sqrt((n(k,l,1)-n(k,l,2))^2 + 4*n(k,l,3)^2));
    n(k,l,5) = 0.5*(n(k,l,1)+n(k,l,2) - sqrt((n(k,l,1)-n(k,l,2))^2 + 4*n(k,l,3)^2));
%     sig(k,l,1:3) = Ds*eps_actual(1:3);
%     sig(k,l,4) = 0.5*(sig(k,l,1)+sig(k,l,2) + sqrt((sig(k,l,1)-sig(k,l,2))^2 + 4*sig(k,l,3)^2));
%     sig(k,l,5) = 0.5*(sig(k,l,1)+sig(k,l,2) - sqrt((sig(k,l,1)-sig(k,l,2))^2 + 4*sig(k,l,3)^2));
    k=k+1;
    u=u+s;
  end
  l=l+1;
  v=v+r;
end
% plot
subplot(2,1,2);
if     (ss==1); surf(P(:,:,1),P(:,:,2),P(:,:,3),eps(:,:,comp));
elseif (ss==2); surf(P(:,:,1),P(:,:,2),P(:,:,3),n(:,:,comp));   end
shading interp;
colormap('default');
% invert default colormap => red=negativ, blue=positive
% COL=colormap;
% invCOL(:,1)=COL(:,3);
% invCOL(:,2)=COL(:,2);
% invCOL(:,3)=COL(:,1);
% colormap(invCOL);
% make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);
% ------
colorbar;
hold on;
create_el_edges(p,q,U,V,CP2);
view(3);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if     (ss==1 && comp==1); title ('strain eps_11');
elseif (ss==1 && comp==2); title ('strain eps_22');
elseif (ss==1 && comp==3); title ('strain 2*eps_12');
elseif (ss==1 && comp==4); title ('principal strain e1');
elseif (ss==1 && comp==5); title ('principal strain e2');
elseif (ss==2 && comp==1); title ('normal force n11');
elseif (ss==2 && comp==2); title ('normal force n22');
elseif (ss==2 && comp==3); title ('inplane shear n12');
elseif (ss==1 && comp==4); title ('principal normal force n1');
elseif (ss==1 && comp==5); title ('principal normal force n2');
end
hold off;