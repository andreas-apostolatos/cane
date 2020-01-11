function plot_gauss_points_on_surface(p,q,U,V,CP)
% Plots the gauss points on a NURBS surface. Additionally, it compares the 
% locations between the real Gauß point locations and the Gauß point 
% positions as they are located in case of assuming non-distorted parameters  

% New Gauß points
% -------------------------------------------------------------------------
der=2;
[GPu,GWu]=get_integrationspunkte(p,U,der);
[GPv,GWv]=get_integrationspunkte(q,V,der);
% GPv=GPu;
k=1; mu = length(U); mv = length(V);
for j = (q+1):(mv-q-1)   %element 1 starts at i=p+1,j=q+1
  for i = (p+1):(mu-p-1)
    % check if element is greater than zero
    if (U(i+1)~=U(i) && V(j+1)~=V(j))
      for kv = 1:GPv(j,1)
        for ku = 1:GPu(i,1)
          % NURBS coordinates u,v from gauss coordinates
          u = GPu(i,ku+1);
          v = GPv(j,kv+1);
          Pt2(k,1:3) = get_point_surf(p,0,u,U,q,0,v,V,CP);
          k = k + 1;
        end
      end    
    end
  end
end

% Original Gauß points 
% -------------------------------------------------------------------------
ngauss = [p+1,q+1];
Pt1 = get_gauss_points(p,q,U,V,CP,ngauss);


% Plot the Gauß points
% -------------------------------------------------------------------------
plot3(Pt1(:,1),Pt1(:,2),Pt1(:,3),'g+');
hold on;
plot3(Pt2(:,1),Pt2(:,2),Pt2(:,3),'rx');
legend(['    Standard Gaußregel: Anzahl =',num2str(length(Pt1))],...
       ['Neue Integrationsregel: Anzahl =',num2str(length(Pt2))],0);

plotNURBS_surf_El_CP(p,q,U,V,CP,'green')
grid on;

title('Lage der Gaußpunkte bei unterschiedlichen Integrationstechniken');

hold off
