function index = plot_referenceAndCurrentConfigurationKirchhoffLoveShellNLinear ...
    (p, q, Xi, Eta, CP, CP_d, homDOFs, technical_parameters, propGraph)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots two windows. The first one shows the deformed geometry with
% loads and supports. The second one shows the strains or stresses on the
% undeformed geometry. The function applies for the visualization of the
% non linear strain resultants.
%
%                Input :
%                  p,q : The polynomial order along the xi- and
%                        eta-parametric coordinates
%               Xi,Eta : The knot vectors along the xi- and eta-parametric 
%                        coordinates
%                   CP : undeformed control points
%                 CP_d : deformed control points
% technical_parameters : structure containing the Young's modulus, the 
%                        thickness of the shell and the Poisson ration
%              homDOFs : The global numbering of the DOFs where
%                        homogeoneous Dirichlet boundary conditions are
%                        applied
%             ropGraph : Structure containing all information on the graphs
%
%               Output :
%                        Graphics
%
% Function layout
%
% 1. First window: Deformed configuration together with supports
%
% 2. Second window: Undeformed configuration together with resultants
% 
%% Function main body

% Read input
mu = length(Xi);
mv = length(Eta);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));
Dm = technical_parameters.E*technical_parameters.t/(1-technical_parameters.nue^2)*[1 technical_parameters.nue 0; technical_parameters.nue 1 0; 0 0 (1-technical_parameters.nue)/2];
Db = technical_parameters.E*technical_parameters.t^3/(12*(1-technical_parameters.nue^2))*[1 technical_parameters.nue 0; technical_parameters.nue 1 0; 0 0 (1-technical_parameters.nue)/2];

%% 1. First window: Deformed configuration together with supports

% On the grid of the graphs
gridu = 49;
gridv = 49;

% deformed geometry, supports and load arrows
[Xp1,Yp1,Zp1] = createBSplineSurface(p,q,Xi,Eta,CP_d,gridu,gridv);
[xs1,ys1,zs1] = createSupports3D(CP,homDOFs);
%[xf1,yf1,zf1] = create_arrows(CP,f);

% Assign the correct numbering to the current graph
figure(propGraph.index)

subplot(2,1,1);

% geometry
surf(Xp1,Yp1,Zp1,'FaceColor','green','EdgeColor','none');
hold;

% element edges
isDeformed = 1;
createEdges(p,q,Xi,Eta,CP_d,isDeformed,gridu,gridv)

% supports
for k =1:length(xs1(:,1))
  plot3(xs1(k,:),ys1(k,:),zs1(k,:),'Linewidth',2,'Color','black');
end
axis equal;

% load arrows
% for k =1:length(xf1(:,1))
%   plot3(xf1(k,:),yf1(k,:),zf1(k,:),'Linewidth',5);
%   plot3(xf1(k,1),yf1(k,1),zf1(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
% end

% control points and polygon
createControlPolygon(CP_d);

camlight left; lighting phong;
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title ('original geometry');
hold off;

%% 2. Second window: Undeformed configuration together with resultants

% Assign a tolerance value
tol=10e-10;

% counting index of lines
l=1;  

%incremental step for v
r=(Eta(mv)-Eta(1))/49;    

% Initial value for the parametric coordinate v
v=Eta(1);

while v <= Eta(mv)+tol
    j = findKnotSpan(v,Eta,nv);
    s = (Xi(mu)-Xi(1))/49;  %incremental step for u
    u = Xi(1);
    k = 1;
  while u <= Xi(mu)+tol
        i = findKnotSpan(u,Xi,nu);
        P(k,l,1:3) = computePointCartesianCoordinatesOnBSplineSurface(p,i,u,Xi,q,j,v,Eta,CP);
        if (propGraph.resultant==1)
            eps_actual = computeVoigtStrainForKirchhoffLoveShellNonLinear(i,p,u,Xi,j,q,v,Eta,CP,CP_d);
            n(k,l,1:3) = Dm*eps_actual(1:3);
            n(k,l,4) = 0.5*(n(k,l,1)+n(k,l,2) + sqrt((n(k,l,1)-n(k,l,2))^2 + 4*n(k,l,3)^2));
            n(k,l,5) = 0.5*(n(k,l,1)+n(k,l,2) - sqrt((n(k,l,1)-n(k,l,2))^2 + 4*n(k,l,3)^2));
        elseif (propGraph.resultant==2)
            kap_actual = computeVoigtCurvatureForKirchhoffLoveShellNonLinear(i,p,u,Xi,j,q,v,Eta,CP,CP_d);
            m(k,l,1:3) = Db*kap_actual(1:3);
            m(k,l,4) = 0.5*(m(k,l,1)+m(k,l,2) + sqrt((m(k,l,1)-m(k,l,2))^2 + 4*m(k,l,3)^2));
            m(k,l,5) = 0.5*(m(k,l,1)+m(k,l,2) - sqrt((m(k,l,1)-m(k,l,2))^2 + 4*m(k,l,3)^2));
        end
        
        % Update counter k
        k=k+1;
        
        % Update u-parametric location
        u=u+s;
  end
  
  % Update counter
  l=l+1;
  
  % Update v-parametric location
  v=v+r;
end

% plot
subplot(2,1,2);
if     (propGraph.resultant==1); surf(P(:,:,1),P(:,:,2),P(:,:,3),n(:,:,propGraph.component));
elseif (propGraph.resultant==2); surf(P(:,:,1),P(:,:,2),P(:,:,3),m(:,:,propGraph.component));   end
shading interp;
colormap('default');

% invert default colormap => red=negativ, blue=positive
COL=colormap;
invCOL(:,1)=COL(:,3);
invCOL(:,2)=COL(:,2);
invCOL(:,3)=COL(:,1);
colormap(invCOL);
% make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);

colorbar;
hold on;
isDeformed = 0;
createEdges(p,q,Xi,Eta,CP,isDeformed,gridu,gridv);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if     (propGraph.resultant==1 && propGraph.component==1)
    title ('normal force n_{11}');
elseif (propGraph.resultant==1 && propGraph.component==2)
    title ('normal force n_{22}');
elseif (propGraph.resultant==1 && propGraph.component==3)
    title ('inplane shear n_{12}');
elseif (propGraph.resultant==1 && propGraph.component==4)
    title ('principal normal force n_1');
elseif (propGraph.resultant==1 && propGraph.component==5)
    title ('principal normal force n_2');
elseif (propGraph.resultant==2 && propGraph.component==1)
    title ('bending moment m_{11}');
elseif (propGraph.resultant==2 && propGraph.component==2)
    title ('bending moment m_{22}');
elseif (propGraph.resultant==2 && propGraph.component==3)
    title ('twisting moment m_{12}');
elseif (propGraph.resultant==2 && propGraph.component==4)
    title ('principal moment m_1');
elseif (propGraph.resultant==2 && propGraph.component==5)
    title ('principal moment m_2');
end

hold off;

% Update the index of the graph
index = propGraph.index + 1;

end
