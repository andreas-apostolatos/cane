function index = plot_two_patches_current(p1,q1,U1,V1,CP1,p2,q2,U2,V2,CP2,d1,d2,rb1,rb2,parameters,fl1,fl2,graph)
%% Function documentation
%
% Plots two NURBS patches in the current configuration
%
%         Input : 
%   p1,q1,p2,q2 : polynomial degrees
%   U1,V1,U2,V2 : knot vectors in U,V-directions
%       CP1,CP2 : set of control points and weights
%       rb1,rb2 : Structures containing information on the supports
%    parameters : Technical and geometrical parameters of the plate
%       fl1,fl2 : Force vectors corresponding to each membrane structure
%         graph : structure containing all information on the plots
%
%        Output :
%         index : the index of the graph
%
%% Function main body

% Compute the control point coordinates for the deformed configuration
CPd1 = displaced_control_points(CP1,d1);
CPd2 = displaced_control_points(CP2,d2);

% Number of knots in u,v-direction
mu1 = length(U1);
mv1 = length(V1);

mu2 = length(U2);
mv2 = length(V2);

% Number of Control Points in u,v-direction
nu1 = length(CPd1(:,1,1));
nv1 = length(CPd1(1,:,1));

nu2 = length(CPd2(:,1,1));
nv2 = length(CPd2(1,:,1));

% Material matrix needed to compute the resultants
D = parameters.E/(1-parameters.nu^2)*[1 parameters.nu 0; parameters.nu 1 0; 0 0 (1-parameters.nu)/2];

% Compute point coordinates for the original geometry, the supports and the
% force arrows
if graph.undeformed==1
    % For the patch 1
    [Xp1,Yp1,Zp1] = create_surface(p1,q1,U1,V1,CP1,50,50);
    [xs1,ys1,zs1] = create_supports(CP1,rb1);
    % For the patch 2
    [Xp2,Yp2,Zp2] = create_surface(p2,q2,U2,V2,CP2,50,50);
    [xs2,ys2,zs2] = create_supports(CP2,rb2);
end
if graph.deformed==1
    % For the patch 1
    [Xp_def1,Yp_def1,Zp_def1] = create_surface(p1,q1,U1,V1,CPd1,50,50);
    [xs_def1,ys_def1,zs_def1] = create_supports(CPd1,rb1);
    % For the patch 2
    [Xp_def2,Yp_def2,Zp_def2] = create_surface(p2,q2,U2,V2,CPd2,50,50);
    [xs_def2,ys_def2,zs_def2] = create_supports(CPd2,rb2);
end

% Compute coordinates of the force arrows
% For patch 1
[xf1,yf1,zf1] = create_force_arrows(CP1,fl1);
% For patch 2
[xf2,yf2,zf2] = create_force_arrows(CP2,fl2);

% Number the current plot
figure(graph.index)

% FIRST WINDOW: UNDEFORMED/DEFORMED
subplot(2,1,1);

% graph counter
access = 0;

% Plot the initial geometry
if graph.undeformed==1
    % For the patch 1
    surf(Xp1,Yp1,Zp1,'FaceColor','green','EdgeColor','none');
    if access == 0
        hold;
    end
    % For the patch 2
    surf(Xp2,Yp2,Zp2,'FaceColor','green','EdgeColor','none');
     
    access = access + 1;
    % Plot the element edges for patch 1
    create_edges(p1,q1,U1,V1,CP1,0,50,50)
    % Plot the element edges for patch 2
    create_edges(p2,q2,U2,V2,CP2,0,50,50)
    
    % control points and control polygon for the patch 1
    create_control_polygon(CP1)
    % Plot the supports for the patch 1
    for k = 1:length(xs1(:,1))
        plot3(xs1(k,:),ys1(k,:),zs1(k,:),'Linewidth',2,'Color','black');
    end
    
    % control points and control polygon for the patch 2
    create_control_polygon(CP2)
    % Plot the supports for the patch 2
    for k = 1:length(xs2(:,1))
        plot3(xs2(k,:),ys2(k,:),zs2(k,:),'Linewidth',2,'Color','black');
    end
    
    axis equal;
end

% Plot the deformed geometry
if graph.deformed==1
    if access == 0
        % For the patch 1
        surf(Xp_def1,Yp_def1,Zp_def1,'FaceColor','green','EdgeColor','none');
        hold;
        % For the patch 2
        surf(Xp_def2,Yp_def2,Zp_def2,'FaceColor','green','EdgeColor','none');
        
        % For the patch 1
        create_edges(p1,q1,U1,V1,CPd1,0,50,50)
        create_control_polygon(CPd1)
        % For the patch 2
        create_edges(p2,q2,U2,V2,CPd2,0,50,50)
        create_control_polygon(CPd2)
        
        % Plot the supports for patch 1
        for k = 1:length(xs_def1(:,1))
            plot3(xs_def1(k,:),ys_def1(k,:),zs_def1(k,:),'Linewidth',2,'Color','black');
        end
        
        % Plot the supports for patch 2
        for k = 1:length(xs_def2(:,1))
            plot3(xs_def2(k,:),ys_def2(k,:),zs_def2(k,:),'Linewidth',2,'Color','black');
        end
        
        axis equal;
    else
        % For the patch 1
        create_edges(p1,q1,U1,V1,CPd1,1,50,50)
        % For the patch 2
        create_edges(p2,q2,U2,V2,CPd2,1,50,50)
    end
end

% Plot the load arrows
% For patch 1
for k =1:length(xf1(:,1))
    plot3(xf1(k,:),yf1(k,:),zf1(k,:),'Linewidth',5);
    plot3(xf1(k,1),yf1(k,1),zf1(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end
% For patch 2
for k =1:length(xf2(:,1))
    plot3(xf2(k,:),yf2(k,:),zf2(k,:),'Linewidth',5);
    plot3(xf2(k,1),yf2(k,1),zf2(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

% Adjust plotting parameters
camlight left; 
lighting phong;
view(2);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title ('Reference and/or current geometry');
hold off;

% SECOND WINDOW: STRAIN/STRESS ON UNDEFORMED STRUCTURE

% element displacement vectors for patch 1
d_el1 = zeros(mu1-p1-1,mv1-q1-1,2*(p1+1)*(q1+1));
for j = (q1+1):(mv1-q1-1)
    for i = (p1+1):(mu1-p1-1)
        k=1; 
        for c = j-q1-1:j-1 
            for b = i-p1:i
                d_el1(i,j,k)   = d1(2*(c*nu1+b)-1);
                d_el1(i,j,k+1) = d1(2*(c*nu1+b));
                k=k+2;
            end
        end
    end
end

% element displacement vectors for patch 2
d_el2 = zeros(mu2-p2-1,mv2-q2-1,2*(p2+1)*(q2+1));
for j = (q2+1):(mv2-q2-1)
    for i = (p2+1):(mu2-p2-1)
        k=1; 
        for c = j-q2-1:j-1 
            for b = i-p2:i
                d_el2(i,j,k)   = d2(2*(c*nu2+b)-1);
                d_el2(i,j,k+1) = d2(2*(c*nu2+b));
                k=k+2;
            end
        end
    end
end

% Assign a tolerance value
tol = 10e-10;

% On the stress/strain resultants

% For the 1st patch

% counting index of lines
l=1;  
%incremental step for v
r = (V1(mv1)-V1(1))/49;    
v = V1(1);
while v <= V1(mv1)+tol
    j = findspan(v,V1,nv1);
    s = (U1(mu1)-U1(1))/49;  %incremental step for u
    u = U1(1);
    k = 1;
    while u <= U1(mu1)+tol
        i = findspan(u,U1,nu1);
        P1(k,l,1:3) = point_on_surface(p1,i,u,U1,q1,j,v,V1,CP1);
        d_actual1(:,1) = d_el1(i,j,:);
        
        % Compute the displacement components 
        disp_actual1(k,l,1:2) = get_displacement_vector(i,p1,u,U1,j,q1,v,V1,CP1,d_actual1);
        
        % Compute the total magnitude of the displacement
        if graph.resultant==1 && graph.component==3
            total_disp_actual1(k,l) = (disp_actual1(k,l,1)^2 + disp_actual1(k,l,2)^2)^(1/2);
        end
        
        % Compute the strain components
        eps_actual1 = get_strain_vector(i,p1,u,U1,j,q1,v,V1,CP1,d_actual1);
        
        eps1(k,l,1:3) = eps_actual1(1:3);
        eps1(k,l,4) = 0.5*(eps1(k,l,1)+eps1(k,l,2) + sqrt((eps1(k,l,1)-eps1(k,l,2))^2 + eps1(k,l,3)^2));
        eps1(k,l,5) = 0.5*(eps1(k,l,1)+eps1(k,l,2) - sqrt((eps1(k,l,1)-eps1(k,l,2))^2 + eps1(k,l,3)^2));
        sig1(k,l,1:3) = D*eps_actual1(1:3);
        sig1(k,l,4) = 0.5*(sig1(k,l,1)+sig1(k,l,2) + sqrt((sig1(k,l,1)-sig1(k,l,2))^2 + 4*sig1(k,l,3)^2));
        sig1(k,l,5) = 0.5*(sig1(k,l,1)+sig1(k,l,2) - sqrt((sig1(k,l,1)-sig1(k,l,2))^2 + 4*sig1(k,l,3)^2));
        k=k+1;
        u=u+s;
    end
    l=l+1;
    v=v+r;
end

% For the 1st patch

% counting index of lines
l=1;  
%incremental step for v
r = (V2(mv2)-V2(1))/49;    
v = V2(1);
while v <= V2(mv2)+tol
    j = findspan(v,V2,nv2);
    s = (U2(mu2)-U2(1))/49;  %incremental step for u
    u = U2(1);
    k = 1;
    while u <= U2(mu2)+tol
        i = findspan(u,U2,nu2);
        P2(k,l,1:3) = point_on_surface(p2,i,u,U2,q2,j,v,V2,CP2);
        d_actual2(:,1) = d_el2(i,j,:);
        
        % Compute the displacement components 
        disp_actual2(k,l,1:2) = get_displacement_vector(i,p2,u,U2,j,q2,v,V2,CP2,d_actual2);
        
        % Compute the total magnitude of the displacement
        if graph.resultant==1 && graph.component==3
            total_disp_actual2(k,l) = (disp_actual2(k,l,1)^2 + disp_actual2(k,l,2)^2)^(1/2);
        end
    
        % Compute the strain components
        eps_actual2 = get_strain_vector(i,p2,u,U2,j,q2,v,V2,CP2,d_actual2);
        
        eps2(k,l,1:3) = eps_actual2(1:3);
        eps2(k,l,4) = 0.5*(eps2(k,l,1)+eps2(k,l,2) + sqrt((eps2(k,l,1)-eps2(k,l,2))^2 + eps2(k,l,3)^2));
        eps2(k,l,5) = 0.5*(eps2(k,l,1)+eps2(k,l,2) - sqrt((eps2(k,l,1)-eps2(k,l,2))^2 + eps2(k,l,3)^2));
        sig2(k,l,1:3) = D*eps_actual2(1:3);
        sig2(k,l,4) = 0.5*(sig2(k,l,1)+sig2(k,l,2) + sqrt((sig2(k,l,1)-sig2(k,l,2))^2 + 4*sig2(k,l,3)^2));
        sig2(k,l,5) = 0.5*(sig2(k,l,1)+sig2(k,l,2) - sqrt((sig2(k,l,1)-sig2(k,l,2))^2 + 4*sig2(k,l,3)^2));
        k=k+1;
        u=u+s;
    end
    l=l+1;
    v=v+r;
end

% plot
subplot(2,1,2);
if graph.resultant==1
    if graph.component == 3
        % For the patch 1
        surf(P1(:,:,1),P1(:,:,2),P1(:,:,3),total_disp_actual1(:,:),'EdgeColor','none');
        hold;
        % For the patch 2
        surf(P2(:,:,1),P2(:,:,2),P2(:,:,3),total_disp_actual2(:,:),'EdgeColor','none');
    else
        % For the patch 1
        surf(P1(:,:,1),P1(:,:,2),P1(:,:,3),disp_actual1(:,:,graph.component),'EdgeColor','none');
        hold;
        % For the patch 2
        surf(P2(:,:,1),P2(:,:,2),P2(:,:,3),disp_actual2(:,:,graph.component),'EdgeColor','none');
    end   
elseif graph.resultant==2
    % For the patch 1
    surf(P1(:,:,1),P1(:,:,2),P1(:,:,3),eps1(:,:,graph.component),'EdgeColor','none');
    hold;
    % For the patch 2
    surf(P2(:,:,1),P2(:,:,2),P2(:,:,3),eps2(:,:,graph.component),'EdgeColor','none');
elseif graph.resultant==3
    % For the patch 1
    surf(P1(:,:,1),P1(:,:,2),P1(:,:,3),sig1(:,:,graph.component),'EdgeColor','none'); 
    hold;
    % For the patch 2
    surf(P2(:,:,1),P2(:,:,2),P2(:,:,3),sig2(:,:,graph.component),'EdgeColor','none'); 
end

% Graphics options
shading interp;
colormap('default');
view(2);
axis equal;
grid off;

% invert default colormap => red=negativ, blue=positive
% COL=colormap;
% invCOL(:,1)=COL(:,3);
% invCOL(:,2)=COL(:,2);
% invCOL(:,3)=COL(:,1);
% colormap(invCOL);
% 
% % make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);

colorbar;
hold on;

% For the patch 1
% create_edges(p1,q1,U1,V1,CP1,0,50,50);
% For the patch 2
% create_edges(p2,q2,U2,V2,CP2,0,50,50);

view(2);
% axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if graph.resultant==1 && graph.component==1 
    title('displacement component u_x');
elseif graph.resultant==1 && graph.component==2 
    title('displacement component u_x');
elseif graph.resultant==1 && graph.component==3 
    title('displacement magnitude u=(u_x^2+u_y^2)^0.5');
elseif graph.resultant==2 && graph.component==1 
    title('strain epsilon xx');
elseif graph.resultant==2 && graph.component==2 
    title('strain epsilon yy');
elseif graph.resultant==2 && graph.component==3 
    title('strain epsilon xy');
elseif graph.resultant==2 && graph.component==4 
    title('principal strain e1');
elseif graph.resultant==2 && graph.component==5 
    title('principal strain e2');
elseif graph.resultant==3 && graph.component==1
    title('stress sigma xx');
elseif graph.resultant==3 && graph.component==2
    title('stress sigma yy');
elseif graph.resultant==3 && graph.component==3 
    title('stress sigma xy');
elseif graph.resultant==3 && graph.component==4
    title('principal stress s1');
elseif graph.resultant==3 && graph.component==5 
    title('principal stress s2');
end

hold off;

% Update plot index by 1
index = graph.index+1;
end

