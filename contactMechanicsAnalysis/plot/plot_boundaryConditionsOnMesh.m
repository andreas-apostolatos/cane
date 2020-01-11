function plot_boundaryConditionsOnMesh(mesh,homDBC,F)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the boundary conditions (Dirichlet/Neumann) onto a given mesh
%
%    Input :
%     mesh : Nodes and elements of the mesh
%   homDBC : Vector containing the the DOFs (via their global numbering)
%            which are prescribed 
%        F : Global load vector
%
%   Output :
%            Graphics
%
%% Function main body

% Create the supports
[xs,ys,zs] = createSupports(mesh.nodes,homDBC);

%supports
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end

% Hold on the graphics
hold on;

% Create the force arrows
[xf,yf,zf] = createForceArrows(mesh.nodes,F);

% load arrows  
for k =1:length(xf(:,1))
    plot3(xf(k,:),yf(k,:),zf(k,:),'Color','blue','Linewidth',5);
    plot3(xf(k,1),yf(k,1),zf(k,1),'Color','blue','Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

end

