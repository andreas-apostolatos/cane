function plot_referenceConfigurationIGAThinStructure ...
    (p, q, Xi, Eta, CP, isNURBS, homDOFs, Fl, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function
%
% Plots the reference configuration for an isogeometric Kirchhoff-Love
% shell which includes the support triangles as well as the load arrows
%
%   Input : 
%     p,q : The polynomial degrees in xi-,eta-directions
%  Xi,Eta : The knot vectors in xi-,eta- directions
%      CP : The set of the Control points and weights
% isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
% homDOFs : The set of boundary conditions
%      Fl : The force vector
%  outMsg : Whether or not to output message on refinement progress
%           'outputEnabled' : enables output information
%
%  Output : 
%           graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Create the arrays containing the Cartesian coordinates of the B-Spline surface grid points, the supports and the arrow vectors
%
% 2. Plot the B-Spline surface
%
% 3. Plot the knots on the B-Spline surface
%
% 4. Plot the control polygon of the B-Spline surface
%
% 5. Plot the supports on the geometry
%
% 6. Plot the load arrows on the geometry
%
% 7. Assign plotting properties
%
% 8. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_______________________________________________________\n');
    fprintf('#######################################################\n');
    fprintf('Plotting the reference configuration for a single patch\n');
    fprintf('isogeometric Kirchhoff-Love shell has been initiated\n');
    fprintf('_______________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Grid point number for the plotting of both the B-Spline surface and the
% knots
xiGrid = 49;
etaGrid = 49;

prestress = 'undefined';
compPrestress = 'undefined';

% Grey color definition
grey = [217 218 219]/255;
% grey = 'none';

%% 1. Create the arrays containing the Cartesian coordinates of the B-Spline surface grid points, the supports and the arrow vectors

% Create the coordinates of the grid points on the surface
[Xp,Yp,Zp] = createBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,prestress,compPrestress,xiGrid,etaGrid);

% Create the coordinates of the vertices of the support triangles
[xs,ys,zs] = createSupports3D(CP,homDOFs);

% Create the start and end points of the arrows representing the loads
[xf,yf,zf] = createForceArrows3D(CP,Fl);

%% 2. Plot the B-Spline surface
surf(Xp,Yp,Zp,'FaceColor',grey,'EdgeColor','none');
hold on;
 
%% 3. Plot the knots on the B-Spline surface
isDeformed = 0;
plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isDeformed,xiGrid,etaGrid);
  
%% 4. Plot the control polygon of the B-Spline surface
plot_ControlPolygonBSplineSurface(CP);

%% 5. Plot the supports on the geometry
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end
  
%% 6. Plot the load arrows on the geometry
for k =1:length(xf(:,1))
    plot3(xf(k,:),yf(k,:),zf(k,:),'Color','blue','Linewidth',5);
    plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end
  
%% 7. Assign plotting properties
axis equal;
camlight(160,50);
camlight left;
lighting phong;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);

%% 8. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Plotting the reference configuration took %.2d seconds \n\n',computationalTime);
    fprintf('_________Plotting Reference Configuration Ended_________\n');
    fprintf('########################################################\n\n\n');
end

end
