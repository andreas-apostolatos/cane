function  plot_postprocCurrentConfigurationIGAThinStructure ...
    (p, q, Xi, Eta, CP, isNURBS, xiGrid, etaGrid, homDOFs, Fl, dHat, graph)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots a window with the reference and/or the current configuration of an
% IGA Kirchhoff-Love shell.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  shell
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
% xiGrid,etaGrid : The grid points used for the plotting of the NURBS
%                  geometry
%        homDOFs : Vector containing information on the supports
%             Fl : The applied load vector
%           dHat : The displacement field of the Control Points
%          graph : Information on the graphics
%
%         Output :
%                  graphics
%
% Function layout :
%
% 1. Create the arrays containing the Cartesian coordinates of the B-Spline surface grid points, the supports and the arrow vectors
%
% 2. Plot the B-Spline surfaces
%
% 3. Plot the knots on the B-Spline surfaces
%
% 4. Plot the control polygon of the B-Spline surfaces
%
% 5. Plot the supports on the geometries
%
% 6. Plot the load arrows on the geometries
%
%% Function main body

%% 0. Read input

% Initialize flags
isUndeformed = 0;
isDeformed = 0;

% Define grey color
color = [217 218 219]/255; %'white'; % 'gray'
% color = 'none';

% Define dummy arrays
prestress = 'undefined';
compPrestress = 'undefined';

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsForIGAKirchhoffLoveShell(CP,dHat);

%% 1. Create the arrays containing the Cartesian coordinates of the B-Spline surface grid points, the supports and the arrow vectors
if strcmp(graph.postprocConfig,'reference')||strcmp(graph.postprocConfig,'referenceCurrent')
    % Create the coordinates of the grid points on the surface
    [XpRef,YpRef,ZpRef] = createBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,prestress,compPrestress,xiGrid,etaGrid);

    % Create the coordinates of the vertices of the support triangles
    [xsRef,ysRef,zsRef] = createSupports3D(CP,homDOFs);

    % Create the start and end points of the arrows representing the loads
    [xfRef,yfRef,zfRef] = createForceArrows3D(CP,Fl);
end
if strcmp(graph.postprocConfig,'current')||strcmp(graph.postprocConfig,'referenceCurrent')
    % Create the coordinates of the grid points on the surface
    [XpCur,YpCur,ZpCur] = createBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CPd,isNURBS,prestress,compPrestress,xiGrid,etaGrid);
    
    % Create the coordinates of the vertices of the support triangles
    [xsCur,ysCur,zsCur] = createSupports3D(CPd,homDOFs);

    % Create the start and end points of the arrows representing the loads
    [xfCur,yfCur,zfCur] = createForceArrows3D(CPd,Fl);
end

%% 2. Plot the B-Spline surfaces
if strcmp(graph.postprocConfig,'reference')
    surf(XpRef, YpRef, ZpRef, FaceColor=color, EdgeColor='none', FaceAlpha=0.3);
    hold on;
elseif strcmp(graph.postprocConfig,'current')
    surf(XpCur, YpCur, ZpCur, FaceColor=color, EdgeColor='none');
    hold on;
elseif strcmp(graph.postprocConfig, 'referenceCurrent')
    surf(XpRef, YpRef, ZpRef, FaceColor='none', EdgeColor='none');
    hold on;
    surf(XpCur, YpCur, ZpCur, FaceColor=color, EdgeColor='none', FaceAlpha=0.3);
end

%% 3. Plot the knots on the B-Spline surfaces
if strcmp(graph.postprocConfig,'reference')
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isUndeformed,xiGrid,etaGrid);
elseif strcmp(graph.postprocConfig,'current')
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CPd,isNURBS,isUndeformed,xiGrid,etaGrid);
elseif strcmp(graph.postprocConfig,'referenceCurrent')
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isUndeformed,xiGrid,etaGrid);
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CPd,isNURBS,isDeformed,xiGrid,etaGrid);
end

%% 4. Plot the control polygon of the B-Spline surfaces
% if strcmp(graph.postprocConfig,'reference')
%     plot_ControlPolygonBSplineSurface(CP);
% elseif strcmp(graph.postprocConfig,'current')
%     plot_ControlPolygonBSplineSurface(CPd);
% elseif strcmp(graph.postprocConfig,'referenceCurrent')
% %     plot_ControlPolygonBSplineSurface(CP);
%     plot_ControlPolygonBSplineSurface(CPd);
% end

%% 5. Plot the supports on the geometries
% if strcmp(graph.postprocConfig,'reference')||strcmp(graph.postprocConfig,'referenceCurrent')
%     for xiCounter = 1:length(xsRef(:,1))
%         plot3(xsRef(xiCounter,:),ysRef(xiCounter,:),zsRef(xiCounter,:),'Linewidth',2,'Color','black');
%     end
% end
% if strcmp(graph.postprocConfig,'current')||strcmp(graph.postprocConfig,'referenceCurrent')
%     for xiCounter = 1:length(xsCur(:,1))
%         plot3(xsCur(xiCounter,:),ysCur(xiCounter,:),zsCur(xiCounter,:),'Linewidth',2,'Color','black');
%     end
% end

%% 6. Plot the load arrows on the geometries
% if strcmp(graph.postprocConfig,'reference')||strcmp(graph.postprocConfig,'referenceCurrent')
%     for xiCounter =1:length(xfRef(:,1))
%         plot3(xfRef(xiCounter,:),yfRef(xiCounter,:),zfRef(xiCounter,:),'Linewidth',5);
%         plot3(xfRef(xiCounter,1),yfRef(xiCounter,1),zfRef(xiCounter,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
%     end
% end
% if strcmp(graph.postprocConfig,'current')||strcmp(graph.postprocConfig,'referenceCurrent')
%     for xiCounter =1:length(xfCur(:,1))
%         plot3(xfCur(xiCounter,:),yfCur(xiCounter,:),zfCur(xiCounter,:),'Linewidth',5);
%         plot3(xfCur(xiCounter,1),yfCur(xiCounter,1),zfCur(xiCounter,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
%     end
% end

end
