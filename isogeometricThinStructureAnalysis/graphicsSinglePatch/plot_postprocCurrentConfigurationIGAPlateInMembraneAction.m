function  plot_postprocCurrentConfigurationIGAPlateInMembraneAction ...
    (p, q, Xi, Eta, CP, isNURBS, xiGrid, etaGrid, homDOFs, Fl, dHat, ...
    propGraph)
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
% IGA plate in membrane action.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
% xiGrid,etaGrid : The grid points used for the plotting of the NURBS
%                  geometry
%        homDOFs : Array of the global numbering of the DOFs where
%                  homogeneous Dirichlet boundary conditions are applied
%             Fl : The applied load vector
%           dHat : The displacement field of the control points
%      propGraph : Structure containing information on the figures,
%                      .index : Index of the current figure
%
%         Output :
%                  graphics
%
% Function layout :
%
% 1. Create the arrays containing the Cartesian coordinates of the B-Spline surface grid points
%
% 2. Create the arrays containing the Cartesian coordinates of the supports and the arrow vectors
%
% 3. Plot the B-Spline surfaces
%
% 4. Plot the knots on the B-Spline surfaces
%
% 5. Plot the control polygon of the B-Spline surfaces
%
% 6. Plot the supports on the geometries
%
% 7. Plot the load arrows on the geometries
%
%% Function main body

%% 0. Read input

% Initialize flags
isUndeformed = 0;
isDeformed = 0;

% Initialize dummy variables
prestress = 'undefined';
compPrestress = 'undefined';

% Face color
FaceColor = [217 218 219]/255;

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsForIGAPlateInMembraneAction(CP, dHat);

%% 1. Create the arrays containing the Cartesian coordinates of the B-Spline surface grid points
if strcmp(propGraph.postprocConfig, 'current') || ...
        strcmp(propGraph.postprocConfig, 'referenceCurrent')
    [XpCur,YpCur,ZpCur] = createBSplineSurfaceOnCartesianSpace...
        (p,q,Xi,Eta,CPd,isNURBS,prestress,compPrestress,xiGrid,etaGrid);
end

%% 2. Create the arrays containing the Cartesian coordinates of the supports and the arrow vectors
if strcmp(propGraph.postprocConfig, 'reference') || ...
        strcmp(propGraph.postprocConfig, 'referenceCurrent')
    % Create the coordinates of the grid points on the surface
    [XpRef, YpRef, ZpRef] = createBSplineSurfaceOnCartesianSpace ...
        (p, q, Xi, Eta, CP, isNURBS, prestress, compPrestress, xiGrid, etaGrid);

    % Create the coordinates of the vertices of the support triangles
    [xsRef, ysRef, zsRef] = createSupports2D(CP, homDOFs);

    % Create the start and end points of the arrows representing the loads
    [xfRef, yfRef, zfRef] = createForceArrows2D(CP, Fl);
end

%% 3. Plot the B-Spline surfaces
if strcmp(propGraph.postprocConfig, 'reference')
    surf(XpRef, YpRef, ZpRef, 'FaceColor', FaceColor, 'EdgeColor', 'none');
    hold on;
elseif strcmp(propGraph.postprocConfig, 'current')
    surf(XpCur, YpCur, ZpCur, 'FaceColor', FaceColor, 'EdgeColor', 'none');
    hold on;
elseif strcmp(propGraph.postprocConfig, 'referenceCurrent')
    surf(XpRef, YpRef, ZpRef, 'FaceColor', 'none', 'EdgeColor', 'none');
    hold on;
    surf(XpCur, YpCur, ZpCur, 'FaceColor', FaceColor, 'EdgeColor', 'none');
end

%% 3. Plot the knots on the B-Spline surfaces
if strcmp(propGraph.postprocConfig, 'reference')
    plot_knotsForBSplineSurfaceOnCartesianSpace ...
        (p, q, Xi, Eta, CP, isNURBS, isUndeformed, xiGrid, etaGrid);
elseif strcmp(propGraph.postprocConfig, 'current')
    plot_knotsForBSplineSurfaceOnCartesianSpace ...
        (p, q, Xi, Eta, CPd, isNURBS, isUndeformed, xiGrid, etaGrid);
elseif strcmp(propGraph.postprocConfig, 'referenceCurrent')
    plot_knotsForBSplineSurfaceOnCartesianSpace ...
        (p, q, Xi, Eta, CP, isNURBS, isUndeformed, xiGrid, etaGrid);
    plot_knotsForBSplineSurfaceOnCartesianSpace ...
        (p, q, Xi, Eta, CPd, isNURBS, isDeformed, xiGrid, etaGrid);
end

%% 4. Plot the control polygon of the B-Spline surfaces
if strcmp(propGraph.postprocConfig, 'reference')
%     plot_ControlPolygonBSplineSurface(CP);
elseif strcmp(propGraph.postprocConfig, 'current')
    plot_ControlPolygonBSplineSurface(CPd);
elseif strcmp(propGraph.postprocConfig, 'referenceCurrent')
%     plot_ControlPolygonBSplineSurface(CP);
    plot_ControlPolygonBSplineSurface(CPd);
end

%% 5. Plot the supports on the geometries
if strcmp(propGraph.postprocConfig, 'reference') || ...
        strcmp(propGraph.postprocConfig, 'referenceCurrent')
    for iXi = 1:length(xsRef(:, 1))
        plot3(xsRef(iXi, :), ysRef(iXi, :), zsRef(iXi, :), ...
            'Linewidth', 2, 'Color', 'black');
    end
end

%% 6. Plot the load arrows on the geometries
if strcmp(propGraph.postprocConfig, 'reference') || ...
        strcmp(propGraph.postprocConfig, 'referenceCurrent')
    for iXi = 1:length(xfRef(:, 1))
        plot3(xfRef(iXi, :), yfRef(iXi, :), zfRef(iXi, :), 'Color', 'blue', ...
            'Linewidth', 5);
        plot3(xfRef(iXi, 1), yfRef(iXi, 1), zfRef(iXi, 1), 'Marker', 'd', ...
            'MarkerFaceColor', 'blue', 'MarkerSize', 10);
    end
end

end