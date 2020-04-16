function index = plot_referenceConfiguration4IGAIncompressibleFlow2D ... 
    (p, q, Xi, Eta, CP, isNURBS, homDOFs, ihnomDOFs, NBC, ...
    t, int, graph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the configuration at a given time instance for an isogeometric 
% incompressible flow in 2D.
%
%               Input :
%                 p,q : The polynomial degrees of the NURBS patch
%              Xi,Eta : The knot vectors of the NURBS patch
%                  CP : The set of the Control Point coordinates and
%                       weights of the NURBS patch
%             isNURBS : Flag on whether the B-Spline patch is a NURBS or
%                       not
%             homDOFs : The global numbering of the DOFs which belong to
%                       the homogeneous Dirichlet boundary
%           inhomDOFs : The global numbering of the DOFs which belong to
%                       the inhomogeneous Dirichlet boundary
%                 NBC : On the Neumann boundary conditions :
%                        .noConditions : Number of segements where those
%                                        conditions are applied
%                              .xiSpan : .noConditions x 2 array containing
%                                        the knot span extension of the
%                                        segments where those conditions
%                                        are applied in xi-direction
%                             .etaSpan : .noConditions x 2 array containing
%                                        the knot span extension of the
%                                        segments where those conditions
%                                        are applied in eta-direction
%                       .loadAmplitude : .noConditions x 2 array containing 
%                                         the rule for computing the load 
%                                         at each segment
%                       .loadDirection : .noConditions x 2 array containing
%                                        the direction of the load at each 
%                                        segment
%                   t : The time instance for which to draw the
%                       configuration
%                 int : On the spatial integration
%                               .type : 'default' or 'manual'
%                              .xiNGP : No. of GPs along xi-direction for 
%                                       stiffness entries
%                             .etaNGP : No. of GPs along eta-direction for 
%                                       stiffness entries
%                       .xiNGPForLoad : No. of GPs along xi-direction for 
%                                       load entries
%                      .etaNGPForLoad : No. of GPs along eta-direction for 
%                                       load entries
%                         .nGPForLoad : No. of GPs along boundary
%               graph : On the graphics
%                                   .index : Index of the current graph
%                       .postProcComponent : Which component to plot in the
%                                            postprocessing
%              outMsg : On printing information during analysis in the
%                       command window
%
%              Output :
%               index : The index of the current graph
%
% Function layout :
%
% 0. Read input
%
% 1. Create the B-Spline surface
%
% 2. Create the supports over the homogeneous DOFs
%
% 3. Create the supports over the inhomogeneous DOFs
%
% 4. Compute the load vector corresponding to a boundary load and to the current time instance
% ->
%    4i. Get the Neumann boundary extension
%
%   4ii. Get the magnitude of the applied load
%
%  4iii. Get the direction of the applied load
%
%   4iv. Compute the load vector and add it to the existing one
% <-
% 5. Create the load arrows
%
% 6. Plot the B-Spline surface
%
% 7. Plot the knots on the B-Spline surface
%
% 8. Plot the Control Polygon
%
% 9. Plot the supports over the homogeneous DOFs
%
% 10. Plot the supports over the inhomogeneous DOFs
%
% 11. Plot the load arrows corresponding to the boundary load vector of the current time instance
%
% 12. Adjust the graphics settings
%
% 13. Update the graphics index
%
% 14. Update the graphics index
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('________________________________________________________\n');
    fprintf('########################################################\n');
    fprintf('Plotting the reference configuration for an isogeometric \n');
    fprintf('incompressible flow in 2D has been initiated.\n\n');
    fprintf('Time instance t = %d (seconds)\n', t);
    fprintf('________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Dummy variables
prestress = 'undefined';
compPrestress = 'undefined';

% Create a B-Spline patch
BSplinePatch.p = p;
BSplinePatch.q = q;
BSplinePatch.Xi = Xi;
BSplinePatch.Eta = Eta;
BSplinePatch.CP = CP;
BSplinePatch.isNURBS = isNURBS;

% Compute number of Control Points
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));
nCPs = nxi*neta;

% Compute the number of degrees of freedom
nDOFs = 3*nCPs;

% Initialize load vector
F = zeros(nDOFs, 1);

% Get grid points for the visualization in both xi- and eta- direction
xiGrid = 49;
etaGrid = 49;

% Face colort
colorFace = [217 218 219]/255;

% Initialize figure handle
figure(graph.index);

%% 1. Create the B-Spline surface
[X,Y,Z] = createBSplineSurfaceOnCartesianSpace ...
    (p, q, Xi, Eta, CP, isNURBS, prestress, compPrestress, xiGrid, etaGrid);

%% 2. Create the supports over the homogeneous DOFs
[XHom, YHom, ZHom] = createSupportsForIncompressibleFlow2D(CP, homDOFs);

%% 3. Create the supports over the inhomogeneous DOFs
[XInhom, YInhom, ZInhom] = createSupportsForIncompressibleFlow2D(CP, ihnomDOFs);

%% 4. Compute the load vector corresponding to a boundary load and to the current time instance
for iNBC = 1:NBC.noCnd
    %% 4i. Get the Neumann boundary extension
    xib = NBC.xiLoadExtension{iNBC};
    etab = NBC.etaLoadExtension{iNBC};

    %% 4ii. Get the magnitude of the applied load
    loadAmplitude = NBC.loadAmplitude{iNBC};

    %% 4iii. Get the direction of the applied load
    loadDirection = NBC.loadDirection{iNBC};
    
    %% 4iv. Get the flag on whether the load is follower
    isFollower = NBC.isFollower(iNBC, 1);

    %% 4vi. Compute the load vector and add it to the existing one
    funcHandle = str2func(NBC.computeLoadVct{iNBC});
    F =  funcHandle ...
        (F, BSplinePatch, xib, etab, loadAmplitude, loadDirection, ...
        isFollower, t, int, '');
end

%% 5. Create the load arrows
[XLoad, YLoad, ZLoad] = createForceArrowsForIncompressibleFlow2D(CP, F);

%% 6. Plot the B-Spline surface
surf(X, Y, Z, 'FaceColor', colorFace, 'EdgeColor', 'none');
hold on;

%% 7. Plot the knots on the B-Spline surface
plot_knotsForBSplineSurfaceOnCartesianSpace ...
    (p, q, Xi, Eta, CP, isNURBS, 0, xiGrid, etaGrid);

%% 8. Plot the Control Polygon
plot_ControlPolygonBSplineSurface(CP);

%% 9. Plot the supports over the homogeneous DOFs
if norm(homDOFs) ~= 0
    for k = 1:length(XHom(:, 1))
        plot3(XHom(k, :), YHom(k, :), ZHom(k, :), 'Linewidth', 2, 'Color', 'black');
    end
end

%% 10. Plot the supports over the inhomogeneous DOFs
if norm(ihnomDOFs) ~= 0
    for k =1:length(XInhom(:,1))
        plot3(XInhom(k, :), YInhom(k, :), ZInhom(k, :), 'Linewidth', 2, 'Color', 'red');
    end
end

%% 11. Plot the load arrows corresponding to the boundary load vector of the current time instance
for k = 1:length(XLoad(:, 1))
    plot3(XLoad(k, :), YLoad(k, :), ZLoad(k, :), 'Linewidth', 5);
    plot3(xXLoadf(k, 1), YLoad(k, 1), ZLoad(k, 1), 'Marker', 'd', ...
        'MarkerFaceColor', 'blue', 'MarkerSize', 10);
end

%% 12. Adjust the graphics settings
view(2);
axis equal;
camlight left; lighting phong;
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
zlabel('z', 'FontSize', 14);
title(sprintf('Configuration of an isogeometric incompressible flow setting at time t = %d',t));

%% 13. Update the graphics index
index = graph.index + 1;

%% 14. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Plotting reference configuration took %.2d seconds \n\n', computationalTime);
    fprintf('__________Plotting Reference Configuration Ended_________\n');
    fprintf('#########################################################\n\n\n');
end

end