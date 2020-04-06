function index = plot_postprocIGAIncompressibleFlow2D ...
    (BSplinePatch, up, homDOFs, inhomDOFs, F, graph, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the geometry together with the visualization of the selected
% resultant along the B-Spline patch for an isogeometric incompressible 
% flow problem in 2D.
%
%        Input :
% BSplinePatch : The polynomial degrees, the knot vectors and the Control 
%                Point coordinates/weights for the given B-Spline patch
%           up : The solution vector for the current time instance
%      homDOFs : The global numbering of the DOFs which belong to the 
%                homogeneous Dirichlet boundary
%    inhomDOFs : The global numbering of the DOFs which belong to the
%                inhomogeneous Dirichlet boundary
%            F : The boundary load vector of the current time instance
%        graph : On the graphics
%                             .index : Index of the current graph
%                 .postProcComponent : Which component to plot in the
%                                      postprocessing
%       outMsg : Enables outputting information in the command window when
%                selected as 'outputEnabled'
%
% Function layout :
%
% 0. Read input
%
% 1. Create the supports over the homogeneous DOFs
%
% 2. Create the supports over the inhomogeneous DOFs
%
% 3. Create the load arrows
%
% 4. Distribute the solution vector into the element solution vectors
%
% 5. Loop over all the parametric locations in the eta-direction
% ->
%    5i. Initialize the coordinate in the xi-direction
%
%   5ii. Find the span in eta-direction
%
%  5iii. Reset the counter in the xi-direction
%
%   5iv. Loop over all the parametric locations in the xi-direction
%   ->
%        5iv.1. Find the span in xi-direction
%
%        5iv.2. Compute the NURBS basis functions at the parametric location
%
%        5iv.3. Compute the Cartesian coordinates of the current parametric location on the reference NURBS patch
%
%        5iv.4. Get the element nodal vector [ux uy p]'
%
%        5iv.5. Compute the actual transport vector at the parametric location
%
%        5iv.6 Update the counter and the parametic location in xi-direction
%   <-
%    5v. Update the counter and the parametic location in eta-direction
% <-
% 
% 6. Plot the B-Spline surface together with the visualization resultant
%
% 7. Plot the supports over the homogeneous DOFs
%
% 8.  Plot the supports over the inhomogeneous DOFs
%
% 9. Plot the load arrows corresponding to the boundary load vector of the current time instance
%
% 10. Plot the knots over the B-Spline surface
%
% 11. Assign a title to the figure
%
% 12. Adjust figure properties
%
% 13. Update the graphics index
%
% 14. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________\n');
    fprintf('###########################################################\n');
    fprintf('Plotting the postprocessing resultants for an isogeometric \n');
    fprintf('incompressible flow in 2D has been initiated.\n\n');
    fprintf('___________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Get the parameters of the B-Spline patch
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Number of knots in xi,eta-direction
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);

% Number of Control Points in xi,eta-direction
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Choose the number of the grid points
xiGrid = 99;
etaGrid = 99;

% Initialize counter along the eta-parametric direction
counter_eta = 1;  

% Define incremental step for xi
if ~strcmp(graph.postProcComponent, 'velocityVectorPlot')
    dxi = (Xi(numKnots_xi) - Xi(1))/xiGrid;
else
    dxi = 2*(Xi(numKnots_xi) - Xi(1))/xiGrid;
end

% Define incremental step for eta
if ~strcmp(graph.postProcComponent, 'velocityVectorPlot')
    deta = (Eta(numKnots_eta) - Eta(1))/etaGrid;
else
    deta = 2*(Eta(numKnots_eta) - Eta(1))/etaGrid;
end

% Get a starting coordinate on the eta-parameter line
eta = Eta(1);

% Assign a tolerance value
tol = 1e-9;

% Initialize plotting arrays
P = zeros(xiGrid + 1, etaGrid + 1, 3);
if ~strcmp(graph.postProcComponent, '2normVelocity')
    upVector = zeros(xiGrid + 1, etaGrid + 1, 3);
else
    upVector = zeros(xiGrid + 1, etaGrid + 1);
end

% Initialize figure handle
figure(graph.index);

%% 1. Create the supports over the homogeneous DOFs
[XHom, YHom, ZHom] = createSupportsForIncompressibleFlow2D ...
    (CP, homDOFs);

%% 2. Create the supports over the inhomogeneous DOFs
[XInhom, YInhom, ZInhom] = createSupportsForIncompressibleFlow2D ...
    (CP, inhomDOFs);

%% 3. Create the load arrows
[XLoad, YLoad, ZLoad] = createForceArrowsForIncompressibleFlow2D ...
    (CP, F);

%% 4. Distribute the solution vector into the element solution vectors
upEl = zeros(numKnots_xi - p - 1, numKnots_eta - q - 1, 3*(p + 1)*(q + 1));
for etaSpan = (q + 1):(numKnots_eta - q - 1)
    for xiSpan = (p + 1):(numKnots_xi - p - 1)
        counter_xi = 1;
        for c = etaSpan - q - 1:etaSpan - 1 
            for b = xiSpan - p:xiSpan
                upEl(xiSpan, etaSpan, counter_xi)   = up(3*(c*numCPs_xi + b) - 2);
                upEl(xiSpan, etaSpan, counter_xi + 1) = up(3*(c*numCPs_xi + b) - 1);
                upEl(xiSpan, etaSpan, counter_xi+2) = up(3*(c*numCPs_xi + b));
                counter_xi = counter_xi + 3;
            end
        end
    end
end

%% 5. Loop over all the parametric locations in the eta-direction
while eta <= Eta(numKnots_eta) + tol
    %% 5i. Initialize the coordinate in the xi-direction
    xi = Xi(1);
    
    %% 5ii. Find the span in eta-direction
    etaSpan = findKnotSpan(eta, Eta, numCPs_eta);
    
    %% 5iii. Reset the counter in the xi-direction
    counter_xi = 1;
    
    %% 5iv. Loop over all the parametric locations in the xi-direction
    while xi <= Xi(numKnots_xi) + tol
        %% 5iv.1. Find the span in xi-direction
        xiSpan = findKnotSpan ...
            (xi, Xi, numCPs_xi);
        
        %% 5iv.2. Compute the NURBS basis functions at the parametric location
        RMtx = computeIGABasisFunctionsAndDerivativesForSurface ...
            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, 0);
        
        %% 5iv.3. Compute the Cartesian coordinates of the current parametric location on the reference NURBS patch
        P(counter_xi, counter_eta, :) = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, RMtx);
        
        %% 5iv.4. Get the element nodal vector [ux uy p]'
        upActual(:, 1) = upEl(xiSpan, etaSpan, :);
        
        %% 5iv.5. Compute the actual transport vector at the parametric location
        up = computeNodalVectorIncompressibleFlow2D ...
            (RMtx, p, q, upActual);
        if ~strcmp(graph.postProcComponent, '2normVelocity')
            upVector(counter_xi, counter_eta, :) = up;
        else
            upVector(counter_xi, counter_eta) = norm(up(1:2, 1));
        end
        
        %% 5iv.6 Update the counter and the parametic location in xi-direction
        counter_xi = counter_xi + 1;
        xi = xi + dxi;
    end
    
    %% 5v. Update the counter and the parametic location in eta-direction
    counter_eta = counter_eta + 1;
    eta = eta + deta;
end

%% 6. Plot the B-Spline surface together with the visualization resultant
if strcmp(graph.postProcComponent,'xVelocity') || ...
        strcmp(graph.postProcComponent,'yVelocity') || ...
        strcmp(graph.postProcComponent,'pressure')
    if strcmp(graph.postProcComponent,'xVelocity')
        surf(P(:, :, 1), P(:, :, 2), P(:, :, 3), upVector(:, :, 1));
    elseif strcmp(graph.postProcComponent, 'yVelocity')
        surf(P(:, :, 1), P(:, :, 2), P(:, :, 3), upVector(:, :, 2));
    elseif strcmp(graph.postProcComponent, 'pressure')
        surf(P(:, :, 1), P(:, :, 2), P(:, :, 3), upVector(:, :, 3));
    end
elseif strcmp(graph.postProcComponent, '2normVelocity')
    surf(P(:, :, 1), P(:, :, 2), P(:, :, 3), upVector);
elseif strcmp(graph.postProcComponent, 'velocityVectorPlot')
    scale = 4;
%     quiver(P(:, :, 1), P(:, :, 2), upVector(:, :, 1), upVector(:, :, 2), 'autoscale', 'off');
    quiver(P(:, :, 1), P(:, :, 2), upVector(:, :, 1), upVector(:, :, 2), scale);
elseif strcmp(graph.postProcComponent, 'velocityContourPlot')
    % Compute the mesh grid
%     [sx, sy] = meshgrid(P(1, 1, 1), min(min(P(:, :, 2))):.1:max(max(P(:, :, 2))));
    [sx, sy] = meshgrid(min(min(P(:, :, 1))):.01:max(max(P(:, :, 1))), ...
        min(min(P(:, :, 2))):.01:max(max(P(:, :, 2))));
    streamline(P(:, :, 1), P(:, :, 2), -P(:, :, 2), P(:, :, 1), sx, sy);
    
    % Visualize the streamline
%     streamline(stream2(P(:, :, 1), P(:, :, 2), uVector(:, :, 1), uVector(:, :, 2), sx, sy));
    streamline(stream2(P(:, :, 1), P(:, :, 2), -P(:, :, 2), P(:, :, 1), sx, sy));
end
hold on;

%% 7. Plot the supports over the homogeneous DOFs
if norm(homDOFs) ~= 0
    for iHom = 1:length(XHom(:, 1))
        plot3(XHom(iHom, :), YHom(iHom, :), ZHom(iHom, :), 'Linewidth', ...
            2, 'Color', 'black');
    end
end

%% 8.  Plot the supports over the inhomogeneous DOFs
if norm(inhomDOFs) ~= 0
    for iHom = 1:length(XInhom(:, 1))
        plot3(XInhom(iHom, :), YInhom(iHom, :), ZInhom(iHom, :), ...
            'Linewidth', 2, 'Color', 'red');
    end
end

%% 9. Plot the load arrows corresponding to the boundary load vector of the current time instance
for iHom = 1:length(XLoad(:, 1))
    plot3(XLoad(iHom, :), YLoad(iHom, :), ZLoad(iHom, :), 'Linewidth', 5);
    plot3(xXLoadf(iHom, 1), YLoad(iHom, 1), ZLoad(iHom, 1), 'Marker', 'd', ...
        'MarkerFaceColor', 'blue', 'MarkerSize', 10);
end

%% 10. Plot the knots over the B-Spline surface
if ~strcmp(graph.postProcComponent, 'velocityVectorPlot') || ...
        ~strcmp(graph.postProcComponent, 'velocityContourPlot')
    plot_knotsForBSplineSurfaceOnCartesianSpace ...
        (p, q, Xi, Eta, CP, isNURBS, 0, xiGrid, etaGrid);
end

%% 11. Assign a title to the figure
if strcmp(graph.postProcComponent, 'xVelocity')
    title(sprintf('Velocity component u_x'));
elseif strcmp(graph.postProcComponent, 'yVelocity')
    title(sprintf('Velocity component u_y'));
elseif strcmp(graph.postProcComponent, 'pressure')
    title(sprintf('Pressure distribution'));
elseif strcmp(graph.postProcComponent, '2normVelocity')
    title(sprintf('Velocity magnitude ||u||_2'));
elseif strcmp(graph.postProcComponent, 'velocityVectorPlot')
    title(sprintf('Vector plot for the velocity field'));
end

%% 12. Adjust figure properties
if ~strcmp(graph.postProcComponent,'velocityVectorPlot') || ...
        ~strcmp(graph.postProcComponent,'velocityContourPlot')
    % Graphics options
    shading interp;
    colormap('default');
    
    % On the color bar
    colorbar;
end
view(2);
axis equal;
xlabel('x','FontSize', 14);
ylabel('y','FontSize', 14);

%% 13. Update the graphics index
index = graph.index + 1;

%% 14. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Plotting of the current configuration took %.2d seconds \n\n', computationalTime);
    fprintf('__Plotting Postprocessing Resultants Configuration Ended__\n');
    fprintf('##########################################################\n\n\n');
end

end