function plot_knotsForBSplineSurfaceOnCartesianSpace...
    (p, q, Xi, Eta, CP, isNURBS, isDeformed, xiGrid, etaGrid)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Draws the element edges for the NURBS surface, i.e the knots on the
% geometry
%
%      Input :
%        p,q : Polynomial degrees
%     Xi,Eta : Knot vectors in u,v-direction
%         CP : Control Point coordinates and weights
%    isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
% isDeformed : Flag on whether the plot is called for the reference or the
%              current configuration
%     xiGrid : Points to use in xi-direction
%    etaGrid : Points to use in eta-direction
%
%     Output : 
%              Graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Plot the edges in xi-direction
%
% 2. Plot the edges in eta-direction
%
% 3. Plot all the element edges
%
%% Function main body

%% 0. Read input

% Number of knots in u,v-direction
mxi = length(Xi);
meta = length(Eta);

% Number of Control Points in u,v-direction
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));

% Make the knot vectors unique
XiUnique = unique(Xi);
EtaUnique = unique(Eta);

% Color of the edge
colorEdge = 'black';
% colorEdge = 'red';
% colorEdge = 'none';

% Assign a tolerance value
eps = 1e-9;

% Initialize counter
l = 1;

% Number of derivatives for the B-Spline basis functions
numDrv = 0;

% Compute step size in xi-direction
dxi = (Xi(mxi) - Xi(1))/xiGrid;

% Compute step size in eta-direction
deta = (Eta(meta) - Eta(1))/etaGrid; 

% Initialize plotting array
P = zeros(xiGrid, etaGrid, 3);

%% 1. Plot the edges in xi-direction
for j2 = 1:length(EtaUnique)
    % Get the starting coordinate in eta-direction
    eta = EtaUnique(j2);
    
    % Find the span in eta-direction
    etaSpan = findKnotSpan(eta, Eta, neta);
    
    % Get the starting coordinate in xi-direction
    xi = XiUnique(1);
    
    % Initialize counter for the edges in xi-direction
    k = 1;
    
    % Loop over all the coordinates in xi-direction
    while xi <= XiUnique(end) + eps
        % Find the span in xi-direction
        xiSpan = findKnotSpan(xi, Xi, nxi);
        
        % Compute the IGA basis functions in xi-direction
        R = computeIGABasisFunctionsAndDerivativesForSurface ...
            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, numDrv);
        
        % Compute the Cartesian image of the parametric point (xi,eta)
        P(k, l, 1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, R);
        
        % Update counter for the edges in xi-direction
        k = k + 1;
        
        % Update the parametric coordinate in xi-direction
        xi = xi + dxi;
    end
    % Update the counter for the edges in eta-direction
    l = l + 1;
end

%% 2. Plot the edges in eta-direction
for i2 = 1:length(XiUnique)
    % Get the starting coordinate in xi-direction
    xi = XiUnique(i2);
    
    % Find the span in xi-direction
    xiSpan = findKnotSpan(xi, Xi, nxi);
    
    % Get the starting coordinate in eta-direction
    eta = EtaUnique(1);
    
    % Initialize counter for the edges in eta-direction
    k = 1;
    
    % Loop over all the coordinates in eta-direction
    while eta <= EtaUnique(end) + eps
        % Find the span in eta-direction
        etaSpan = findKnotSpan(eta, Eta, neta);
        
        % Compute the IGA basis functions in xi-direction
        R = computeIGABasisFunctionsAndDerivativesForSurface ...
            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, isNURBS, numDrv);
        
        % Compute the Cartesian image of the parametric point (xi,eta)
        P(k,l,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface ...
            (xiSpan, p, xi, Xi, etaSpan, q, eta, Eta, CP, R);
        
        % Update counter for the edges in eta-direction
        k = k + 1;
        
        % Update the parametric coordinate in eta-direction
        eta = eta + deta;
    end
    % Update the counter for the edges in xi-direction
    l = l + 1;
end

%% 3. Plot all the element edges
if ~isDeformed
    plot3(P(:, :, 1), P(:, :, 2) ,P(:, :, 3), 'Color', colorEdge, 'LineWidth', .01);
else
    plot3(P(:, :, 1), P(:, :, 2), P(:, :, 3), 'Color', colorEdge, 'LineStyle', '-.');
end
% axis equal;
grid on;
xlabel('x', 'FontSize', 18);
ylabel('y', 'FontSize', 18);
zlabel('z', 'FontSize', 18);
    
end
