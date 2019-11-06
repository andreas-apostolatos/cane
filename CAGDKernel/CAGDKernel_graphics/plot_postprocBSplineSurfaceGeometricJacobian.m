%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_postprocBSplineSurfaceGeometricJacobian(p,q,Xi,Eta,CP,isNURBS)
%% Function documentation
%
% Plots the reference configuration of an isogeometric Kirchhoff-Love shell
% together with the selected postprocessing resultant.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
%
%         Output :
%                  graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the resultant array to be used for the visualization
%
% 2. Visualize the Jacobian determinant and the knots over the domain
%
%% Function main body

%% 0. Read input

% Define the grid points over which to compute the determinant of the
% Jacobian
xiGrid = 50;
etaGrid = 50;

% Number of knots in xi,eta-direction
mxi = length(Xi);
meta = length(Eta);

% Assign a tolerance value
tol = 10e-10;

% Compute incremental steps for the resultant computation over the domain
% incremental step for eta:
deta = (Eta(meta)-Eta(1))/(xiGrid-1);

% incremental step for xi:
dxi = (Xi(mxi)-Xi(1))/(etaGrid-1);

% Initialize array of the Jacobian determinant to be visualized
detJ = zeros(xiGrid,etaGrid);

% Cartesian image of the parameter space
P = zeros(xiGrid,etaGrid,3);

% Number of Control Points in xi-,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

%% 1. Compute the resultant array to be used for the visualization

% counting index in eta-direction
etaCounter = 1;  

% Initialize coordinate in eta-direction
eta = Eta(1);

% Loop over all the parametric coordinates in eta-direction
while eta <= Eta(meta)+tol
    % Find the span in the eta-direction
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    % Initialize coordinate in xi-direction
    xi = Xi(1);
    
    % Initialize counter in xi-direction
    xiCounter = 1;
    
    % Loop over all the parametric coordinates in xi-direction
    while xi <= Xi(mxi)+tol
        % Find the span in xi-direction
        xiSpan = findKnotSpan(xi,Xi,nxi);
        
        % Compute the basis functions and their first derivatives
        nDrv = 1;
        dR = computeIGABasisFunctionsAndDerivativesForSurface...
            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
        
        % Compute the base vectors
        nDrvBaseVct = 0;
        [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
            (xiSpan,p,etaSpan,q,CP,nDrvBaseVct,dR);

        % Compute the surface normal (third covariant base vector not 
        % normalized) for the reference configuration
        A3Tilde = cross(A1(:,1),A2(:,1));

        % Compute the surface normal vector
        A3 = A3Tilde/norm(A3Tilde);
        
        % Compute the Cartesian image of the paratric point
        P(xiCounter,etaCounter,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
        
        % Compute the determinant of the Jacobian determinant on the Gauss 
        % point
        detJ(xiCounter,etaCounter) = norm(A3Tilde);
        
        % Update the counter in xi-direction
        xiCounter = xiCounter + 1;
        
        % Update the parametric coordinate in xi-direction
        xi = xi + dxi;
    end
    % Update the counter in eta-direction
    etaCounter = etaCounter + 1;
    
    % Update the parametric coordinate in eta-direction
    eta = eta + deta;
end

%% 2. Visualize the Jacobian determinant and the knots over the domain

% Plot the resultant over the reference configuration
surf(P(:,:,1),P(:,:,2),P(:,:,3),detJ(:,:),'EdgeColor','none');
hold on;

% Plot the element boundaries on the undeformed geometry
isUndeformed = 0;
plot_knotsForBSplineSurfaceOnCartesianSpace...
    (p,q,Xi,Eta,CP,isNURBS,isUndeformed,xiGrid,etaGrid);
hold off;

end