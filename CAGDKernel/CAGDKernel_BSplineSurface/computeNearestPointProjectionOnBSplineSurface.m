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
function [xi,eta,Projected,isProjected,noIter] = ...
    computeNearestPointProjectionOnBSplineSurface...
    (P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,propNewtonRaphson)
%% Function documentation
% 
% Returns the nearest point projection of a given point P onto a given
% NURBS-parametrized surface. The applied method for the solution of the
% non-linear system is the Newton-Rapson iterations.
%
% Source : Tianyang Wang (2010),'Investigation and Implementation of
%          Non-Matching Grids Data Transfer', Master Thesis, Technische
%          Universität München
%
%             Input :
%                 P : The point to be projected on the B-Spline surface
%               p,q : The polynomial degrees of the NURBS surface
%            Xi,Eta : The knot vectors of the NURBS surface
%                CP : Set of Control Point coordinates and weights
%           isNURBS : Flag on whether the patch is a NURBS or a B-Spline
%          xi0,eta0 : Initial guess for the surface parameters
% propNewtonRaphson : Properties of the Newton-Rapshon for the projection
%                     of a node over the surface B-Spline patch 
%                           .eps : Residual tolerance
%                         .maxIt : Maximum number of iterations
%
%            Output :
%            xi,eta : The computed surface parameters for the closest point 
%                     projection onto the NURBS surface
%         Projected : The closest point projection onto the NURBS surface
%       isProjected : Flag on wether the algorithm has converged
%             noIter : Number of iterations for convergence
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the Newton-Rapson iterations
%
%    1i. Update the iteration counter
%
%   1ii. Find the respective span
%
%  1iii. Compute the NURBS basis functions and their first and second derivatives
%
%   1iv. Compute the point on the NURBS surface
%
%    1v. Check for point coincidence
%
%   1vi. Compute the base vectors and their derivatives for the NURBS surface at (xi,eta)
%
%  1vii. Check the condition of zero cosine
%
% 1viii. Compute the Jacobian of the nonlinear equation system
%
%   1ix. Compute the residual of the nonlinear equation system
%
%    1x. Compute the solution increment corresponding to the Newton-Raphson method
%
%   1xi. Modify (xi,eta) if they are out of their intervals of restriction
%
%% Function main body

%% 0. Read input

% Initialize counter
counter = 0;

% Initialize the surface parameters
xi = xi0;
eta = eta0;

% Initialize flag to true
isProjected = 1;

% Tolerance up to which the Jacobian is assumed to be singular
tolJacobian = 1e-6;

% Clamp the parametric coordinates
clampXi = false;
clampEta = false;

% Compute the number of Control Points at each parametric direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

%% 1. Loop over all the Newton-Rapson iterations
while counter <= propNewtonRaphson.maxIt
    %% 1i. Update the iteration counter
    counter = counter + 1;
    
    %% 1ii. Find the respective span
    xiSpan = findKnotSpan(xi,Xi,nxi);
    etaSpan = findKnotSpan(eta,Eta,neta); 
    
    %% 1iii. Compute the NURBS basis functions and their first and second derivatives
    dR = computeIGABasisFunctionsAndDerivativesForSurface...
        (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,2);
    
    %% 1iv. Compute the point on the NURBS surface
    Projected = computeCartesianCoordinatesOfAPointOnBSplineSurface...
        (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
    
    %% 1v. Check for point coincidence
    distance = (Projected - P);
    if norm(distance) <= propNewtonRaphson.eps
        break;
    end
    
    %% 1vi. Compute the base vectors and their derivatives for the NURBS surface at (xi,eta)
    [dG1,dG2] = computeBaseVectorsAndDerivativesForBSplineSurface(xiSpan,p,etaSpan,q,CP,1,dR);
    
    %% 1vii. Check the condition of zero cosine
    
    % Compute the cosine with respect to the u-parametric coordinate
    xiCos = abs(dG1(:,1)'*distance)/norm(dG1(:,1))/norm(Projected-P);
    
    % Compute the cosine with respect to the u-parametric coordinate
    etaCos = abs(dG2(:,1)'*distance)/norm(dG2(:,1))/norm(Projected-P);
    
    % Check the orthogonality condition
    if xiCos <= propNewtonRaphson.eps && etaCos <= propNewtonRaphson.eps
       break;
    end    
    
    %% 1viii. Compute the Jacobian of the nonlinear equation system
    J = [norm(dG1(:,1))^2 + distance'*dG1(:,2)    dG1(:,1)'*dG2(:,1) + distance'*dG1(:,3)
         dG2(:,1)'*dG1(:,1) + distance'*dG1(:,3)  norm(dG2(:,1))^2 + distance'*dG2(:,2)];
     
    %% 1ix. Compute the residual of the nonlinear equation system
    residual =  [distance'*dG1(:,1) 
                 distance'*dG2(:,1)];
    
    %% 1x. Compute the solution increment corresponding to the Newton-Raphson method
    if abs(J(1,1)) < tolJacobian || clampXi
        delta = zeros(2,1);
        delta(1,1) = 0;
        delta(2,1) = residual(2,1)/J(2,2);
        clampXi = false;
        clampEta = true;
    elseif abs(J(2,2)) < tolJacobian || clampEta
        delta = zeros(2,1);
        delta(1,1) = residual(1,1)/J(1,1);
        delta(2,1) = 0;
        clampXi = true;
        clampEta = false;
    else
        delta = J\-residual;
    end
    
    % Check if the step size is too small
    condition = norm(delta(1,1)*dG1(:,1) + delta(2,1)*dG2(:,1));
    if condition <= propNewtonRaphson.eps
        break;
    end
    
    % Update the parametric locations
    xi = xi + delta(1,1);
    eta = eta + delta(2,1);
    
    %% 1xi. Modify (xi,eta) if they are out of their intervals of restriction
    if xi > Xi(length(Xi)) || isnan(xi)
        xi = Xi(length(Xi)); 
    end 
    if xi < Xi(1)
        xi = Xi(1); 
    end
    if eta > Eta(length(Eta)) || isnan(eta)
        eta = Eta(length(Eta)); 
    end
    if eta < Eta(1)
        eta = Eta(1); 
    end
end

%% 2. Assign the output values
noIter = counter;
if noIter > propNewtonRaphson.maxIt
    isProjected = 0;
end

end