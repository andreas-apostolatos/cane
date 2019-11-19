function [xi,eta,Projected,isProjected,noIter] = ...
    computeNearestPointProjectionOnBSplineSurface...
    (P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,propNewtonRaphson)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
%% Function documentation
% 
% Returns the nearest point projection of a given point P onto a given
% NURBS-parametrized surface. The applied method for the solution of the
% non-linear system is the Newton-Rapson iterations.
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
%    1x. Compute the solution increment corresponding to the Newton-Raphson method and take into consideration special cases around singularities
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
isProjected = true;

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
    
    %% 1vii. Compute the residual and check the condition for convergence
    
    % Compute the residual based on the orthogonality condition
    residualXi = - distance'*dG1(:,1);
    residualEta = - distance'*dG2(:,1);
    fprintf("\nresidualXi: %.15d\n", residualXi);
    fprintf("residualEta: %.15d\n\n", residualEta);
    
    % Compute the cosine with respect to the xi-parametric coordinate
    xiCos = abs(residualXi)/norm(dG1(:,1))/norm(distance);
    
    % Compute the cosine with respect to the eta-parametric coordinate
    etaCos = abs(residualEta)/norm(dG2(:,1))/norm(distance);
    
    % Check the orthogonality condition
    if xiCos <= propNewtonRaphson.eps && etaCos <= propNewtonRaphson.eps
       break;
    end    
    
    %% 1viii. Compute the Jacobian of the nonlinear equation system
    J_00 = norm(dG1(:,1))^2 + distance'*dG1(:,2);
    J_01 = dG1(:,1)'*dG2(:,1) + distance'*dG1(:,3);
    J_11 = norm(dG2(:,1))^2 + distance'*dG2(:,2);   
    
    %% 1x. Compute the solution increment corresponding to the Newton-Raphson method and take into consideration special cases around singularities
    
    % Check the conditioning of the Jacobian matrix
    conditionFirstRowZero = false;
    if abs(J_00) < propNewtonRaphson.eps && J_01 < propNewtonRaphson.eps
        conditionFirstRowZero = true;
    end
    conditionSecondRowZero = false;
    if abs(J_01) < propNewtonRaphson.eps && abs(J_11) < propNewtonRaphson.eps
        conditionSecondRowZero = true;
    end
    conditionFirstColumnZero = false;
    if abs(J_00) < propNewtonRaphson.eps && abs(J_01) < propNewtonRaphson.eps
        conditionFirstColumnZero = true;
    end
    conditionSecondColumnZero = false;
    if abs(J_01) < propNewtonRaphson.eps && abs(J_11) < propNewtonRaphson.eps
        conditionSecondColumnZero = true;
    end
    
    % Check if the system is solvable by checking the condition of the diagonal entries
    isSystemInvertible = true;
    if conditionFirstRowZero || conditionSecondRowZero || conditionFirstColumnZero || conditionSecondColumnZero
        isSystemInvertible = false;
    end
    
    % Solve the 2x2 linear equation system and take into account special cases where singularities occur
    if isSystemInvertible
        det_J = J_00 * J_11 - J_01 * J_01;
        delta_xi = - (residualEta * J_01 - residualXi * J_11)/det_J;
        delta_eta = - (residualXi * J_01 - residualEta * J_00)/det_J;
    else
        if conditionFirstRowZero
            delta_xi = residualEta/J_11;
            delta_eta = 0.0;
        elseif conditionSecondRowZero
            delta_xi = residualXi/J_00;
            delta_eta = 0.0;
        elseif conditionFirstColumnZero
            delta_xi = 0.0;
            delta_eta = (residualXi + residualEta)/(J_01 + J_11);
        elseif conditionSecondColumnZero
            delta_xi = (residualXi + residualEta)/(J_00 + J_01);
            delta_eta = 0.0;
        end
    end
    
    % Check if the step size is too small
    condition = norm(delta_xi*dG1(:,1) + delta_eta*dG2(:,1));
    if condition <= propNewtonRaphson.eps
        break;
    end
    
    % Update the parametric locations
    xi = xi + delta_xi;
    eta = eta + delta_eta;
    
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
    isProjected = false;
end

end
