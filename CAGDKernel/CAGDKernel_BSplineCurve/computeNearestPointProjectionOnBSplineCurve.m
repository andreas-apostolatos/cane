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
function [xi,Projected,isProjected,noIter] = ...
    computeNearestPointProjectionOnBSplineCurve...
    (P,p,Xi,CP,isNURBS,xi0,propNewtonRaphson)
%% Function documentation
% 
% Returns the nearest point projection of a given point P onto a given
% NURBS-parametrized curve. The applied method for the solution of the
% non-linear system is the Newton-Rapson iterations.
%
%             Input :
%                 P : The point to be projected on the B-Spline surface
%                 q : The polynomial degrees of the NURBS surface
%                Xi : The knot vectors of the NURBS surface
%                CP : Set of Control Point coordinates and weights
%           isNURBS : Flag on whether the patch is a NURBS or a B-Spline
%               xi0 : Initial guess for the curve parameters
% propNewtonRaphson : Properties of the Newton-Rapshon for the projection
%                     of a node over the curve B-Spline patch 
%                           .eps : Residual tolerance
%                         .maxIt : Maximum number of iterations
%
%            Output :
%                xi : The computed curve parameter for the closest point 
%                     projection onto the NURBS curve
%         Projected : The closest point projection onto the NURBS curve
%       isProjected : Flag on wether the algorithm has converged
%            noIter : Number of iterations for convergence
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the Newton-Rapson iterations
% ->
%    1i. Update the iteration counter
%
%   1ii. Find the respective span
%
%  1iii. Compute the NURBS basis functions and their first and second derivatives
%
%   1vi. Compute the position vector and its derivatives for the nurbs curve
%
%  1vii. Check for point coincidence
%
% 1viii. Compute the Jacobian of the nonlinear equation system
%
%   1ix. Compute the residual of the nonlinear equation system
%
%    1x. Compute the solution increment corresponding to the Newton-Raphson method
%
%   1xi. Check if the step size is too small
%
%  1xii. Update the parametric locations
%
% 1xiii. Modify xi if it is out of its interval of definition
%
% <-
%
% 2. Assign the output values
%
%% Function main body

%% 0. Read input

% Initialize counter
counter = 0;

% Initialize the curve parameters
xi = xi0;

% Initialize flag to true
isProjected = true;

% Compute the number of Control Points at each parametric direction
nxi = length(CP(:,1));

%% 1. Loop over all the Newton-Rapson iterations
while counter <= propNewtonRaphson.maxIt
    %% 1i. Update the iteration counter
    counter = counter + 1;
    
    %% 1ii. Find the respective span
    xiSpan = findKnotSpan(xi,Xi,nxi);
    
    %% 1iii. Compute the NURBS basis functions and their first and second derivatives
    noDerivsBasis = 2;
    dR = computeIGABasisFunctionsAndDerivativesForCurve...
        (xiSpan,p,xi,Xi,CP,isNURBS,noDerivsBasis);
    
    %% 1vi. Compute the position vector and its derivatives for the nurbs curve
    noDerivsBaseVct = length(dR(1,:));
    dG = computeBaseVectorsAndDerivativesForBSplineCurve...
        (xiSpan,p,CP,noDerivsBaseVct,dR);
    Projected = dG(:,1);
    
    %% 1vii. Check for point coincidence
    distance = (Projected - P);
    if norm(distance) <= propNewtonRaphson.eps
        break;
    end
    
    %% 1viii. Compute the Jacobian of the nonlinear equation system
    J = dG(:,2)'*dG(:,2) + (Projected - P)'*dG(:,3);
     
    %% 1ix. Compute the residual of the nonlinear equation system
    residual =  (Projected - P)'*dG(:,2);
    if abs(residual) <= propNewtonRaphson.eps
       break;
    end
    
    %% 1x. Compute the solution increment corresponding to the Newton-Raphson method
    deltaXi = J\-residual;
    
    %% 1xi. Check if the step size is too small
    if norm(deltaXi*dG(:,2)) <= propNewtonRaphson.eps
        break;
    end
    
    %% 1xii. Update the parametric locations
    xi = xi + deltaXi;
    
    %% 1xiii. Modify xi if it is out of its interval of definition
    if xi > Xi(length(Xi)) || isnan(xi)
        xi = Xi(length(Xi)); 
    end 
    if xi < Xi(1)
        xi = Xi(1); 
    end
end

%% 2. Assign the output values
noIter = counter;
if noIter > propNewtonRaphson.maxIt
    isProjected = false;
end

end