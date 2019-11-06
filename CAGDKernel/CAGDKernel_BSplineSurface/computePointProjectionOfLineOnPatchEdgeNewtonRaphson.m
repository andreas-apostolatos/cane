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
function [t,lambda,distance,noNRIter,isConvergent] = ...
    computePointProjectionOfLineOnPatchEdge...
    (t0,P1,P2,p,Xi,q,Eta,CP,isNURBS,xiB,etaB,newtonRaphson)
%% Function documentation
%
% Solves the nonlinear problem of the projection of a straight line in 3D
% over the boundary of a B-Spline patch. The projection is based on the
% orthogonality condition of the connecting line from the patch boundary
% onto the straight line over the B-Spline surface.
%
%         Input :
%            t0 : The initial guess for the curve parameter over the 
%                 B-Spline boundary
%         P1,P2 : The corners of the straight line to be projected over the
%                 B-Spline boundary
%           p,q : The polynomial orders of the B-Spline patch in xi-,eta-
%                 parametric directions
%        Xi,Eta : The knot vectors of the B-Spline patch in xi-,eta-
%                 parametric directions
%            CP : The set of the Control Point coordinates and weights for
%                 the B-Spline patch
%       isNURBS : Flag on whether the basis of the patch is a NURBS or not
%      xiB,etaB : The extension of the B-Spline boundary
% newtonRaphson : Properties of the Newton-Rapshon iterations for the
%                 projection of the line onto the B-Spline boundary
%                       .eps : The tolerance for convergence to the
%                              Newton-Rapshon iterations
%                     .maxIt : Maximum number of iterations
%   
%        Output :
%             t : The computed parameter onto the B-Spline boundary
%        lambda : The computed parameter onto the straight line between P1
%                 and P2 measured from point P1
%      distance : The computed orthogonal distance between the straight
%                 line and the patch boundary
%      noNRIter : Number of Newton-Raphson iterations
%  isConvergent : Flag on whether the Newton-Raphson iterations converged
%                 or not
%
% Function layout :
%
% 0. Read input
%
% 1. Initialize the parameters on the B-Spline boundary
%
% 2. Loop over all the Newton-Rapshon iterations
% ->
%    2i. Update the counter
%
%   2ii. Find the knot spans
%
%  2iii. Compute the basis funtions and the base vectors as well as their derivatives
%
%   2iv. Compute the Cartesian coordinates of the currently found point
%
%    2v. Compute the base vectors and their derivatives
%
%   2vi. Compute the normal of the B-Spline surface
%
%  2vii. Compute the normal between the line to be cut by the edge and the normal to the B-Spline surface
%
% 2viii. Compute the distance vector
%
%   2ix. Compute the residual
%
%    2x. Compute the derivative of the residual in the direction of the patch boundary
%
%   2xi. Check the termination condition on the residual
%
%  2xii. Update the parametric coordinates
%
% 2xiii. Clamp the knots if they jump out of the knot vector interval
% <-
% 
% 3. Compute the parameter on the straight segment connecting the vertivces P1 and P2
%
% 4. Write out the output
%
%% Function main body

%% 0. Read input

% Initialize counter for the Newton-Raphson iterations
noNRIter = 0;

% Initialize convergence flag
isConvergent = false;

% Define an interval tolerance for the lambda parameters
tol = 1e-3;

%% 1. Initialize the parameters on the B-Spline boundary
if etaB(1) == etaB(2)
    xi = t0;
    eta = etaB(1);
else
    xi = xiB(1);
    eta = t0;
end

%% 2. Loop over all the Newton-Rapshon iterations
while ~isConvergent && noNRIter < newtonRaphson.maxIt
    %% 2i. Update the counter
    noNRIter = noNRIter + 1;
    
    %% 2ii. Find the knot spans
    xiSpan = findKnotSpan(xi,Xi,length(CP(:,1,1)));
    etaSpan = findKnotSpan(eta,Eta,length(CP(1,:,1)));
    
    %% 2iii. Compute the basis funtions and the base vectors as well as their derivatives
    dR = computeIGABasisFunctionsAndDerivativesForSurface...
        (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,2);
    
    %% 2iv. Compute the Cartesian coordinates of the currently found point
    Q = computeCartesianCoordinatesOfAPointOnBSplineSurface...
        (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
    
    %% 2v. Compute the base vectors and their derivatives
    [dGXi,dGEta] = computeBaseVectorsAndDerivativesForBSplineSurface...
        (xiSpan,p,etaSpan,q,CP,1,dR);
    
    %% 2vi. Compute the normal of the B-Spline surface
    G3 = cross(dGXi(:,1),dGEta(:,1));
    
    %% 2vii. Compute the normal between the line to be cut by the edge and the normal to the B-Spline surface
    P1P2 = P2 - P1;
    n = cross(P1P2,G3);
    
    %% 2viii. Compute the distance vector
    P1Q = Q - P1;
    
    %% 2ix. Compute the residual
    Residual = n'*P1Q;
    
    %% 2x. Compute the derivative of the residual in the direction of the patch boundary
    if etaB(1) == etaB(2)
        product1 = cross(dGXi(:,3),dGEta(:,1));
        product2 = cross(dGXi(:,1),dGXi(:,2));
        dResidual = dGXi(:,1)'*n;
    else
        product1 = cross(dGXi(:,2),dGEta(:,1));
        product2 = cross(dGXi(:,1),dGEta(:,2));
        dResidual = dGEta(:,1)'*n;
    end
    product1 = product1 + product2;
    product3 = cross(P1P2,product1);
    dResidual = dResidual + P1Q'*product3;
    
    %% 2xi. Check the termination condition on the residual
    if norm(Residual) < newtonRaphson.eps
        isConvergent = true;
        break;
    end
    
    %% 2xii. Update the parametric coordinates
    if etaB(1) == etaB(2)
        xi = xi - Residual/dResidual;
    else
        eta = eta - Residual/dResidual;
    end
    
    %% 2xiii. Clamp the knots if they jump out of the knot vector interval
    if xi < Xi(1);
        xi = Xi(1);
    end
    if xi > Xi(end);
        xi = Xi(end);
    end
    if eta < Eta(1);
        eta = Eta(1);
    end
    if eta > Eta(end);
        eta = Eta(end);
    end
end

%% 3. Compute the parameter on the straight segment connecting the vertivces P1 and P2
N = G3/norm(G3);
h10 = P1Q'*N;
P1P = P1Q - h10*N;
lambda = norm(P1P)/norm(P1P2);
if lambda < 0 && lambda > 0 - tol
    lambda = 0;
end
if lambda > 1 && lambda < 1 + tol
    lambda = 1;
end

%% 4. Write out the output

% Compute the distance vector between P and Q
QP = P1 + lambda*(P2 - P1) - Q;

% The orthogonal distance between the straight line and the patch boundary
distance = norm(QP);

% The parameter on the patch boundary
if etaB(1) == etaB(2)
    t = xi;
else
    t = eta;
end

end