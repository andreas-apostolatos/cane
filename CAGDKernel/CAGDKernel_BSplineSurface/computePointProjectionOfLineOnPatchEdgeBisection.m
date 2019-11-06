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
function [xi,eta,Q,lambda,minDist,isConvergent] = computePointProjectionOfLineOnPatchEdgeBisection...
    (xiIn,etaIn,PIn,POut,p,Xi,q,Eta,CP,isNURBS,propNewtonRaphson,propBisection)
%% Function documentation
%
% Returns the parametric coordinates and the physical coordinates of the
% projection of a line over the patch boundary, the parameter on the line 
% to be projected, the distance between the line and the patch boundary and
% a flag on whether the iterations have converged. The employed method is
% the bisection method.
%
%             Input :
%        xiIn,etaIn : The parametric coordinates of the edge corner which 
%                     is projected in the current patch
%               PIn : The physical coordinates of the edge corner which is
%                     projected in the current patch
%              POut : The physical coordinates of the edge corner which is
%                     projected outside the patch
%               p,q : The polynomial order of the B-Spline surface patch 
%                     at both parametric directions
%            Xi,Eta : The knot vectors of the B-Spline surface patch in 
%                     both parametric directions
%                CP : The set of Control Point coordinates and weights of 
%                     the B-Spline surface patch
%           isNURBS : Flag on whether the given B-Spline surface patch is a
%                     B-Spline or a NURBS
% propNewtonRaphson : Properties of the Newton-Rapshon for the projection
%                     of a node over the surface B-Spline patch 
%                           .eps : Residual tolerance
%                         .maxIt : Maximum number of iterations
%     propBisection : Properties of the bisection algorithm for the
%                     projection of the line segment over the patch
%                     boundary
%
%            Output :
%            xi,eta : The parametric coordinates on the B-Spline patch 
%                     where the projection of the line segment is found
%                 Q : The physical coordinates on the B-Spline patch 
%                     where the projection of the line segment is found
%            lambda : The line parameter on the line segment where the
%                     closest point onto the B-Spline boundary is found
%           minDist : The distance of the line segment to the patch
%                     boundary
%      isConvergent : Flag on whether the bisection algorithm has converged
%                     or not
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over the bisection iterations
% ->
%    1i. Update the iteration counter
%
%   1ii. Compute the bisection point
%
%  1iii. Computer the orientation vector from the node projected into the patch to the bisection node
%
%   1iv. Compute the projection of the bisection point to the surface B-Spline patch
%
%    1v. Update the points projected in and out the patch
%
%   1vi. Check the convergence criterion
% <-
% 
% 2. Write the output
%
%% Function main body

%% 0. Read input

% Initialize the segment points which are projected inside or outside the
% B-Spline surface patch
PInSaved = PIn;
POutSaved = POut;

% Initialize the iteration counter
counterIter = 0;

% Initialize the initial guess for the Newton-Raphson algorithm for the
% projection of a node on the B-Spline patch
xi0 = xiIn;
eta0 = etaIn;

% Initialize the convergence flag to false
isConvergent = false;

%% 1. Loop over the bisection iterations
while counterIter < propBisection.maxIt && ~isConvergent
    %% 1i. Update the iteration counter
    counterIter = counterIter + 1;
    
    %% 1ii. Compute the bisection point
    P = (PIn + POut)/2;
    
    %% 1iii. Computer the orientation vector from the node projected into the patch to the bisection node
    PInP = P - PIn;
    
    %% 1iv. Compute the projection of the bisection point to the surface B-Spline patch
    [xi,eta,Q,isProjected,noIter] = computeNearestPointProjectionOnBSplineSurface...
        (P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,propNewtonRaphson);
    
    %% 1v. Update the points projected in and out the patch
    if isProjected
        PIn = P;
        QIn = Q;
        xi0 = xi;
        eta0 = eta;
    else
        POut = P;
    end
    
    %% 1vi. Check the convergence criterion
    isConvergent = norm(PInP) < propBisection.eps;
end

%% 2. Write the output
lambda = norm(PIn - PInSaved)/norm(POutSaved - PInSaved);
minDist = norm(PIn - QIn);

end