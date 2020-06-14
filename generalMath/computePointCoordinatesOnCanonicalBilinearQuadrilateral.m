function [uv, isConverged] = ...
    computePointCoordinatesOnCanonicalBilinearQuadrilateral...
    (x, xCorner, newtonRapshon)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the coordinates of a point with coordinates x within the 2D 
% quadrilateral x1-x2-x3-x4 in an arbitrary Cartesian space to the
% canonical bilinear quadrilateral. The applied method for solving the
% non-linear system is Newton-Rapson.
%
%         Input :
%             x : The coordinates of the point into the Cartesian system
%       xCorner : The coordinates of the corner nodes of the quadrilateral
%                 as follows :
%                       xCorner = [x1 x2 x3 x4]
% newtonRapshon : Parameters for the Newton-Raphson iterations :
%                    .eps : Termination criterion for the 2-norm of the
%                           residual
%                  .maxIt : Maximum number of Newton-Rapshon iterations
%
%       Output :
%           uv : The coordinates of the point in the canonical space of the
%                bilinear quadrilateral uv = [u v]
%  isConverged : Flag on the convergence of the Newton iterations
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the Newton iterations
%
%    1i. Compute the bilinear basis functions and its derivatives at u
%
%   1ii. Compute the right-hand side residual
%
%  1iii. Compute the Jacobian
%
%   1iv. Solve the linear system
%
%    1v. Check condition for convergence
%
%   1vi. Update the iteration counter
%
% 2. Check if convergence has been achieved
%
%% Function main body

%% 0. Read input

% Coordinates of point x into the canonical space. Initialize them into the
% center of the canonical bilinear quadrilateral
uv = zeros(2,1);

% Initialize counter
counter = -1;

% Initialize convergence flag to true
isConverged = 1;

%% 1. Loop over all the Newton iterations
while counter <= newtonRapshon.maxIt
    %% 1i. Compute the bilinear basis functions and its derivatives at u
    dN = computeBilinearBasisFunctionsAndFirstDerivatives(uv(1,1),uv(2,1));
    
    %% 1ii. Compute the right-hand side residual
    
    % Compute the physical coordinates of the point with canonical
    % coordinates u
    xPredicted = zeros(2,1);
    DxPredictedDxi = zeros(2,1);
    DxPredictedDeta = zeros(2,1); 
    
    % Loop over all the basis functions and add contributions
    for i = 1:length(dN(:,1))
        % The coordinates of the point
        xPredicted = xPredicted + dN(i,1)*xCorner(:,i);
        
        % The derivatives of the position vector w.r.t. tp xi
        DxPredictedDxi = DxPredictedDxi + dN(i,2)*xCorner(:,i);
        
        % The derivatives of the position vector w.r.t. tp eta
        DxPredictedDeta = DxPredictedDeta + dN(i,3)*xCorner(:,i);
    end
    
    % Compute the residual vector
    residual = x - xPredicted;
    
    %% 1iii. Compute the Jacobian
    J = [DxPredictedDxi(1,1) DxPredictedDeta(1,1)
         DxPredictedDxi(2,1) DxPredictedDeta(2,1)];
    
    %% 1iv. Solve the linear system
    delta = J\residual;
    
    % Update the parametric locations
    uv = uv + delta;
    
    %% 1v. Check condition for convergence
    if norm(residual) <= newtonRapshon.eps
        break;
    end
    
    %% 1vi. Update the iteration counter
    counter = counter + 1;
end

%% 2. Check if convergence has been achieved
if counter == newtonRapshon.maxIt
    isConverged = 0;
end

end
