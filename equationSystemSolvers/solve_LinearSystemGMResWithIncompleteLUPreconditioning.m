function [x, isConverged] = solve_LinearSystemGMResWithIncompleteLUPreconditioning...
    (A, b, x)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the solution to a system of linear equations Ax=b using the GMRES 
% (Generalized Minimal Residual) solver with an incomplete LU factorization
%
%       Input :
%           A : The left hand side matrix
%           b : The right hand side vector
%           x : The initial guess of the GMRes iterations
%
%  Output :
%           x : The solution vector
% isConverged : Flag on whether the linear solver has converged
%
%   Function layout :
%
% 0. Read input
%
% 1. Perform the LU factorization
%
% 2. Loop over the GMRes iterations
% ->
%    2i. Solve the GMRes step
%
%   2ii. Check convergence criterion
% <-
% 
% 3. Check if the GMRes solver conveged
%
%% Function main body

%% 0. Read input

% set the parameters of the iterations
tolerance = 1e-12;
maxIt = 1e2;

% Initialize the iteration counter
counterIter = 1;

% Initialize convergence flag
isConverged = true;

%% 1. Perform the LU factorization
[L, U] = ilu(sparse(A));

%% 2. Loop over the GMRes iterations
while counterIter < maxIt
    %% 2i. Solve the GMRes step
    [x, ~, relResidual] = gmres(A, b, [], tolerance, [], L, U, x);
    counterIter = counterIter + 1;
    
    %% 2ii. Check convergence criterion
    if relResidual < tolerance
        break;
    end
end

%% 3. Check if the GMRes solver conveged
if counterIter == maxIt
    warning( 'GMRes solver did not converge in the given number of iterations!' );
    isConverged = false;
end

end
