function unknownVector = computeNodalVectorIncompressibleFlow2D ... 
    (R, p, q, elementNodalVector)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the element nodal vector vector u = [ux uy p]' at the given 
% parametric location for a 2D incompressible flow problem.
%
%         Input : 
%          RMtx : Matrix containing the basis function values
%           p,q : Polynomial orders of the B-Spline basis along the xi- and
%                 eta- parametric coordinates
%        scalar : element transport vector
%
%        Output :
% unknownVector : The unknown vector u = [ux uy p]' at the given parametric
%                 location
%
% Function layout :
%
% 1. Compute the shape function matrix
%
% 2. Compute the unknown vector as a product u = Rmatrix*d
%
%% Function main body

% Number of basis functions per element
numCPs = (p + 1)*(q + 1);  

% Number of degrees of freedom per element
numDOFs = 3*numCPs;

% Compute the shape function matrix
Rmatrix = zeros(3, numDOFs);

for i=1:numCPs
    Rmatrix(1, 3*i - 2) = R(i);
    Rmatrix(2, 3*i - 1) = R(i);
    Rmatrix(3, 3*i) = R(i);
end

% displacement
unknownVector = Rmatrix*elementNodalVector;

end