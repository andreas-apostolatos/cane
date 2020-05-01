function valuesFromInhomDOFs = ...
    computeIGAVctPresrcribedDOFs ...
    (inhomDOFs, valuesFromInhomDOFs, prescribedMotion, CP, parameters, t, ...
    isUniqueOnBoundary)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the vector of the prescribed values on the respective Dirichlet 
% DOFs which are defined by the array prb and the function or value pointer
% prescribedMotion.
%
%              Input :
%                irb : The vector of the DOFs on which inhomogeneous
%                      boundary conditions are applied
%                vrb : The array of the values on the DOFs where
%                      inhomogeneous boundary conditions are applied
%   prescribedMotion : The function or value pointer which defines the 
%                      values of the prescribed DOFs
%                 CP : The set of the Control Point coordinates and weights
%         parameters : Parameters of the underlying field
%                  t : The current time
% isUniqueOnBoundary : Flag on whether the applied prescribed value is
%                      unique, i.e. superior to possible other prescribed
%                      values from intersecting boundaries
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the prescribed DOFs
%
%% Function main body

%% 0. Read input

% number of control points in xi-direction
nxi = length(CP(:,1,1));

%% 1. Loop over all the prescribed DOFs
for i=1:length(inhomDOFs)
    % Get the corresponding Control Point number p and indices CP(i,j)
    p = ceil(inhomDOFs(i)/3);
    jIndex = ceil(p/nxi);
    iIndex = p - (jIndex-1)*nxi;
    
    % Get the Cartesian coordinates of the Control Point
    x = CP(iIndex,jIndex,1);
    y = CP(iIndex,jIndex,2);
    z = CP(iIndex,jIndex,3);
    
    % Assign the prescribed value corresponding to a steady-state regime
    if ~isnumeric(prescribedMotion)
        if ~isUniqueOnBoundary
            valuesFromInhomDOFs(inhomDOFs(i),1) = valuesFromInhomDOFs(inhomDOFs(i),1) + prescribedMotion(parameters,x,y,z,t);
        else
            valuesFromInhomDOFs(inhomDOFs(i),1) = prescribedMotion(parameters,x,y,z,t);
        end
    else
        if ~isUniqueOnBoundary
            valuesFromInhomDOFs(inhomDOFs(i),1) = valuesFromInhomDOFs(inhomDOFs(i),1) + prescribedMotion;
        else
            valuesFromInhomDOFs(inhomDOFs(i),1) = prescribedMotion;
        end
    end
end

end

