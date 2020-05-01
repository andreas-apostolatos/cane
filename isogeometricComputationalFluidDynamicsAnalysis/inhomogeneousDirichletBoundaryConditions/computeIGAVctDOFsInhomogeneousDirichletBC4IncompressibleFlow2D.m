function inhomDOFs = ...
    computeIGAVctDOFsInhomogeneousDirichletBC4IncompressibleFlow2D ...
    (inhomDOFs, xi, eta, direction, CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function main body
%
% Returns the global numbering of the DOFs on which inhomogeneous Dirichlet
% boundary conditions are prescribed.
%
%     Input :
% inhomDOFs : The existing vector with the global numbering of the DOFs where
%             inhomogeneous Dirichlet boundary conditions are applied
%    xi,eta : The surface parameters on the NURBS patch
% direction : The direction of the DOF in the global coordinate system
%        CP : The set of Control Point coordinates and weights
%
%    Output :
% inhomDOFs : The updated vector with the global numbering of the DOFs 
%             where inhomogeneous Dirichlet boundary conditions are applied
%
% Function layout :
%
% 0. Read input
%
% 1. Iterate and add new supports preserving the old ones
%
% 2. Sort the output array and delete double entries
%
%% Function main body

%% 0. Read input

% counters
prbCounter = length(inhomDOFs) + 1;
    
% number of control points in u,v-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

%% 1. Iterate and add new supports preserving the old ones
for j = eta(1)*(neta-1)+1:eta(2)*(neta-1)+1
    for i = xi(1)*(nxi-1)+1:xi(2)*(nxi-1)+1
        % Get the DoF to be prescribed
        inhomDOFs(prbCounter) = 3*((j-1)*nxi + i-1) + direction;
        
        % Round to nearest integer
        inhomDOFs(prbCounter) = round(inhomDOFs(prbCounter));
        
        % Update counter
        prbCounter = prbCounter + 1;
    end
end

%% 2. Sort the output array and delete double entries

% Sort the values
inhomDOFs = sort(inhomDOFs);

% Delete double entries
inhomDOFs = unique(inhomDOFs);

end

