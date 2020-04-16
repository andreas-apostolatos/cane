function [inhomDOFs, valuesInhomDOFs] = ...
    computeVectorOfPreMotionOnInhomDoFsForIncompressibleFlow2D ... 
    (inhomDOFs, valuesInhomDOFs, xi, eta, direction, prescribedMotion, ...
    isUniqueOnBoundary, CP, parameters, t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the prescribed velocity field on the inhomogeneous Dirichlet
% boundary for the case of a 2D incompressible flow problem.
%
%               Input :
%           inhomDOFs : The outdated global numbering of the DoFs with 
%                       inhomogeneous Dirichlet boundary conditions
%     valuesInhomDOFs : The outdated vector of the prescribed motion on the
%                       DoFs with inhomogeneous Dirichlet boundary
%                       conditions
%              xi,eta : region to be fixed (e.g. xi = [0 1], eta = [1 1])
%           direction : Direction of the prescribed motion
%    prescribedMotion : Prescribed motion in the chosen direction
%  isUniqueOnBoundary : Flag determining whether this inhomogeneous
%                       Dirichlet boundary condition is unique on the given
%                       boundary or it exists in combination with the
%                       previous ones
%                  CP : Set of control point coordinates and weigths
%          parameters : Flow parameters
%                   t : The current simulation time
%
%              Output : 
%                 vrb : The updated vector of the prescribed motion on the
%                       DoFs with inhomogeneous Dirichlet boundary
%                       conditions
%                 prb : The updated global numbering of the DoFs with
%                       inhomogeneous Dirichlet boundary conditions
%
% Function layout :
%
% 0. Read input
%
% 1. Iterate and add new supports preserving the old ones
%
% 2. Loop over the supports and delete double entries if any
%
%% Function main body

%% 0. Read input

% counters
prbCounter = length(inhomDOFs) + 1;
vrbCounter = length(valuesInhomDOFs) + 1;

% Check if the given input is consistent
if (prbCounter ~= vrbCounter); error('The given vectors do not match'); end
    
% number of control points in u,v-direction
numCPs_xi = length(CP(:,1,1));
numCPs_eta = length(CP(1,:,1));

%% 1. Iterate and add new supports preserving the old ones
for j = eta(1)*(numCPs_eta-1)+1:eta(2)*(numCPs_eta-1)+1
    for i = xi(1)*(numCPs_xi-1)+1:xi(2)*(numCPs_xi-1)+1
        % Get the DoF to be prescribed
        inhomDOFs(prbCounter) = 3*((j-1)*numCPs_xi + i-1) + direction;
        
        % Round to nearest integer
        inhomDOFs(prbCounter) = round(inhomDOFs(prbCounter));
        
        % Add the prescribed to the DoF value
        if isnumeric(prescribedMotion)
            valuesInhomDOFs(prbCounter) = prescribedMotion;
        else
            % get the corresponding Control Point number p and indices CP(i,j)
            p = ceil(inhomDOFs(prbCounter)/3);
            jIndex = ceil(p/numCPs_xi);
            iIndex = p - (jIndex-1)*numCPs_xi;
            
            % Get the Cartesian coordinates of the Control Point
            x = CP(iIndex,jIndex,1);
            y = CP(iIndex,jIndex,2);
            z = CP(iIndex,jIndex,3);
            
            % Assign the prescribed value corresponding to a steady-state
            % regime
            valuesInhomDOFs(prbCounter) = prescribedMotion ...
                (parameters, x, y, z, t);
        end
        
        % Update counter
        prbCounter = prbCounter + 1;
    end
end

%% sort rb and delete double entries 

% Sort the values
[inhomDOFs,indexAray] = sort(inhomDOFs);
valuesInhomDOFs = valuesInhomDOFs(indexAray);

% Initialize counter
i=1;
% Loop over the supports and delete double entries if any
while i < length(inhomDOFs)
    if inhomDOFs(i)==inhomDOFs(i+1)
        % Delete the double entry in the prb array
        inhomDOFs(i+1) = [];  
        
        % Delete the double entry in the vrb array but add its value to the
        % prescribed motion
        if ~isUniqueOnBoundary
            valuesInhomDOFs(i) = valuesInhomDOFs(i) + valuesInhomDOFs(i+1);
        end
        valuesInhomDOFs(i+1) = [];
        
        % Decrease counter
        i = i - 1;  
    end
    % Update counter
    i = i + 1;
end
    
end