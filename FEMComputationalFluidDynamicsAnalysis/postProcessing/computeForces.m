function Force = computeForces (FComplete, parameters, nodesDomain, noDimensions)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% Calculate force on a single body
%
%               Input :
%           FComplete : The complete force vector
%              fldMsh : Nodes and elements for the fluid mesh
%              homDBC : The global numbering of the DOFs where homogeneous
%                       Dirichlet boundary conditions are applied
%          parameters : Flow parameters
%        noDimensions : dimensionality of the problem 
%                       2 for 2D problems, 3 for 3D problems
%   
%              Output :
%            Force(1) : Force on the body in x-direction
%            Force(2) : Force on the body in y-direction
%            Force(3) : Force on the body in z-direction
%
%% Function main body
    
    % Reshape the input forces to [Fx, Fy, Fz]
    FComplete = -reshape(FComplete, [(noDimensions +1),                 ...
                   length(FComplete)/(noDimensions +1)])';
    
    % Select forces on body nodes
    F = FComplete(nodesDomain, 1:noDimensions) * parameters.rho;
    
    % Sum the individual components
    Force(1) = sum(F(:,1));
    Force(2) = sum(F(:,2));
    Force(3) = 0;
    
    % Sum Fz if we have a 3D problem
    if noDimensions == 3
        Force(3) = sum(F(:,3));
    end
    
end