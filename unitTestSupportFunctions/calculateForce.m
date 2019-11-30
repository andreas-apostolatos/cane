function [Fx, Fy, Fz] = calculateForce (FComplete, fldMsh, homDBC, parameters, noDimensions)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
% Andreas Apostolatos
%
%% Function documentation
%
% Calculate force on a single body with no-slip boundary conditions 
% inside a rectangular domain
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
%                  Fx : Force on the body in x-direction
%                  Fy : Force on the body in y-direction
%                  Fz : Force on the body in z-direction
%
%% Function main body
    
    % Get node indices on the body
    idNodesOnBody = getBodyNodesOnRectangularDomain(fldMsh, homDBC);
    
    % Reshape the input forces to [Fx, Fy, Fz]
    FComplete = -reshape(FComplete, [(noDimensions +1),                 ...
                   length(FComplete)/(noDimensions +1)])';
    
    % Select forces on body nodes
    F = FComplete(idNodesOnBody, 1:noDimensions) * parameters.rho;
    
    % Sum the individual components
    Fx = sum(F(:,1));
    Fy = sum(F(:,2));
    Fz = 0;
    
    % Sum Fz if we have a 3D problem
    if noDimensions == 3
        Fz = sum(F(:,3));
    end
    
end