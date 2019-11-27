%% Function documentation
%
% Calculate drag force on a single body with no-slip boundary conditions in a rectangular domain
%
%               Input :
%           FComplete : The complete force vector
%              fldMsh : Nodes and elements for the fluid mesh
%              homDBC : The global numbering of the DOFs where homogeneous
%                       Dirichlet boundary conditions are applied
%          parameters : Flow parameters
%
%
%              Output :
%                  Fx : Force (drag) on the body in x-direction
%
%% Function main body
function Fx = calculateDrag (FComplete, fldMsh, homDBC, parameters)

    % Get node indices on the body
    indices = getBodyNodes(fldMsh, homDBC);
    
    % Reshape the input forces to [Fx, Fy, Fz]
    FComplete = -reshape(FComplete, [3, length(FComplete)/3])';
    
    % Select forces on body nodes
    F = FComplete(indices, 1:2) * parameters.rho;

    Fx = sum(F(:,1));
    
end