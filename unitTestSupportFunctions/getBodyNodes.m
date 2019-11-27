%% Function documentation
%
% Returns all non-boundary nodes that have boundary conditions set on them
% Works only for 2D rectangular meshes 
%
%               Input :
%              fldMsh : Nodes and elements for the fluid mesh
%              homDBC : The global numbering of the DOFs where homogeneous
%                       Dirichlet boundary conditions are applied
%
%
%              Output :
%             homDOFs : Force (drag) on the body in x-direction
%    nodalCoordinates : fill this
%
%% Function main body
function [homDOFs, nodalCoordinates] = getBodyNodes (fldMsh, homDBC)

    % Initialize tolerance
    eps  = 1e-10;
    
    % Preallocate
    homDOFs              = homDBC;
    nodalCoordinates     = zeros(length(homDOFs), 2);
    bodyNodeCounter      = 0;
    
    % Get mesh boundaries (edges or the rectangular domain)
    xmin = min(fldMsh.nodes(:,1));
    xmax = max(fldMsh.nodes(:,1));
    ymin = min(fldMsh.nodes(:,2));
    ymax = max(fldMsh.nodes(:,2));
    
    % Find body nodes
    for k=1:length(homDBC)
        idNode = ceil(homDBC(k)/3);
        if idNode == int64(idNode)
            
            idNode = int64(idNode);
           
            % Check if point is on the boundary
            isBodyNode = true;
            
            if abs(xmin - fldMsh.nodes(idNode, 1)) < eps
                isBodyNode = false;
            elseif abs(xmax - fldMsh.nodes(idNode, 1)) < eps
                isBodyNode = false;
            elseif abs(ymin - fldMsh.nodes(idNode, 2)) < eps
                isBodyNode = false;
            elseif abs(ymax - fldMsh.nodes(idNode, 2)) < eps
                isBodyNode = false;
            end
            
            % Update outputs
            if isBodyNode
                bodyNodeCounter = bodyNodeCounter + 1;
                homDOFs(bodyNodeCounter) = idNode;
                nodalCoordinates(bodyNodeCounter,:) = fldMsh.nodes(idNode, 1:2);
            end
        end
    end
    
    % Chop unnecessary/redundant entries
    [homDOFs,uniqueIndices]  = unique(homDOFs(1:bodyNodeCounter));
    nodalCoordinates         = nodalCoordinates(uniqueIndices,:);

end