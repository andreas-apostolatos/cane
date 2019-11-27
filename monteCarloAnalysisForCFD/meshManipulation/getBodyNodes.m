function [homDOFs, nodalCoordinates] = getBodyNodes (mesh, homDBC)
    % [nodes, coordinates] = getBodyNodes (mesh, BCNodes)
    % Returns all non-boundary nodes that have boundary conditions set on
    % them
    % Works on for 2D rectangular meshes

    % Init
    eps  = 1e-10;
    
    % Prealloc
    homDOFs           = homDBC;
    nodalCoordinates     = zeros(length(homDOFs), 2);
    bodyNodeCounter = 0;
    
    % Get mesh boundaries
    xmin = min(mesh.nodes(:,1));
    xmax = max(mesh.nodes(:,1));
    ymin = min(mesh.nodes(:,2));
    ymax = max(mesh.nodes(:,2));
    
    % Find body nodes
    for k=1:length(homDBC)
        idNode = ceil(homDBC(k)/3);
        if idNode == int64(idNode)
            
            idNode = int64(idNode);
           
            % Check if point is on the boundary
            isBodyNode = true;
            if abs(xmin - mesh.nodes(idNode, 1)) < eps
                isBodyNode = false;
            elseif abs(xmax - mesh.nodes(idNode, 1)) < eps
                isBodyNode = false;
            elseif abs(ymin - mesh.nodes(idNode, 2)) < eps
                isBodyNode = false;
            elseif abs(ymax - mesh.nodes(idNode, 2)) < eps
                isBodyNode = false;
            end
            
            % Update outputs
            if isBodyNode
                bodyNodeCounter = bodyNodeCounter + 1;
                homDOFs(bodyNodeCounter) = idNode;
                nodalCoordinates(bodyNodeCounter,:) = mesh.nodes(idNode, 1:2);
            end
       
        end
    end
    
    % Chop unnecessary/redundant entries
    [homDOFs,uniqueIndices]   = unique(homDOFs(1:bodyNodeCounter));
    nodalCoordinates             = nodalCoordinates(uniqueIndices,:);

end