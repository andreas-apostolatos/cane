function [nodes, coordinates] = getBodyNodes (mesh, BCNodes)
    % [nodes, coordinates] = getBodyNodes (mesh, BCNodes)
    % Returns all non-boundary nodes that have boundary conditions set on
    % them
    % Works on for 2D rectangular meshes

    % Init
    tolerance  = 1e-10;
    
    % Prealloc
    nodes           = BCNodes;
    coordinates     = zeros(length(nodes), 2);
    bodyNodeCounter = 0;
    
    % Get mesh boundaries
    xmin = min(mesh.nodes(:,1));
    xmax = max(mesh.nodes(:,1));
    ymin = min(mesh.nodes(:,2));
    ymax = max(mesh.nodes(:,2));
    
    % Find body nodes
    for k=1:length(BCNodes)
        index = ceil(BCNodes(k)/3);
        if index == int64(index)
            
            index = int64(index);
           
            % Check if point is on the boundary
            bodyNode = true;
            if abs(xmin - mesh.nodes(index, 1)) < tolerance
                bodyNode = false;
            elseif abs(xmax - mesh.nodes(index, 1)) < tolerance
                bodyNode = false;
            elseif abs(ymin - mesh.nodes(index, 2)) < tolerance
                bodyNode = false;
            elseif abs(ymax - mesh.nodes(index, 2)) < tolerance
                bodyNode = false;
            end
            
            % Update outputs
            if bodyNode
                bodyNodeCounter = bodyNodeCounter + 1;
                nodes(bodyNodeCounter) = index;
                coordinates(bodyNodeCounter,:) = mesh.nodes(index, 1:2);
            end
       
        end
    end
    
    % Chop unnecessary/redundant entries
    [nodes,uniqueIndices]   = unique(nodes(1:bodyNodeCounter));
    coordinates             = coordinates(uniqueIndices,:);

end