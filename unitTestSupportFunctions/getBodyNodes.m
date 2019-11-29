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
%              Output :
%    nodalCoordinates : Nodes and elements from the fluid mesh that lie on
%                       a body 
%             homDOFs : The global numbering of the DOFs that lie on a body
%                       where homogeneous Dirichlet boundary conditions are
%                       applied
%
%% Function main body
function homDOFs = getBodyNodes (fldMsh, homDBC)

    % Initialize tolerance for node finding
    eps  = 1e-10;
    
    % Preallocate variables
    homDOFs              = homDBC;
    bodyNodeCounter      = 0;
    
    % Get mesh boundaries (edges or the rectangular domain)
    xmin = min(fldMsh.nodes(:,1));
    xmax = max(fldMsh.nodes(:,1));
    ymin = min(fldMsh.nodes(:,2));
    ymax = max(fldMsh.nodes(:,2));
    
    % Find body nodes
    for k=1:length(homDBC)
        
        idNode = ceil(homDBC(k)/3);   

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
        end
        
    end
    
    % Chop unnecessary/redundant entries
    homDOFs  = unique(homDOFs(1:bodyNodeCounter));

end