function homDOFs = getBodyNodesOnRectangularDomain (fldMsh, homDBC)
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
% Returns all nodes that do not lie on the rectangular boundary and have 
% no-slip boundary conditions set on them, works only for 2D od 3D  
% rectangular meshes where body of interest does not touch the boundary
% Example:  
%
%               Input :
%              fldMsh : Nodes and elements for the fluid mesh
%              homDBC : The global numbering of the DOFs where homogeneous
%                       Dirichlet boundary conditions are applied
%
%              Output :
%             homDOFs : The global numbering of the DOFs that lie on a body
%                       where homogeneous Dirichlet boundary conditions are
%                       applied
%
%% Function main body

    % Initialize tolerance for node finding
    eps  = 1e-10;
    
    % Preallocate variables
    homDOFs = homDBC;
    bodyNodeCounter = 0;
    
    % Get mesh boundaries 2D (edges or the rectangular domain)
    xmin = min(fldMsh.nodes(:,1));
    xmax = max(fldMsh.nodes(:,1));
    ymin = min(fldMsh.nodes(:,2));
    ymax = max(fldMsh.nodes(:,2));
    zmin = min(fldMsh.nodes(:,3));
    zmax = max(fldMsh.nodes(:,3));
    
    % check if it is a 2D or 3D problem
    is3D = true;
    if abs(zmin - zmax) < eps
        is3D = false;
    end
    
    % Find body nodes
    for k=1:length(homDBC)
        
        % find the id of a node
        idNode = ceil(homDBC(k)/3);   
        % Check if point is on the 2D boundary. Returns true if it is
        cnd2D = abs(xmin - fldMsh.nodes(idNode, 1)) < eps ||            ...
                abs(xmax - fldMsh.nodes(idNode, 1)) < eps ||            ...
                abs(ymin - fldMsh.nodes(idNode, 2)) < eps ||            ...
                abs(ymax - fldMsh.nodes(idNode, 2)) < eps;
       
        % Additional check for 3D boundary. Returns true if it is
        cnd3D = false;
        if is3D
            cnd3D = abs(zmin - fldMsh.nodes(idNode, 3)) < eps ||        ...
                    abs(zmax - fldMsh.nodes(idNode, 3)) < eps;          ... 
        end   
         
        % Combine 2D and 3D boundary check
        cnd = cnd2D || cnd3D;
    
        % Update the homDOFs if a node is not on the boundary
        if cnd == false
            bodyNodeCounter = bodyNodeCounter + 1;
            homDOFs(bodyNodeCounter) = idNode;
        end
    end
    
    % Chop unnecessary/redundant entries
    homDOFs  = unique(homDOFs(1:bodyNodeCounter));

end