function cn = getContactNodes(p,U,q,V,CP,boundaryConditions,mesh)
%% Function documentation
%
% Selectes nodes on the contact boundary
%
%              Input :
%                p,q : NURBS polynomial degree
%                U,V : Knot vectors in u,v-direction
%                 CP : set of control point coordinates and weights
% boundaryConditions : contains information on the selected contact boundaries:
%                             increment : Increment on the search algorithm
%                      tolerance_points : Tolerance on the search algorithm
%                
%               mesh : Contains nodes and elements of the mesh
%
%             Output :
%                 cn : Vector 1xn containing the numbering of the nodes
%                      on the prescribed boundaries
%
%% Function main body 

% Number of Control Points in u-direction
nu = length(CP(:,1,1));

% Number of Control Points in v-direction
nv = length(CP(1,:,1));

% Initialize counter for the nodes on the boundary
node_counter = 1;

% Loop over all boundary segments with a contact contraint
for i=1:boundaryConditions.contact.number_of_segments
    
    % Get the parametric domain of the segment i
    segmentu = boundaryConditions.contact.segmentu(i,:);
    segmentv = boundaryConditions.contact.segmentv(i,:);
    
    if isscalar(segmentu) || (~isscalar(segmentv) && boundaryConditions.contact.isOnU(i)==0)
        % Set up the starting coordinates
        u = segmentu(1);
        v = segmentv(1);
        
        % Fixed coordinate
        is_on_u = 0;
        
        % Running coordinate
        s = v;
        S = segmentv;
        
    elseif isscalar(segmentv) || (~isscalar(segmentu) && boundaryConditions.contact.isOnU(i)==1)
        % Set up the starting coordinates
        u = segmentu(1);
        v = segmentv(1);
        
        % Fixed coordinate
        is_on_u = 1;
        
        % Running coordinate
        s = u;
        S = segmentu;
    end
    
    % Initialize contact point counter
    contact_point_counter = 1;
    
    % Loop over all the test points on the contact boundary domain
    while s<=S(length(S))+boundaryConditions.tolerance_search
        
        % Find the span index in u-direction
        spanu = findspan(u,U,nu);

        % Find the span index in u-direction
        spanv = findspan(v,V,nv);
        
        % Get the point on contact boundary domain
        point_on_contactBound(contact_point_counter,:) = point_on_surface(p,spanu,u,U,q,spanv,v,V,CP);
        
        % update the parametric location
        s = s + boundaryConditions.increment;
        
        if is_on_u
            u = s;
        else
            v = s;
        end
        
        % Update counter
        contact_point_counter = contact_point_counter + 1;
    end
    
    % Find the contact boundary node numbering
    for j=1:length(mesh.boundaryNodes) 
        for k=1:length(point_on_contactBound)
            
            % Select node on the boundary 
            node_on_mesh_boundary = mesh.nodes(mesh.boundaryNodes(j),:);
            point_on_boundary = point_on_contactBound(k,:);
            
            % Get their distance in the L2-norm
            L2normDistance = norm(node_on_mesh_boundary-point_on_boundary);
            
            if L2normDistance<=boundaryConditions.tolerance_points
                % Assign the corresponding node number
                    cn(node_counter) =mesh.boundaryNodes(j);
                    
                    % Update node counter
                    node_counter = node_counter+1;    
                
            end
            
        end        
    end
    
    % Delete dublicate elements from the selected node number array
    cn = unique(cn);
    
    % Adjust accordingly the counter
    node_counter = length(cn)+1;
end


end

