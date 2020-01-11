function mesh = meshCADgeometry_Hertz(p,U,q,V,CP,meshAttributes,radius)
%% Function documentation
%
% Returns the elements, the nodes and the boundary keypoints of a mesh
% related to a given surface representation using NURBS-based CAD geometry
%
%              Input :
%                p,q : NURBS polynomial degree
%                U,V : Knot vectors in u,v-direction
%                 CP : set of control point coordinates and weights
%     meshAttributes : Mesh options :
%                      1) boundary : 
%                              uv0 : Meshing parametric line (u,v)=(:,0)
%                              u1v : Meshing parametric line (u,v)=(1,:)
%                              uv1 : Meshing parametric line (u,v)=(:,1)
%                              u0v : Meshing parametric line (u,v)=(0,:)
%                      2) hdata :
%                          hmax : Maximum edge size
%                      3) options : 
%                           maxit : Maximum number of iterations (def. 20)
%                            mlim : The convergence tolerance. The maximum 
%                                   percentage change in edge length per 
%                                   iteration must be less than (def. 0.02)
%                           dhmax : The maximum allowable (relative) gradient 
%                                   in the size function (def.0.3, 30.0%)
%                          output : Displays the mesh and the mesh statistics 
%   `                               upon completion ( def. TRUE)
%             radius : Radius for refinement on circle
%
%             Output :
%               mesh : Elements, nodes, boundary nodes and keypoints of 
%                      the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Get the polygonal description of the CAD boundary
%
% 2. Meshing the interior using the algorithm mesh2d
%
% 3. Get the boundary nodes of the mesh
%
% 4. Add the z-component of the nodes taken by the CAD representation
%
%% Function main body

%% 0. Read input

% Assign a tolerance value
tolerance = 1e-8;

%% 1. Get the polygonal description of the CAD boundary
mesh.keypoints = polygonalizeCADBoundary(p,U,q,V,CP,meshAttributes);

%% 2. Meshing the interior using the algorithm mesh2d
[mesh.nodes,mesh.elements] = mesh2d(mesh.keypoints(:,1:2),[],meshAttributes.hdata,meshAttributes.options);
for j=7:9
M=logical(zeros(length(mesh.elements),1));
for i=1:length(mesh.elements)
    if(mesh.nodes(mesh.elements(i,1),1) > 0.1*j*radius || mesh.nodes(mesh.elements(i,2),1) > 0.1*j*radius || mesh.nodes(mesh.elements(i,3),1) > 0.1*j*radius)
        M(i)=true;
    end
end
[mesh.nodes,mesh.elements]=refine(mesh.nodes,mesh.elements,M);
smoothmesh(mesh.nodes,mesh.elements,meshAttributes.options.maxit*10,tolerance);
end
%% 3. Get the boundary nodes of the mesh

% Find the boundary nodes by logical arrays
[~,nodes_on_boundary_logical] = inpoly(mesh.nodes(:,1:2),mesh.keypoints(:,1:2),[],tolerance);

% Find the boundary nodes by indices
nodes_on_boundary = find(nodes_on_boundary_logical);

% Add the global numbering of the boundary nodes onto the mesh
mesh.boundaryNodes = nodes_on_boundary;

%% 4. Add the z-component of the nodes taken by the CAD representation
for i=1:length(mesh.nodes)
    mesh.nodes(i,3) = mesh.keypoints(1,3);
end

end

