function  propContact = computeGapFunction(mesh,propContact,segmentsContact)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
% Date: 07.03.2020
%
%% Function documentation
%
% Expands the data structure propContact. The gap function is the distance
% between the nodes of the deformable structure and the points on the 
% segments defining the boundary of the rigid wall
%
%             Input :
%              mesh : Finite element mesh of the deformable body
%                       .nodes : Nodes of the mesh
%                    .elements : Elements of the mesh
%       propContact : Data structure containing the contact properties,
%                           .nodeIds : global numbering of contact nodes
%                     .numberOfNodes : number of nodes 
%   segmentsContact : Data sturcture containing information about the 
%                     boundaries of the rigid wall :
%                           .points : a list of 2x2 matrices containing
%                                      end points of the segment(s)
%                            .number : total number of segments
%                           .normals : normal vector of each segment
%            Output :
%       propContact : Updated data structure of the contact properties with
%                      .gap : Normal distance of every node to each segment
%
% Function layout :
%
% 1. Loop over all segments
% ->
%    1i. Loop over all contact nodes
%    ->
%        1i.1. Get the nodal coordinates
%
%        1i.2. Get the vertices of the segment
%
%        1i.3. Project the node on the segment
%
%        1i.4. Compute the normal distance of the node to the segment
%    <-
% <-
%
%% Function main body

%% 1. Loop over all segments
for iSeg = 1:segmentsContact.number
    %% 1i. Loop over all contact nodes
    for iCN = 1:size(propContact.nodeIDs,1)
        %% 1i.1. Get the nodal coordinates
        nodeContact = mesh.nodes(propContact.nodeIDs(iCN),1:2);
        
        %% 1i.2. Get the vertices of the segment
        vertexA = segmentsContact.points(1,:,iSeg);
        vertexB = segmentsContact.points(2,:,iSeg);
        
        %% 1i.3. Project the node on the segment
        lambda = ((vertexA - vertexB)*(nodeContact-vertexA)')/((vertexA - vertexB)*(vertexB - vertexA)');
        nodeContact_proj = (1 - lambda)*vertexA + lambda*vertexB;
        
        %% 1i.4. Compute the normal distance of the node to the segment
        propContact.gap(iCN,iSeg) = segmentsContact.normals(iSeg,:)*(nodeContact - nodeContact_proj)';
    end
end

end