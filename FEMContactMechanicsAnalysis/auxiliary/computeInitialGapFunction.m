function  gapFunctionInitial = computeInitialGapFunction ...
    (mesh, propContact, segmentsContact)
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
% Returns the initial nodal gap function, that is the normal distance of 
% each potential contact node to all the rigid segments.
%
%              Input :
%               mesh : Finite element mesh of the deformable body
%                        .nodes : Nodes of the mesh
%                     .elements : Elements of the mesh
%        propContact : Data structure containing the contact properties,
%                           .nodeIds : global numbering of contact nodes
%                          .numNodes : number of nodes 
%    segmentsContact : Data sturcture containing information about the 
%                      boundaries of the rigid wall :
%                            .points : a list of 2x2 matrices containing
%                                       end points of the rigid segments
%                       .numSegments : total number of segments
%                           .normals : normal vector of each segment
%             Output :
% gapFunctionInitial : Normal distance of every node to each segment
%
% Function layout :
%
% 0. Initialize gap function for each node to segment combination
%
% 1. Loop over all rigid segments
% ->
%   1i. Loop over all contact nodes
%   ->
%       1i.1 Get the id of the corresponding Lagrange Multipliers DOF
%
%       1i.2. Get the nodal coordinates
%
%       1i.3. Get the vertices of the segment
%
%       1i.4. Project the node on the segment
%
%       1i.5. Compute the normal distance of the node to the segment
%   <-
% <-
%
%% Function main body

%% 0. Initialize gap function for each node to segment combination
noDOFsLM = propContact.numNodes*segmentsContact.numSegments;
gapFunctionInitial = zeros(noDOFsLM, 1);

%% 1. Loop over all rigid segments
for iSeg = 1:segmentsContact.numSegments
    %% 1i. Loop over all contact nodes
    for iNode = 1:propContact.numNodes
        %% 1i.1 Get the id of the corresponding Lagrange Multipliers DOF
        idLM = (iSeg - 1)*propContact.numNodes + iNode;
        
        %% 1i.2. Get the nodal coordinates
        nodeContact = mesh.nodes(propContact.nodeIDs(iNode), 2:3);
        
        %% 1i.3. Get the vertices of the segment
        vertexA = segmentsContact.points(iSeg, 1:2);
        vertexB = segmentsContact.points(iSeg, 3:4);
        
        %% 1i.4. Project the node on the segment
        [nodeContact_proj, ~] = computePointProjectionOnLine2D...
            (nodeContact, vertexA, vertexB);
        
        %% 1i.5. Compute the normal distance of the node to the segment
        gapFunctionInitial(idLM, 1) = ...
            - segmentsContact.normals(iSeg,:)*(nodeContact - nodeContact_proj)';
    end
end

end
