function  contactNodes = computeGapFunction(contactNodes,segments)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% Expands a data structure contactNodes. Gap function is the distance
% between nodes of the structure and the segments points (walls)
%
%              Input :
%       contactNodes :
%           .indices : global numbering of the contact canditate-nodes
%         .positions : coordinates of the candidate nodes
%
%           segments : data stucture containing informations about the
%                      rigid wall segments (normal vector, parallel vector,
%                      position)
%            .points : points(:,:,i) is a list of 2x2 matrices containing
%                      end points of the segment(s)
%
%             Output :
%   contactNodes.gap : Matrix containing distance of every node to each
%                      segment (normal and parallel to both ends)
%
%% %% Function main body

% loop through the number of segments
for j=1:segments.number
    % loop through the contact nodes
    for i=1:size(contactNodes.indices,1)
        % get node of interest
        nodeOfInterest = contactNodes.positions(i,1:2);
        
        % Normal distance to the segment i
        contactNodes.gap(i,2,j)= dot( nodeOfInterest , segments.normals(j,:) ) + segments.constants(j);
        
        % Parallel distance to the left point of the segment i
        contactNodes.gap(i,1,j)= dot( (nodeOfInterest-segments.points(1,:,j)) , segments.directors(j,:) );
        
        % Parallel distance to the right point of the segment i
        contactNodes.gap(i,3,j)= dot( (nodeOfInterest-segments.points(2,:,j)) , segments.directors(j,:) );    
    end
end

end