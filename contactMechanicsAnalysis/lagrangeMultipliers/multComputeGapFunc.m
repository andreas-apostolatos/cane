function  contactNodes = multComputeGapFunc(contactNodes,segments)
%MULTICOMPUTEGAPFUNC Computes gap function
% Returns the gap function which is the distance between nodes of the
% structure and the segments points.
% APPLY THIS FUNCTION IF MORE THAN ONE RIGID WALL SEGMENT EXISTS
%
%              Input :
%       contactNodes : STRUCTURE ARRAY 'cn(j=1..n).indices' 
%                      containing the global numbering of the canditate-nodes 
%                      for contact to segments(j) 
%                      in the field 'indices' and the position in the field 'position'
%           segments : data stucture containing informations about the
%                      rigid wall segments (normal vector, parallel vector,
%                      position)
%             points : points(:,:,i) is a list of 2x2 matrices containing
%                      extremities points of the segment(s)
%      
%             Output :
%                gap : Matrix containing distance of every node to the
%                      segment (normal, and parallel to both
%                      extremities of the segment)
%%
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