function  contactNodes = multComputeGapFunc(contactNodes,segments,segmentsPoints )
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
%
for j=1:size(contactNodes,2)
    for i=1:size(contactNodes(j).indices,1)
        %Normal distance to the segment i
        contactNodes(j).gap(i,2)= dot( contactNodes(j).positions(i,1:2) , segments.normals(j,:) ) + segments.constants(j);
        %Parallel distance to the left point of the segment i
        contactNodes(j).gap(i,1)= dot( (contactNodes(j).positions(i,1:2)-segmentsPoints(1,:,j)) , segments.directors(j,:) );
        %Parallel distance to the right point of the segment i
        contactNodes(j).gap(i,3)= dot( (contactNodes(j).positions(i,1:2)-segmentsPoints(2,:,j)) , segments.directors(j,:) );    
    end
end

end