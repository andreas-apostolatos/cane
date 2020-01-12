function [ gap ] = computeGapFunction( nodes , segments, points )
%COMPUTEGAPFUNC Computes gap function
% Returns the gap function which is the distance between nodes of the
% structure and the segments points.
%              Input :
%              nodes : Structure containing infos about selected nodes
%                      (e.g. all nodes on the boundaries), especially
%                      index and position 
%           segments : data stucture containing informations about the
%                      rigid wall segments (normal vector, parallel vector,
%                      position)
%             points : points(i,:,:) is a list of 2x2 matrices containing
%                      extremities points of the segment(s)
%      
%             Output :
%                gap : Matrix containing distance of every node to the
%                      segment (normal, and parallel to both
%                      extremities of the segment)
%
gap=zeros(length(nodes.index),3);
for i=1:length(nodes.index)
    %Normal distance to the segment i
    gap(i,2)=dot(nodes.positions(i,1:2),segments.normals(1,:)) + segments.constants(1);
    %Parallel distance to the left point of the segment i
    gap(i,1)=dot( (nodes.positions(i,1:2)'-squeeze(points(1,1,:))) , segments.directors(1,:));
    %Parallel distance to the right point of the segment i
    gap(i,3)=dot( (nodes.positions(i,1:2)'-squeeze(points(1,2,:))) , segments.directors(1,:));    
end
end

 