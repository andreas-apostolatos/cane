function [ active_nodes ] = multDetectActiveNodes( cn, displacement, segments )
%MULTIDETECTACTIVENODES Detect the active nodes for the current segment
% Loops over all segments to check whether the displacement of the node put
% in an active position according to the segment
% APPLY THIS FUNCTION IF MORE THAN ONE RIGID WALL SEGMENT EXISTS
%
%              Input :
%                 cn : STRUCTURE ARRAY 'cn(j=1..n).indices' 
%                      containing the global numbering of the canditate-nodes 
%                      for contact to segments(j) 
%                      in the field 'indices' and the gap in the field 'gap'
%       displacement : Vector of the FULL displacement field
%           segments : Structure containing infos about constraint segments
%      
%             Output :
%        active_node : The resulting vector containing index of the 
%         restricted vector of mesh.boundaryNodes containing global indices
%
active_nodes=[];
k=1;
for j=1:size(cn,2)
for i=1:size(cn(j).indices,2)%loop over every node
    tmp_normal=dot(displacement(2*cn(j).indices(i)-1:1:2*cn(j).indices(i)),segments.normals(j,:));%computation for 2D analysis with 2disp DoFs per node
    tmp_parallel=dot(displacement(2*cn(j).indices(i)-1:1:2*cn(j).indices(i)),segments.directors(j,:));

    if (tmp_normal+cn(j).gap(i,2)<sqrt(eps) && tmp_parallel>min(cn(j).gap(i,1),cn(j).gap(i,3)) && tmp_parallel<max(cn(j).gap(i,1),cn(j).gap(i,3)))
     active_nodes=[active_nodes,k];
    end
    k=k+1;
end
end
end

 