function [ active_node ] = detectActiveNodes( nodes, displacement, segments, gap, nonlateral )
%DETECTACTIVENODES Detect the active nodes for the current segment
% Loops over all segments to check whether the displacement of the node put
% in an active position according to the segment´
%
%              Input :
%              nodes : Structure containing infos about nodes on boundaries
%       displacement : Vector of the FULL displacement field
%           segments : Structure containing infos about constraint segments
%                gap : Vector of distances of the nodes to the segments
%         nonlateral : boolean value allowing lateral check of the segments
%      
%             Output :
%        active_node : The resulting vector containing index of the 
%         restricted vector of mesh.boundaryNodes containing global indices
%
if(isempty(nonlateral))
    nonlateral=0;
end
j=0;
active_node=[];
for i=1:length(nodes.index)%loop over every node
    idx=nodes.index(i);
    tmp=dot(displacement(2*idx-1:1:2*idx),segments.normals(1,:));%computation for 2D analysis with 2disp DoFs per node
    if (tmp+gap(i,2)<sqrt(eps))%sqrt(eps) is a little hack to avoid rounding error to cause problems of convergence
        tmp=dot(displacement(2*idx-1:1:2*idx),segments.directors(1,:));
        if(nonlateral || (tmp+gap(i,1)<0 && tmp+gap(i,3)>0) || (tmp+gap(i,3)<0 && tmp+gap(i,1)>0))
        j=j+1;
        active_node(j)=i;
        end
    end
end

end
 
