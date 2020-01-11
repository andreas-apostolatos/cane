function [ inactive_nodes ] = detectInactiveNodes(DOF, cn, dexp, segments, gap )
%DETECTINACTIVENODES Detect the inactive nodes for the current segment
% Loops over all segments to check whether the displacement of the node put
% in an active position according to the segment, or whether the Lagrange
% multipliers are non-compressevive
%              Input :
%                DOF : Number of DoFs in the system
%                 cn : List of indices of the nodes which shall be
%                       evaluated (e.g. all nodes on the boundary)
%               dexp : Vector of the FULL displacement field and Lagrange multipliers
%           segments : Structure containing infos about constraint segments
%                gap : Vector of distances of the nodes to the segments
%      
%             Output :
%      inactive_node : The resulting vector containing index of the 
%         restricted vector of mesh.boundaryNodes containing global indices
%
inactive_nodes=[];
for i=1:length(cn)%loop over every node
    tmp_normal=dot(dexp(2*cn(i)-1:1:2*cn(i)),segments.normals(1,:));%computation for 2D analysis with 2disp DoFs per node
    tmp_parallel=dot(dexp(2*cn(i)-1:1:2*cn(i)),segments.directors(1,:));

    if (tmp_normal+gap(i,2)>sqrt(eps)|| tmp_parallel<min(gap(i,1),gap(i,3)) || tmp_parallel>max(gap(i,1),gap(i,3)) )
     inactive_nodes=[inactive_nodes,DOF+i];
    elseif (dexp(DOF+i)>0)
     inactive_nodes=[inactive_nodes,DOF+i];
    end
end

end 