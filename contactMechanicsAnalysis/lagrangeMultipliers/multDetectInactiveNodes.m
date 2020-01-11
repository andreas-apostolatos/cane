function inactive_nodes = multDetectInactiveNodes(DOF, cn, dexp, segments )
%MULTIDETECTINACTIVENODES Detect the inactive nodes for the current segment
% Loops over all segments to check whether the displacement of the node put
% in an active position according to the segment, or whether the Lagrange
% multipliers are non-compressevive
% APPLY THIS FUNCTION IF MORE THAN ONE RIGID WALL SEGMENT EXISTS
%
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
    inactive_nodes=[];
k=1;
for j=1:size(cn,2)
for i=1:size(cn(j).indices,2)%loop over every node
    tmp_normal=dot(dexp(2*cn(j).indices(i)-1:1:2*cn(j).indices(i)),segments.normals(j,:));%computation for 2D analysis with 2disp DoFs per node
    tmp_parallel=dot(dexp(2*cn(j).indices(i)-1:1:2*cn(j).indices(i)),segments.directors(j,:));

    if (tmp_normal+cn(j).gap(i,2)>sqrt(eps) || tmp_parallel<min(cn(j).gap(i,1),cn(j).gap(i,3)) || tmp_parallel>max(cn(j).gap(i,1),cn(j).gap(i,3)))
     inactive_nodes=[inactive_nodes,DOF+k];
    elseif (dexp(DOF+i)>0)
     inactive_nodes=[inactive_nodes,DOF+k];
    end
    k=k+1;
end
end
end

 