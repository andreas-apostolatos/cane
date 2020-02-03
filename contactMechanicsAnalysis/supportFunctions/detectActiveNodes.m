function activeNodes = detectActiveNodes...
    (contactNodes,displacement,segments)
%MULTIDETECTACTIVENODES Detect the active nodes for the current segment
% Loops over all segments to check whether the displacement of the node put
% in an active position according to the segment
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
%%
% initialize variables
activeNodes=[];
k=1;

% loop over displacement_exp vector
% loop over the number of contact nodes segments
for j=1:segments.number
    % loop over every node in that segment
    for i=1:size(contactNodes.indices,1)
        
        index = 2*contactNodes.indices(i)-1 :1: 2*contactNodes.indices(i);
        tmp = displacement(index);
        tmp_normal=  dot(tmp,segments.normals(j,:));
        tmp_parallel=dot(tmp,segments.directors(j,:));
        
        % get distances to the segment i
        leftGap = contactNodes.gap(i,1,j);
        normalGap = contactNodes.gap(i,2,j);
        rightGap = contactNodes.gap(i,3,j);
        
        % conditions for geometry (non-penetration)
        cnd1 = tmp_normal + normalGap < sqrt(eps);
        cnd2 = tmp_parallel > min(leftGap,rightGap);
        cnd3 = tmp_parallel < max(leftGap,rightGap);
        
        % if all conditions hold
        if (cnd1 && cnd2 && cnd3)
            activeNodes = [activeNodes,k];
        end
        
        % update counter
        k=k+1;
    end
end

end