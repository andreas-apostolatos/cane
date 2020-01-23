function inactive_nodes = multDetectInactiveNodes...
    (nDOF,contactNodes,displacement_exp,segments)
%MULTIDETECTINACTIVENODES Detect the inactive nodes for the current segment
% Loops over all segments to check whether the displacement of the node put
% in an active position according to the segment, or whether the Lagrange
% multipliers are non-compressevive
% APPLY THIS FUNCTION IF MORE THAN ONE RIGID WALL SEGMENT EXISTS
%
%              Input :
%               nDOF : total number of DOFs in the system
%       contactNodes : List of indices of the nodes which shall be evaluated
%   displacement_exp : Vector of the FULL displacement field and Lagrange multipliers
%           segments : Structure containing infos about constraint segments
%      
%             Output :
%      inactive_node : The resulting vector containing index of the 
%         restricted vector of mesh.boundaryNodes containing global indices

% initialize variables
inactive_nodes=[];
k=1;
% loop over the number of contact nodes segments
for j=1:size(contactNodes,2)
    % loop over every node in that segment
    for i=1:size(contactNodes(j).indices,1)
        
        tmp = displacement_exp(2*contactNodes(j).indices(i)-1:1:2*contactNodes(j).indices(i));
        tmp_normal = dot(tmp,segments.normals(j,:));
        tmp_parallel = dot(tmp,segments.directors(j,:));

        % get distances to the segment i
        leftGap = contactNodes(j).gap(i,1);
        normalGap = contactNodes(j).gap(i,2);
        rightGap = contactNodes(j).gap(i,3);
        
        % conditions for geometry (non-penetration)
        cnd1 = tmp_normal + normalGap > sqrt(eps);
        cnd2 = tmp_parallel < min(leftGap,rightGap);
        cnd3 = tmp_parallel > max(leftGap,rightGap);
        
        % condition for Lagrange multipliers (non-compressive)
        cnd4 = displacement_exp(nDOF+i) > 0;
        
        % if any of the conditions hold
        if (cnd1 || cnd2 || cnd3 || cnd4)
            inactive_nodes = [inactive_nodes,nDOF+k];
        end
        
        % update counter
        k=k+1;
    end
end

end