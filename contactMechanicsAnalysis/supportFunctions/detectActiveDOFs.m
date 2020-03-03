function activeNodes = detectActiveDOFs(contactNodes,displacement,segments)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
% Date : 04.02.2020
%
%% Function documentation
%
% Detect the active nodes for the current segment. Loops over all segments
% to check whether the displacement of the node put in an active position
% according to the segment
%
%          Input :
%   contactNodes : structure containing the global numbering of contact
%                  canditate-nodes coordinates of the candidate nodes
%   displacement : Vector of the FULL displacement field
%       segments : data stucture containing informations about the
%                  rigid wall segments (normal vector, parallel vector,
%                  position)
%      
%         Output :
%    active_node : Vector containing global indices of active contact nodes
%
%% Function main body

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