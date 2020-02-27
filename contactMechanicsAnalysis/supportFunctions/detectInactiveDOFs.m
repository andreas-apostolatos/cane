function inactive_DOFs = detectInactiveDOFs...
    (nDOF,contactNodes,displacement_exp,segments)
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
% Detect the inactive nodes for the current segment Loops over all segments
% to check whether the displacement of the node put in an active position
% according to the segment, or whether the Lagrange multipliers are valid
%
%              Input :
%               nDOF : total number of DOFs in the system
%       contactNodes : structure containing the global numbering of contact
%                      canditate-nodes coordinates of the candidate nodes
%   displacement_exp : Vector of the FULL displacement and Lagrange multipliers
%           segments : data stucture containing informations about the
%                      rigid wall segments (normal vector, parallel vector,
%                      position)
%      
%             Output :
%      inactive_DOFs : The resulting vector containing index of the 
%                      restricted vector of mesh.boundaryNodes containing 
%                      the global indices
%
%% Function main body

% initialize variables
inactive_DOFs=[];
k=1;

% loop over displacement_exp vector
% loop over the number of contact nodes segments
for j=1:segments.number
    % loop over every node in that segment
    for i=1:size(contactNodes.indices,1)
        
        index = 2*contactNodes.indices(i)-1 :1: 2*contactNodes.indices(i);
        tmp_displacement = displacement_exp(index);
        tmp_normal = dot(tmp_displacement,segments.normals(j,:));
        tmp_parallel = dot(tmp_displacement,segments.directors(j,:));

        % get distances to the segment i
        leftGap = contactNodes.gap(i,1,j);
        normalGap = contactNodes.gap(i,2,j);
        rightGap = contactNodes.gap(i,3,j);
        
        % conditions for geometry (non-penetration)
        isCnd1 = tmp_normal + normalGap > sqrt(eps);
        isCnd2 = tmp_parallel > leftGap;
        isCnd3 = tmp_parallel < rightGap;
        
        % condition for Lagrange multipliers (non-compressive)
        isCnd4 = displacement_exp(nDOF+k) > 0;
        
        % if any of the conditions hold then the node is inactive
        if (isCnd1 || isCnd2 || isCnd3 || isCnd4)
            inactive_DOFs = [inactive_DOFs,nDOF+k];
        end
        
        % update counter
        k=k+1;
    end
end

end