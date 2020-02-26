function C = buildConstraintMatrix(nDOF,contactNodes,activeNodes,segments)
%
% Build the constraint matrix to be appended to K
% The constraint matrix is built with these dimensions :
% nDOF x number of active Lagrange multipliers
% The matrix is filled with the normal vectors of the segments applied on
% the right couple of displacement.
%
%             Input :
%              nDOF : Number of DoF of the system
%      contactNodes : structure containing the global numbering of the 
%                     contact canditate nodes for contact to segments
%      active_nodes : List of indices of the nodes for which the matrix 
%                     should be built (e.g. the set of all currently active
%                     nodes).If active_nodes is empty all nodes defined in
%                     contactNodes are taken
%          segments : data stucture containing informations about the
%                     rigid wall segments (field for the normal vector is
%                     required)
%
%            Output :
%                 C : The Constraint matrix
%
%% Function main body

C = zeros(nDOF,1);
k=1;
l=1;
% loop through segments
for j=1:segments.number
    % loop through every contact node in each segment
    for i=1:size(contactNodes.indices,1)
        if isempty(activeNodes) || max(ismember(activeNodes,l))
            % find the index of constrain
            index = 2*contactNodes.indices(i)-1 : 2*contactNodes.indices(i);
            C(index,k)=segments.normals(j,:);
            k=k+1;
        end
        l=l+1;
    end 
end

end 