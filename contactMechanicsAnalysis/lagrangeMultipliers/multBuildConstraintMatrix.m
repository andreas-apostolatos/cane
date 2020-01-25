function C = multBuildConstraintMatrix(DOF,contactNodes,activeNodes,segments)
%
% MULTBUILDCONSTRAINTMATRIX Build the constraint matrix to be appended to K
% The constraint matrix is built with these dimensions :
% number of displacement unknowns X number of Lagractive_nodesge multipliers
% The matrix is filled with the normal vectors of the segments applied on
% the right couple of displacement.
% APPLY THIS FUNCTION IF MORE THAN ONE RIGID WALL SEGMENT EXISTS
%
%              Input :
%                DOF : Number of DoF of the system
%       contactNodes : STRUCTURE ARRAY 'cn(j=1..n).indices' 
%                      containing the global numbering of the canditate-nodes 
%                      for contact to segments(j) 
%                      in the field 'indices'
%       active_nodes : List of indices of the nodes for which the matrix should be built
%                      (e.g. the set of all currently active nodes).If
%                      active_nodes is empty all nodes defined in cn are
%                      taken
%           segments : data stucture containing informations about the
%                      rigid wall segments (field for the normal vector is
%                      required)
%
%             Output :
%                  C : The Constraint matrix
%
%%

C = zeros(DOF,1);
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