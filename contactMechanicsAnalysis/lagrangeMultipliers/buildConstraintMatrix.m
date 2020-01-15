function C = buildConstraintMatrix( N_unknown , active_nodes, segments )
%% Function documentation
%
% BUILDCONSTRAINTMATRIX Build the constraint matrix to be appended to K
% The constraint matrix is built with these dimensions :
% number of displacement unknowns X number of Lagractive_nodesge multipliers
% The matrix is filled with the normal vectors of the segments applied on
% the right couple of displacement.
% 
%              Input :
%          N_unknown : Number of non prescribed DoF
%       active_nodes : List of indices of the nodes for which the matrix should be built
%                      (e.g. the set of all currently active nodes)
%           segments : data stucture containing informations about the
%                      rigid wall segments (field for the normal vector is
%                      required)
%
%             Output :
%                  C : The Constraint matrix
%
%%

N_active_nodes=length(active_nodes);
C=zeros(N_unknown,N_active_nodes);
for i=1:N_active_nodes
    C(2*active_nodes(i)-1:2*active_nodes(i),i)=segments.normals(1,:);
end
       
end

 