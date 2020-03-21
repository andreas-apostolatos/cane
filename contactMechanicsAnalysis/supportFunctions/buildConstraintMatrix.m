function C = buildConstraintMatrix(nDOF,propContact,activeNodes,segments)
%
% Build the constraint matrix to be appended to K. The constraint matrix is
% built with dimensions: nDOF x number of active Lagrange multipliers. The
% matrix is filled with the normal vectors of the segments applied on the
% right couple of displacement
%
%             Input :
%              nDOF : Number of DoF of the system
%       propContact : structure containing the global numbering of the 
%                     contact canditate nodes for contact to segments
%      active_nodes : List of indices of the nodes for which the matrix 
%                     should be built (e.g. the set of all currently active
%                     nodes).If activeNodes is empty all nodes defined in
%                     contactNodes are taken
%          segments : data stucture containing informations about the
%                     rigid wall segments (field for the normal vector is
%                     required)
%
%            Output :
%                 C : The Constraint matrix
%
%% Function main body

% intialize C matrix
C = zeros(nDOF,1);

k=1;
l=1;
% loop through segments
for m=1:segments.number
    % loop through every contact node in each segment
    for n=1:propContact.numberOfNodes
        if isempty(activeNodes) || max(ismember(activeNodes,l))
            % find the index of constrain
            DOF = 2*propContact.nodeIDs(n)-1 : 2*propContact.nodeIDs(n);
            C(DOF,k) = segments.normals(m,:);
            k=k+1;
        end
        l=l+1;
    end 
end

end 