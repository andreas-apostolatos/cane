function F_exp = multBuildRHS(F,contactNodes,activeNodes,segments)
%% Function documentation
% BUILDRHS add the gap values to the right hand side
% Add the gap function of every Lagrange multiplier to the load vector
% APPLY THIS FUNCTION IF MORE THAN ONE RIGID WALL SEGMENT EXISTS
% 
%              Input :
%                Fin : Number of non prescribed DoF
%       active_nodes : List of indices of the nodes for which the matrix should be built
%                      (e.g. the set of all currently active nodes)
%                 cn : STRUCTURE ARRAY 'cn(j=1..n).indices' 
%                      containing the global numbering of the canditate-nodes 
%                      for contact to segments(j) 
%                      in the field 'indices' and the gap in the field 'gap'
%
%             Output :
%                  F : Vector containing in its first entries the force 
%                      on every DoF  and in its last entries   
%                      the gap constant for all nodes defined in the 
%                      list 'active_nodes' and the 
%
%%
% get number of DOFs
nDOF = length(F);

% initialize F_exp and add zeros
zero_vector = zeros(size(activeNodes,2),1);
F_exp = [F;zero_vector];

% initialize counters
k=1;
l=1;
% loop through the number of segments
for j=1:segments.number
    % loop through contact nodes
    for i=1:size(contactNodes.indices,1)

        if (isempty(activeNodes) || max(ismember(activeNodes,l)))
            F_exp(nDOF+k) = -contactNodes.gap(i,2,j);
            
            k=k+1;
        end
        l=l+1;
    end
end

end

 