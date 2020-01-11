function [ F ] = buildRHS( Fin, active_node, gap )
%% Function documentation
%
% BUILDRHS add the gap values to the right hand side
% Add the gap function of every Lagrange multiplier to the load vector
% 
%              Input :
%                Fin : Number of non prescribed DoF
%       active_nodes : List of indices of the nodes for which the matrix should be built
%                      (e.g. the set of all currently active nodes)
%                gap : Vector containing distance of every node to the
%                      segment
%
%             Output :
%                  F : Vector containing in its first entries the force 
%                      on every DoF  and in its last entries   
%                      the gap constant for all nodes defined in the 
%                      list 'active_nodes' and the 
%
%%
N_active=length(active_node);        
F=[Fin ; zeros(N_active,1)];
for i=1:N_active
    F(length(Fin)+i)=-gap(active_node(i),2);
end

end

 