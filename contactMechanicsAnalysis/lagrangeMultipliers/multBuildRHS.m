function [ F ] = multBuildRHS( F,active_nodes, cn )
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
DOF=length(F);

N_active_node=0;
for j=1:size(cn,2)
    N_active_node=N_active_node+size(cn(j).indices,2);       
end
F(N_active_node)=0;

k=1;
l=1;
for j=1:size(cn,2)
    N_node=size(cn(j).indices,2);
for i=1:N_node

    if isempty(active_nodes) || max(ismember(active_nodes,l))
    F(DOF+k)=squeeze(-cn(j).gap(i,2));
    
    k=k+1;
    end
    l=l+1;
end
end
end

 