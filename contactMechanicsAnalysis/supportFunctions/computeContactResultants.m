function [contactLength,contactForce,maxContactPressure] = computeContactResultants...
    (mesh,displacement,lagrange,parameters)
% Computes the length of contact, its reaction force and the maximal
% pressure. Used for the comparison with analytical values according to
% Hertz contact theory. 
%
%               Input :
%                mesh : Elements and nodes of the mesh
%        displacement : The displacement field
%        active_nodes : List of indices of the nodes which shall be
%                        evaluated (e.g. the set of all currently active nodes)
%           lagr_mult : Vector with the Lagrange multiplier for every node 
%                       contained in 'active_nodes'
%          parameters : The material properties of the structure
%
%              Output :
%       contactLength : Length of contact area
%        contactForce : The total reaction force on the contact 
%  maxContactPressure : The maximal compressive pressure on the contact 
%
%%

% Initialize the array of the displaced nodes
nodes_displaced = zeros(length(mesh.nodes),3);

% Initialize pseudocounter
counter = 1;

% Update the mesh with displacements
for i=1:length(mesh.nodes)
    % Add the x and y components of the displacement field to mesh nodes
    nodes_displaced(i,1) = mesh.nodes(i,1) + displacement(2*counter-1);
    nodes_displaced(i,2) = mesh.nodes(i,2) + displacement(2*counter);
    
    % Update counter
    counter = counter + 1;
end

k=1;
l=1;
for i=1:length(lagrange.multipliers)
   
    if lagrange.multipliers(i)<0
        
        active_node(k,l) = lagrange.active_nodes(i);
        active_lagr(k,l) = lagrange.multipliers(i);
        k=k+1;
    elseif i>1 && lagrange.multipliers(i-1)<0 && lagrange.multipliers(i)>=0
        l=l+1;
        
    end


end

contactLength(size(active_node,2))=0;
contactForce(size(active_node,2))=0;

for i=1:size(active_node,2)

    contactForce(i)=sum(active_lagr(:,i));
    for j=2:size(active_node,1)
        
        if active_node(j,i)~=0
           x0=nodes_displaced(active_node(j-1,i),1);
           y0=nodes_displaced(active_node(j-1,i),2);
                
           x1=nodes_displaced(active_node(j,i),1);
           y1=nodes_displaced(active_node(j,i),2);
           
           f0=-active_lagr(j-1,i);
           f1=-active_lagr(j,i);
           
           contactLength(i)=contactLength(i)+sqrt((x0-x1)^2+(y0-y1)^2);
          
           contactPressure(j-1,i)=(f0+f1)/(2*parameters.t*(sqrt((x0-x1)^2+(y0-y1)^2)));

        end
    
    end
    
end

maxContactPressure = max(max(contactPressure));

end