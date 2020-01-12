%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %
%   Authors                                                               %
%   _______                                                               %
%                                                                         %
%   Fabien Pean, Andreas Hauso, Georgios Koroniotis                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [displacement,lagrange] = solveSignoriniLagrange1...
    (mesh,homDBC,contactNodes,F,segmentsPoints,materialProperties,analysis,maxIteration)
%% Function documentation
%
% Returns the displacement field and the Lagrange multipliers corresponding 
% to a plain stress/strain analysis for the given mesh of the geometry 
% together with its Dirichlet and Neumann boundary conditions and the 
% contact constraints for a SINGLE rigid wall by applying the Lagrange multiplier method
% 
%              Input :
%               mesh : Elements and nodes of the mesh
%                 rb : Vector of the prescribed DoFs (by global numbering)
%                 cn : VECTOR containing the global numbering of the canditate-nodes for contact 
%                  F : Global load vector
%     segmentsPoints : Matrix with the coordinates of two wall determining points
% materialProperties : The material properties of the structure
%              graph : Structure containing information on the graphics
%
%             Output :
%       displacement : The resulting displacement field
%           lagrange : The resulting values of the Lagrange multipliers (*.multipliers)
%                      and the node numbers of the active nodes (*.active_nodes)
%
% Function layout :
%
% 0. Remove fully constrained nodes
%
% 1. Compute the gap function
%
% 2. Compute the master stiffness matrix of the structure
%
% 3. Create the expanded system of equations
%
% 4. Solve system iteratively
% |->
%    4.1 Determine the inactive_nodes nodes
%  
%    4.2 Reduce the system of equations according to the constraints
%  
%    4.3 Solve the reduced system of equations
%  
%    4.3 Assemble to the expanded displacement/Lagrange  multiplier vector
% <-|
% 5. Get the values for the displacement and the Lagrange multipliers
%
% 6.  Print info
%
%% Function main body
if strcmp(analysis.type,'planeStress')
    fprintf('Plane stress analysis has been initiated \n');
elseif strcmp(analysis.type,'planeStrain')
    fprintf('Plane strain analysis has been initiated \n');
end
fprintf('\n');


%% 0. Remove fully constrained nodes
nodes.index=contactNodes;
nodes.positions=mesh.nodes(contactNodes,:);
% Remove fully constrained nodes from the tests
for i=length(nodes.index):-1:1
    idx=nodes.index(i);
    if(ismember(2*idx-1:2*idx,homDBC))
        nodes.index(i)=[];
        nodes.positions(i,:)=[];
    end
end
contactNodes=nodes.index;

%% 1. Compute the gap function

% Compute normal, parallel and position vectors of the segments:
segments=buildSegmentsData(segmentsPoints);
% Compute for all nodes the specific gap:
gap=computeGapFunc(nodes, segments, segmentsPoints);

%% 2. Compute the master stiffness matrix of the structure

fprintf('\t Computing master stiffness matrix... \n');

% Master stiffness matrix
K = computeStiffnessMatrixPlateInMembraneActionLinear(mesh,materialProperties,analysis);

%% 3. Create the expanded system of equations

fprintf('\t Creating the expanded system of equations... \n');

% Assemble the values of the normal vector of segments to the constraint matrix:
C=buildConstraintMatrix(length(F),contactNodes,segments);

% Create a zero matrix for the bellow right side of the equation system:
zero_matrix=zeros(length(contactNodes),length(contactNodes));

% Expand the C matrix with the zero matrix:
Ctmp=[C;zero_matrix];

% Expand the K matrix with the C matrix below:
Kexp=[K;C'];

% Expand the K matrix with the Ctmp matrix on the right side:
Kexp=[Kexp,Ctmp];

% Expand the F vector with the gap constants for every node:
Fexp=[F;-gap(:,2)];

clear C;
clear Ctmp;
clear zero_matrix;

% Initial values for the itaration:
it=0;
dexp=zeros(length(Fexp));
inactive_nodes=[];
inactive_old=[];

% Counts the number of equations which are solved during iteration
equations_counter=0;

%% 4. Solve the system iteratively

% Repeat loop until all Lagrange multipliers are valid and system converges 
% to stable solution:
% Check whether solution is found
% If the Lagrange multiplier of the current node are valid and we had 
% no change in the set of inactive_nodes during the last iteration 
% the solution has been found and the loop is terminated
while (it==0 || (it<maxIteration && ~(isequal(inactive_old,inactive_nodes) && max(dexp(length(F)+1:length(Fexp)))<=0)))
 
inactive_old=inactive_nodes;
% update iteration counter
it=it+1;
    
 %% 4.1 Determine the inactive_nodes nodes
 
% Detect non-penetrating nodes and nodes with non-compressive Lagr. multipliers
inactive_nodes= detectInactiveNodes( length(F),contactNodes, dexp, segments, gap );

 
%% 4.2 Reduce the system of equations according to the constraints

% Merge the vectors with equation numbers which will be deleted due to a
% Dirichlet boundary condition or a contact constraint:
rb_cb=[homDBC,inactive_nodes];
rb_cb=unique(rb_cb);
rb_cb=sort(rb_cb);


Kred = Kexp;
Fred = Fexp;

% Remove constrained DOFs and inactive_nodes equations:
for i = length(rb_cb):-1:1
    Kred(:,rb_cb(i)) = [];
    Kred(rb_cb(i),:) = [];
    Fred(rb_cb(i)) = [];
end
 
%% 4.3 Solve the reduced system of equations

equations_counter=equations_counter+length(Kred);
fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(Kred), cond(Kred));

dred=Kred\Fred;


%% 4.3 Assemble to the expanded displacement/Lagrange multiplier vector

% Build dexp out of dred:
dexp = buildFullDisplacement( length(Fexp) ,rb_cb, dred );

 
end

%% 5. Get the values for the displacement and the Lagrange multipliers

% Select and save node numbers of active nodes :
lagrange.active_nodes=setdiff(contactNodes,contactNodes(inactive_nodes-length(F)));

% The first entries of dexp correspond to the displacement
displacement=dexp(1:length(F));

% The last entries of dexp correspond to the Lagrange multipliers:
lagrange.multipliers=dexp(length(F)+1:length(Fexp));

%Keep only Lagr. multipliers of the active nodes:
lagrange.multipliers=lagrange.multipliers(lagrange.multipliers<0);


%% 6. Print info

fprintf('\n');
fprintf('Output informations...\n');
fprintf('\t Constraints solved in %d iterations.A total of %d equations were solved. \n',it,equations_counter);
fprintf('\t %d active nodes found.\n',length(lagrange.active_nodes));
energy=displacement'*K*displacement;
fprintf('\t #DOF: %d .Energy norm of the structure : %4.2f\n',length(F),energy);

end