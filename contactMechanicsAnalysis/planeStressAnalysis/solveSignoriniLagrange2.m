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
function [displacement,lagrange] = solveSignoriniLagrange2(mesh,rb,cn,F,segmentsPoints,materialProperties,analysis,maxIter)
%% Function documentation
% 
% Returns the displacement field and the Lagrange multipliers corresponding 
% to a plain stress/strain analysis for the given mesh of the geometry 
% together with its Dirichlet and Neumann boundary conditions and the 
% contact constraints for a single rigid wall by applying the Lagrange multiplier method
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
% 3. Reduce the system according to the given constraints
%|->
% 4. Solve the system according to Lagrange multipliers
%   4.1 Assemble to the complete displacement vector
%   4.2 Detect active nodes
%   4.3 Rebuild system if new active nodes found
%       (4.3.1 Update number of active nodes
%       (4.3.2 Build constraint matrix C and rhs F
%       (4.3.3 Build master system matrix
%       (4.3.4 Reduce the system according to the BCs
%   4.4 Relax the system until ONLY valid Lagrange multipliers are computed
%   |->
%   4.4.1 Compute the displacement & Lagrange multipliers
%   4.4.2 Detect and delete non-valid Lagrange multipliers and related rows
%   <-|
%<-|
% 5. compute the complete load vector and verify the results
%
%% Function main body
if strcmp(analysis.physics,'plain_stress')
    fprintf('Plain stress analysis has been initiated \n');
elseif strcmp(analysis.physics,'plain_strain')
    fprintf('Plain strain analysis has been initiated \n');
end
fprintf('\n');


%%  0. Remove fully constrained nodes
nodes.index=cn;
nodes.positions=mesh.nodes(cn,:);
%Remove fully constrained nodes from the tests
for i=length(nodes.index):-1:1
    idx=nodes.index(i);
    if(ismember(2*idx-1:2*idx,rb))
        nodes.index(i)=[];
        nodes.positions(i,:)=[];
    end
end
cn=nodes.index;

%% 1. Compute the gap function

% Compute normal, parallel and position vectors of the segments:
segments=buildSegmentsData(segmentsPoints);
% Compute for all nodes the specific gap:
gap=computeGapFunc(nodes, segments, segmentsPoints);

%% 2. Compute the master stiffness matrix of the structure

fprintf('\t Computing master stiffness matrix... \n');

% Master stiffness matrix
K = computeStiffnessMatrixPlateMembraneActionLinear(mesh,materialProperties,analysis);

%% 3. Reduce the system according to the given constraints
fprintf('\t Reducing the system according to the constraints... \n');

Kred = K;
Fred = F;

for i = length(rb):-1:1
    Kred(:,rb(i)) = [];
    Kred(rb(i),:) = [];
    Fred(rb(i)) = [];
end
%% 4. Solve the system according to Lagrange multipliers

% Initial values for the itaration:
    dred=zeros(length(Fred),1);
    flag_active=1;
    N_active=0;
    active_node=[];
    it=0;
    
% Counts the number of equations which are solved during iteration    
    equations_counter=0;
    
% Iterate until no more invalid Lagrange multipliers AND no new active node added in the pool
while( flag_active && it<maxIter )    
%% 4.1 Assemble to the complete displacement vector

    displacement = buildFullDisplacement(length(F), rb, dred);
    
%% 4.2 Detect active nodes
    active_node_tmp=detectActiveNodes(nodes,displacement,segments,gap,0);
%     flag_active=0;
%     if(~isempty(setxor(active_node_tmp,active_node)))
%         flag_active=1;
%     end
    if(isequaln(active_node_tmp,active_node)&& it~=0)
        flag_active=0;
    end
     
%% 4.3 Rebuild system matrices if new active nodes found
    if(~isempty(active_node_tmp) && flag_active)
        %% 4.3.1 Update active nodes
        active_node=active_node_tmp;
        N_active=length(active_node);
        %% 4.3.2 Build constraint matrix C and rhs F
        C=buildConstraintMatrix(length(F),cn(active_node),segments);
        Fred=buildRHS(F,active_node,gap);
        %% 4.3.3 Build master system matrix
        Kred=[K C ; C' zeros(N_active)];
        %% 4.3.4 Reduce the system according to the BCs
        for i = length(rb):-1:1
            Kred(:,rb(i)) = [];
            Kred(rb(i),:) = [];
            Fred(rb(i)) = [];
        end
    end
%% 4.4 Relax the system until ONLY valid Lagrange multipliers are computed
    flag_lagrange=1;
    while(flag_lagrange && it<maxIter )
        it=it+1;
        flag_lagrange=0;
        %% 4.4.1 Compute the displacement & Lagrange multipliers
        equations_counter=equations_counter+length(Kred);
        fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(Kred), cond(Kred));
        dred = Kred\Fred;
        %% 4.4.2 Detect and delete non-valid Lagrange multipliers and related rows
        for i=length(dred):-1:(length(dred)-N_active+1)
            if(dred(i)>=0)
                Kred(:,i) = [];
                Kred(i,:) = [];
                Fred(i) = [];
                active_node(i-length(dred)+N_active)=[];
                flag_lagrange=1;
                flag_active=1;
            end
        end
        N_active=length(active_node);
    end
    
end
displacement = buildFullDisplacement(length(F), rb, dred);
lagrange.multipliers=dred(length(dred)-N_active+1:length(dred));
fprintf('\n');
fprintf('Output informations...\n');
fprintf('\t Constraints solved in %d iterations.A total of %d equations were solved. \n',it,equations_counter);
fprintf('\t %d active nodes found.\n',N_active);
energy=displacement'*K*displacement;
fprintf('\t #DOF: %d .Energy norm of the structure : %4.2f\n',length(F),energy);


%% Build back the active nodes indexes into mesh.boundaryNodes compatible
offset = 0;
j=1;
for i=1:length((cn))
    if(j>length(active_node))
        break;
    end
    idx=cn(i);
    if(ismember(2*idx-1:2*idx,rb))
        offset=offset+1;
    end
    if(idx==nodes.index(active_node(j)))
        active_node(j)=active_node(j)+offset;
        j=j+1;
    end
end

lagrange.active_nodes=cn(active_node);







end