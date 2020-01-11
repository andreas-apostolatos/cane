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
function [displacement, lagrange] = multSolveSignoriniLagrange2(mesh,rb,cn,F, segmentsPoints,materialProperties,analysis,maxIter)
%% Function documentation
%
% Returns the displacement field and the Lagrange multipliers corresponding 
% to a plain stress/strain analysis for the given mesh of the geometry 
% together with its Dirichlet and Neumann boundary conditions and the 
% contact constraints for MULTIPLE rigid walls by applying the Lagrange
% multiplier method.
% In the structure array cn(j) can be specified which canditates for 
% contact nodes are related to a certain wall segment (j)
% 
%              Input :
%               mesh : Elements and nodes of the mesh
%                 rb : Vector of the prescribed DoFs (by global numbering)
%                 cn : STRUCTURE ARRAY 'cn(j=1..n).indices' 
%                      containing the global numbering of the canditate-nodes 
%                      for contact to segment j [segmentPoints(j,:,:)] 
%                      in the field 'indices'
%                  F : Global load vector
%     segmentsPoints : Matrix with the coordinates of two wall determining
%                      points, for every segment j=1..n
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


%% 0. Remove fully constrained nodes
for j=1:size(cn,2)

%Remove fully constrained nodes from the tests
for i=size(cn(j).indices,2):-1:1
    % Determine how many Dirichlet conditions correspond to the node:  
    nodeHasDirichlet=ismember(floor((rb+1)/2),cn(j).indices(i));
    numberOfDirichlet=length(nodeHasDirichlet(nodeHasDirichlet==1));
    % If the 2D node has at least two Dirichlet conditions exclude it from the contact canditates :  
    if (numberOfDirichlet>=2)
       cn(j).indices(i)=[];
    end

end
cn(j).positions=mesh.nodes(cn(j).indices,:);
end

%% 1. Compute the gap function

% Compute normal, parallel and position vectors of the segments:
segments=buildSegmentsData(segmentsPoints);

% Compute for all nodes the specific gap and save it in the field 'cn.gap':
cn=multComputeGapFunc(cn, segments, segmentsPoints);

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
%% 3.3 Assemble to the complete displacement vector

    displacement = buildFullDisplacement(length(F), rb, dred);
    
%% 3.4 Detect active nodes
    active_node_tmp=multDetectActiveNodes(cn,displacement, segments );

    if(isequaln(active_node_tmp,active_node)&& it~=0)
        flag_active=0;
    end
     
%% 3.5 Rebuild system matrices if new active nodes found
    if(~isempty(active_node_tmp) && flag_active)
        %% 3.5.1 Update active nodes
        active_node=active_node_tmp;
        N_active=length(active_node);
        %% 3.5.2 Build constraint matrix C and rhs F
        C=multBuildConstraintMatrix(length(F),cn,active_node,segments);
        Fred=multBuildRHS(F,active_node,cn);
        %% 3.5.3 Build master system matrix
        Kred=[K C ; C' zeros(N_active)];
        %% 3.5.4 Reduce the system according to the BCs
        for i = length(rb):-1:1
            Kred(:,rb(i)) = [];
            Kred(rb(i),:) = [];
            Fred(rb(i)) = [];
        end
    end
%% 3.0 Relax the system until ONLY valid Lagrange multipliers are computed
    flag_lagrange=1;
    while(flag_lagrange && it<maxIter )
        it=it+1;
        flag_lagrange=0;
%% 3.1 Compute the displacement & Lagrange multipliers
        equations_counter=equations_counter+length(Kred);
fprintf('\t Solving the linear system of %d equations, condition number %e ... \n',length(Kred), cond(Kred));
        dred = Kred\Fred;
%% 3.2 Detect and delete non-valid Lagrange multipliers and related rows
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



 allnodes=[];
for j=1:size(cn,2)
    allnodes=[allnodes,cn(j).indices];
end  

lagrange.active_nodes=allnodes(active_node);


end