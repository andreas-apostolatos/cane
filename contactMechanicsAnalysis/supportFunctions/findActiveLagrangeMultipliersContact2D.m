function DOFsLMActive = findActiveLagrangeMultipliersContact2D...
    (mesh, dHat_stiffMtxLM, segmentsContact, propContact)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
% Date : 04.02.2020
%
%% Function documentation
%
% Detects the nodes which are not penetrating the rigid segments and stores
% the IDs of the corresponding Lagrange Multipliers DOFs into an array
%
%              Input :
%             strMsh : Nodes and elements in the mesh
%    dHat_stiffMtxLM : Vector of DOFs containing both displacement and
%                      Lagrange Multipliers
%    segmentsContact : Data sturcture containing information about the 
%                      boundaries of the rigid wall :
%                        .points : a list of 2x2 matrices containing end 
%                                 points of the segment(s)
%                        .number : total number of segments
%                       .normals : normal vector of each segment
%        propContact : Data structure containing the contact properties
%                          .nodeIds : global numbering of contact nodes
%                    .numberOfNodes : number of nodes
%
%             Output :
%     DOFsLMInactive : IDs of the inactive Lagrange Multipliers DOFs
%
%% Function layout :
%
% 0. Read input
%
% 1. Loop over all rigid segments
% ->
%    1i. Loop over all potential contact nodes
%    ->
%        1i.1 Get the DOF numbering of the current node
%
%        1i.2. Get the displacement of the current node
%
%        1i.3. Get the coordinates of the current node
%
%        1i.4. Compute the displaced coordinates of the current node
%
%        1i.5. Get the end vertices of the current rigid segment
%
%        1i.6. Project the current displaced node on the rigid segment
%
%        1i.7. Compute the penetration of the current node with respect to the current segment
%
%        1i.8. Compute non-penetration conditions (displacement)
%
%        1i.9. Check whether the node is active depending on whether any of the above-defined conditions is valid
%
%        1i.10. Update counter
%    <-
% <-
%
% 2. Keep only non-zero entries of the inactive_nodes
%
%% Function main body

%% 0. Read input

% initialize output variable
DOFsLMActive = zeros(1, segmentsContact.number*propContact.numberOfNodes);

% Set contact tolerance
tolerance = sqrt(eps);

% Initialize counters
counterLM = 1;
counterActiveNodes = 1;

%% 1. Loop over all rigid segments
for iSeg = 1:segmentsContact.number
    %% 1i. Loop over all potential contact nodes
    for iCN = 1:propContact.numberOfNodes
        %% 1i.1 Get the DOF numbering of the current node
        DOF = 2*propContact.nodeIDs(iCN) - 1 : 2*propContact.nodeIDs(iCN);
        
        %% 1i.2. Get the displacement of the current node
        u = dHat_stiffMtxLM(DOF);
        
        %% 1i.3. Get the coordinates of the current node
        node = mesh.nodes(propContact.nodeIDs(iCN),1:2);
        
        %% 1i.4. Compute the displaced coordinates of the current node
        nodeDisp = node + u';
        
        %% 1i.5. Get the end vertices of the current rigid segment
        vertexA = segmentsContact.points(1,:,iSeg);
        vertexB = segmentsContact.points(2,:,iSeg);
        
        %% 1i.6. Project the current displaced node on the rigid segment
        lambda = ((vertexA - vertexB)*(nodeDisp - vertexA)')/((vertexA - vertexB)*(vertexB - vertexA)');
        nodeDisp_proj = (1 - lambda)*vertexA + lambda*vertexB;
        
        %% 1i.7. Compute the penetration of the current node with respect to the current segment
        penetration = segmentsContact.normals(iSeg,:)*(nodeDisp - nodeDisp_proj)';
        
        %% 1i.8. Compute non-penetration conditions (displacement)
        isCnd1 = penetration < tolerance;
        isCnd2 = lambda > tolerance;
        isCnd3 = lambda <= 1;
        
        %% 1i.9. Check whether the node is active depending on whether any of the above-defined conditions is valid
        if (isCnd1 && isCnd2 && isCnd3)
            DOFsLMActive(counterActiveNodes) = counterLM;
            counterActiveNodes = counterActiveNodes+1;
        end
        
        %% 1i.10 Update counter
        counterLM = counterLM + 1;
    end
end

%% 2. Keep only non-zero entries of the inactive_nodes
DOFsLMActive = DOFsLMActive(1:counterActiveNodes - 1);

end