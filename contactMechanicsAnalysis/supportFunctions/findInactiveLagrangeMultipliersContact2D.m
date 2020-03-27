function homDOFsLM = findInactiveLagrangeMultipliersContact2D...
    (homDOFsLM, mesh, noDOFs, dHat_stiffMtxLM, gapFunctionInitial, ...
    segmentsContact, propContact)
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
% Returns the updated array of constrained Lagrange Multipliers DOFs based
% on the complementary contact conditions. For each inactive node it is
% computed the gap function and if it is positive then the node is 
% deactivated, i.e. its corresponding Lagrange Multiplier is removed from 
% the array of the constrained Lagrange Multipliers.
%
% Alternatively, for each active node it is found whether the corresponding 
% Lagrange Multiplier is valid, that is negative, and if not, meaning it is 
% positive, the Lagrange Multipliers is added to the array of the
% constrained Lagrange Multipliers.
%
%              Input :
%          homDOFsLM : array containing the IDs of the inactive 
%                      Lagrange Multipliers DOFs
%             strMsh : Nodes and elements in the mesh
%             noDOFs : Number of primal DOFs
%    dHat_stiffMtxLM : Vector of DOFs containing both displacement and
%                      Lagrange Multipliers
% gapFunctionInitial : Array containing the intial gap function between the
%                      potentially contact nodes and the rigid segments
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
%          homDOFsLM : Updated array containing the IDs of the inactive 
%                      Lagrange Multipliers DOFs
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all rigid segments
% ->
%    1i. Loop over all potential contact nodes
%    ->
%        1i.1. Get the id of the corresponding Lagrange Multiplier DOF
%
%        1i.2. Get the index of the Lagrange Multiplier DOF in the array of prescribed Lagrange Multipliers DOFs
%
%        1i.3. Check if the potential contact node is active or inactive using the corresponding Lagrange Multipliers DOF
%
%        1i.4. Update the Lagrange Multipliers DOF counter
%    <-
% <-
%
%
%% Function main body

%% 0. Read input

% Set contact tolerance
tolerance = sqrt(eps);

% Initialize counters
counterLM = 1;
counterActiveNodes = 1;

%% 1. Loop over all rigid segments
for iSeg = 1:segmentsContact.numSegments
    %% 1i. Loop over all potential contact nodes
    for iNode = 1:propContact.numNodes
        %% 1i.1. Get the id of the corresponding Lagrange Multiplier DOF
        idLM = noDOFs + counterLM;
        
        %% 1i.2. Get the index of the Lagrange Multiplier DOF in the array of prescribed Lagrange Multipliers DOFs
        indexDOFLM = find(idLM == homDOFsLM, 1);
        
        %% 1i.3. Check if the potential contact node is active or inactive using the corresponding Lagrange Multipliers DOF
        if dHat_stiffMtxLM(idLM) == 0
            % Get the DOF numbering of the current node
            DOF = 2*propContact.nodeIDs(iNode) - 1 : 2*propContact.nodeIDs(iNode);

            % Get the displacement of the current node
            u = dHat_stiffMtxLM(DOF);

            % Get the coordinates of the current node
            node = mesh.nodes(propContact.nodeIDs(iNode),1:2);

            % Compute the displaced coordinates of the current node
            nodeDisp = node + u';

            % Get the end vertices of the current rigid segment
            vertexA = segmentsContact.points(1,:,iSeg);
            vertexB = segmentsContact.points(2,:,iSeg);

            % Compute the gap function
            penetration = gapFunctionInitial(counterLM,1) - segmentsContact.normals(iSeg,:)*u;

            % Find out if penetration did occur
            [~,isIntersection] = computeIntersectionBetweenStraightLines...
                (node,nodeDisp,vertexA,vertexB);

            % Decide whether the node is active upon the different conditions
            if isIntersection && penetration > tolerance
                homDOFsLM(indexDOFLM) = [];
            else
                homDOFsLM(counterActiveNodes) = noDOFs + counterLM;
                counterActiveNodes = counterActiveNodes + 1;
            end
        else
            % Deactivate the Lagrange Multipliers DOF if its value is
            % positive
            if dHat_stiffMtxLM(idLM) > 0
                homDOFsLM(counterActiveNodes) = noDOFs + counterLM;
                counterActiveNodes = counterActiveNodes + 1;
            end
        end
        
        %% 1i.4. Update the Lagrange Multipliers DOF counter
        counterLM = counterLM + 1;
    end
end

end