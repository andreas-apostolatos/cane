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
%                      boundaries of the rigid wall:
%                   .numSegments : total number of segments
%                        .points : Array containing the coordinates of the
%                                  vertices for each rigid segment
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
%    <-
% <-
%
% 2. Clean up from duplicate entries the array of the Lagrange Multipliers DOF IDs which have to be enabled
%
% 3. Clean up from duplicate entries the array of the Lagrange Multipliers DOF IDs which have to be disabled
%
% 4. Remove the DOF IDs of the Lagrange Multipliers DOFs which have to be enabled from the array of the Lagrange Multipliers DOF IDs which are disabled
%
%% Function main body

%% 0. Read input

% warning('Make sure that no two Lagrange Multipliers for a candidate contact node are enabled at the same time');

% Set contact tolerance
tolerance = sqrt(eps);

idLM_to_be_removed = [];

%% 1. Loop over all rigid segments
for iSeg = 1:segmentsContact.numSegments
    %% 1i. Loop over all potential contact nodes
    for iNode = 1:propContact.numNodes
        %% 1i.1. Get the id of the corresponding Lagrange Multiplier DOF
        idLM = noDOFs + (iSeg - 1)*propContact.numNodes + iNode;
        
        %% 1i.2. Get the index of the Lagrange Multiplier DOF in the array of prescribed Lagrange Multipliers DOFs
        indexDOFLM = find(idLM == homDOFsLM, 1);
        
        %% 1i.3. Check if the potential contact node is active or inactive using the corresponding Lagrange Multipliers DOF
        if dHat_stiffMtxLM(idLM) == 0
            % Get the DOF numbering of the current node
            DOF = 2*propContact.nodeIDs(iNode) - 1 : 2*propContact.nodeIDs(iNode);

            % Get the displacement of the current node
            u = dHat_stiffMtxLM(DOF);

            % Get the coordinates of the current node
            node = mesh.nodes(propContact.nodeIDs(iNode),2:3);

            % Compute the displaced coordinates of the current node
            nodeDisp = node + u';

            % Get the end vertices of the current rigid segment
            vertexA = segmentsContact.points(iSeg,1:2);
            vertexB = segmentsContact.points(iSeg,3:4);

            % Compute the gap function
            penetration = gapFunctionInitial(idLM - noDOFs,1) - segmentsContact.normals(iSeg,:)*u;

            % Find out if penetration did occur
            [~,isIntersection] = computeIntersectionBetweenStraightLines...
                (node,nodeDisp,vertexA,vertexB);
%             if isIntersection && gapFunctionInitial(idLM - noDOFs,1) < tolerance && penetration < tolerance
%                 error('Node (%d,%d) with displacement (%d,%d) intersects segment with vertices (%d,%d) and (%d,%d) but the penetration is found to be %d',...
%                     node(1,1), node(1,2), u(1,1), u(2,1), vertexA(1,1), vertexA(1,2), vertexB(1,1), vertexB(1,2), penetration);
%             end

            % Decide whether the node is active upon the different conditions
%             if (isIntersection && penetration > tolerance) || abs(penetration) < tolerance
            if (isIntersection && penetration > tolerance) || abs(penetration) < tolerance
                homDOFsLM(indexDOFLM) = [];
            else
                homDOFsLM = [homDOFsLM idLM];
            end
        else
            % Deactivate the Lagrange Multipliers DOF if its value is
            % negative
            if dHat_stiffMtxLM(idLM) < tolerance
                homDOFsLM = [homDOFsLM idLM];
            else
                % Get the DOF numbering of the current node
                DOF = 2*propContact.nodeIDs(iNode) - 1 : 2*propContact.nodeIDs(iNode);

                % Get the displacement of the current node
                u = dHat_stiffMtxLM(DOF);
                
                % Get the coordinates of the current node
                node = mesh.nodes(propContact.nodeIDs(iNode),2:3); 
                
                % Compute the displaced coordinates of the current node
                nodeDisp = node + u';
                
                % Get the end vertices of the current rigid segment
                vertexA = segmentsContact.points(iSeg,1:2);
                vertexB = segmentsContact.points(iSeg,3:4);
                
                % Compute the projection of the displaced node onto the
                % segment with respect to which it is constrained
                [~, lineParameter] = computePointProjectionOnLine2D...
                    (nodeDisp, vertexA, vertexB);
                
                % Check if the node comes into contact with the segment
                if lineParameter < 0 - tolerance || lineParameter > 1 + tolerance
                    % Disable the node
                    homDOFsLM = [homDOFsLM idLM];
                    
                    % Loop over all segments to find onto which one the
                    % node has a projection
                    for iSeg2 = 1:segmentsContact.numSegments
                        % Get the end vertices of the current rigid segment
                        vertexA = segmentsContact.points(iSeg2,1:2);
                        vertexB = segmentsContact.points(iSeg2,3:4);
                        
                        % Project the displaced node onto the segment
                        [~, lineParameter] = computePointProjectionOnLine2D ...
                            (nodeDisp, vertexA, vertexB);
                        
                        % Check if the projection is valid
                        if lineParameter > 0 - tolerance && lineParameter < 1 + tolerance
                            % Get the index of the corresponding Lagrange
                            % Multiplier
                            idLM = noDOFs + (iSeg2 - 1)*propContact.numNodes + iNode;
                            
                            % Save the value to remove it at the end
                            idLM_to_be_removed = [idLM_to_be_removed idLM];
                        end
                    end
                end
            end
        end
    end
end

%% 2. Clean up from duplicate entries the array of the Lagrange Multipliers DOF IDs which have to be enabled
idLM_to_be_removed = unique(idLM_to_be_removed);

%% 3. Clean up from duplicate entries the array of the Lagrange Multipliers DOF IDs which have to be disabled
homDOFsLM = unique(homDOFsLM);

%% 4. Remove the DOF IDs of the Lagrange Multipliers DOFs which have to be enabled from the array of the Lagrange Multipliers DOF IDs which are disabled
for i = 1:length(idLM_to_be_removed)
    homDOFsLM(idLM_to_be_removed(1, i) == homDOFsLM) = [];
end

end
