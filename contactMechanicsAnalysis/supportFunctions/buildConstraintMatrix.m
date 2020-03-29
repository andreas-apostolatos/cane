function C = buildConstraintMatrix(mesh,propContact,segmentsContact)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%                  Marko Leskovar
%
% Date: 07.03.2020
%
%% Function documentation
%
% Returns the coupling matrix for a set of active nodes corresponding to a
% Lagrange Multipliers formulation in a 2D contact problem. If the set of
% active nodes is empty it builds the full coupling matrix.
%
%             Input :
%            noDOFs : Number of Degrees of Freedom for the primal problem
%       propContact : Data structure containing the contact properties,
%                           .nodeIds : global numbering of contact nodes
%                     .numberOfNodes : number of nodes
%   segmentsContact : Data sturcture containing information about the 
%                     boundaries of the rigid wall :
%                           .points : a list of 2x2 matrices containing
%                                      end points of the segment(s)
%                            .number : total number of segments
%                           .normals : normal vector of each segment
%
%            Output :
%                 C : The assembled coupling matrix corresponding to the
%                     Lagrange Multipliers method
%
%% Function layout :
%
% 0. Read input
%
% 1. Loop over all rigid segments
% ->
%    1i. Loop over all potential contact nodes
%    ->
%        1i.1. Check if the current DOFs is a member of active nodes
%    
%        1i.2. Find the DOFs of the node
%
%        1i.3. Assemble the entry to the coupling matrix
%
%        1i.4. Update counter which counts the Lagrange Multipliers DOFs
%
%        1i.5. Update counter which counts total number of nodes
%    <-
% <-
%
%% Function main body

%% 0. Read input

% Number of primal DOFs
noDOFs = 2*length(mesh.nodes(:,1));

% Number of Lagrange Multipliers DOFs
noDOFsLM = numel(propContact.isActive);

% Initialize coupling matrix
C = zeros(noDOFs, noDOFsLM);

% Initialize counter
counterLM = 1;
counterNodes = 1;

%% 1. Loop over all rigid segments
for iSeg = 1:segmentsContact.numSegments
    %% 1i. Loop over all potential contact nodes
    for iNode = 1:propContact.numNodes
        %% 1i.2. Find the DOF ids at the node
        DOFs = 2*propContact.nodeIDs(iNode) - 1 : 2*propContact.nodeIDs(iNode);

        %% 1i.3. Assemble the entry to the coupling matrix
        C(DOFs,counterLM) = - segmentsContact.normals(iSeg,:);

        %% 1i.4. Update counter which counts the Lagrange Multipliers DOFs
        counterLM = counterLM + 1;

        %% 1i.5. Update counter which counts total number of nodes
        counterNodes = counterNodes + 1;
    end 
end

end 