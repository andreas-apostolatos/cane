function C = computeLagrangeMultipliersMatrixContact2D ... 
    (numDOFs, numDOFsLM, propContact, segmentsContact)
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
% Returns the complete Lagrange Multipliers matrix for enforcing the
% contact constraints corresponding to a 2D contact problem.
%
%             Input :
%           numDOFs : Number of DOFs for the primal problem
%         numDOFsLM : Numner of DOFs for the Lagrange Multipliers field
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
% Function layout :
%
% 0. Read input
%
% 1. Loop over all rigid segments
% ->
%    1i. Loop over all potential contact nodes
%    ->
%        1i.1. Get the id of the corresponding Lagrange Multipliers DOF
%    
%        1i.2. Find the DOFs of the node
%
%        1i.3. Assemble the entry to the coupling matrix
%    <-
% <-
%
%% Function main body

%% 0. Read input

% Initialize coupling matrix
C = zeros(numDOFs, numDOFsLM);

%% 1. Loop over all potential contact nodes
for iSeg = 1:segmentsContact.numSegments
    %% 1i. Loop over all rigid segments
    for iNode = 1:propContact.numNodes
        %% 1i.1 Get the id of the corresponding Lagrange Multipliers DOF
        idLM = (iSeg - 1)*propContact.numNodes + iNode;
        
        %% 1i.2. Find the DOF ids at the node
        DOFs = 2*propContact.nodeIDs(iNode) - 1 : 2*propContact.nodeIDs(iNode);

        %% 1i.3. Assemble the entry to the coupling matrix
        C(DOFs,idLM) = - segmentsContact.normals(iSeg,:);
    end 
end

end 