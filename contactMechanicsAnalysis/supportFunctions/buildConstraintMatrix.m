function C = buildConstraintMatrix(noDOFs,propContact,nodesActive,segmentsContact)
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
% Returns the coupling matrix corresponding to a Lagrange Multipliers
% formulation for a 2D contact problem. 
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
% Function layout :
%
% 0. Read input
%
% 1. Loop over all rigid segments
% ->
%    1i. Loop over all potential contact nodes
%    ->
%        1i.1. Find the DOFs of the node
%
%        1i.2. Assemble the entry to the coupling matrix
%
%        1i.3. Update counter which counts the Lagrange Multipliers DOFs
%    <-
% <-
%
%% Function main body

%% 0. Read input

% Number of Lagrange Multipliers DOFs
noDOFsLM = propContact.numberOfNodes;

% Number of active nodes
noActiveDOFs = length(nodesActive);

% Initialize activeNodes conditions
isCnd_empty = isempty(nodesActive);

% Initialize coupling matrix
if(isCnd_empty)
    C = zeros(noDOFs, noActiveDOFs);
else
    C = zeros(noDOFs, noDOFsLM);
end

% Initialize counter
counterLM = 1;
l = 1;

%% 1. Loop over all rigid segments
for iSeg = 1:segmentsContact.number
    %% 1i. Loop over all potential contact nodes
    for iCN = 1:propContact.numberOfNodes
        if isCnd_empty || max(ismember(nodesActive,l))
        %% 1i.1. Find the DOFs of the node
        DOF = 2*propContact.nodeIDs(iCN) - 1 : 2*propContact.nodeIDs(iCN);
        
        %% 1i.2. Assemble the entry to the coupling matrix
        C(DOF,counterLM) = segmentsContact.normals(iSeg,:);
        
        %% 1i.3. Update counter which counts the Lagrange Multipliers DOFs
        counterLM = counterLM + 1;
        end
        l=l+1;
    end 
end

end 