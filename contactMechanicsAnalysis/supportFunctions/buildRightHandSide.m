function F_exp = buildRightHandSide(F,propContact,nodesActive,segmentsContact)
%% Function documentation
%
% Add the gap function of every Lagrange multiplier to the RHS load vector
% 
%             Input :
%                 F : Global force vector
%      contactNodes : structure containing the global numbering of contact
%                     canditate-nodes coordinates of the candidate nodes
%       activeNodes : List of indices of the nodes for which the matrix
%                     should be built (the set of all currently active nodes)
%   segmentsContact : Data sturcture containing information about the 
%                     boundaries of the rigid wall :
%                           .points : a list of 2x2 matrices containing
%                                      end points of the segment(s)
%                            .number : total number of segments
%                           .normals : normal vector of each segment
%
%            Output :
%             F_exp : Vector containing in its first entries the force on
%                     every DOF  and in its last entries the gap constant
%                     for all nodes defined in the activeNodes
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
%        1i.2. Assemble the entry to the expanded force vector
%
%        1i.3. Update counter which counts the Lagrange Multipliers DOFs
%
%        1i.4. Update counter which counts total number of nodes
%    <-
% <-
%
%% Function main body

%% 0. Read input

% get number of DOFs
noDOFs = length(F);

% Initialize activeNodes conditions
isCnd_empty = isempty(nodesActive);

% Initialize F_exp
F_exp = [F
         zeros(length(nodesActive),1)];

% Initialize counters
counterLM = 1;
counterNodes = 1;

%% 1. Loop over all rigid segments
for iSeg = 1:segmentsContact.number
    %% 1i. Loop over all potential contact nodes
    for iCN = 1:propContact.numberOfNodes
        %% 1i.1. Check if the current DOFs is a member of active nodes
        if isCnd_empty || max(ismember(nodesActive,counterNodes))
            %% 1i.2. Assemble the entry to the expanded force vector
            F_exp(noDOFs+counterLM) = -propContact.gap(iCN,iSeg);
            
            %% 1i.3. Update counter which counts the Lagrange Multipliers DOFs
            counterLM = counterLM + 1;
        end
         %% 1i.4. Update counter which counts total number of nodes
        counterNodes = counterNodes + 1;
    end
end

end