function propContact = removeFullyConstrainedNodes(homDBC,propContact)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
% Date: 07.03.2020
%
%% Function documentation
%
%   Removes the fully constrained nodes from the set of possible contact
%   nodes (propContact). Nodes that have both DOFs fixed (x,y) are rigid
%   and do not need to be tested for possible contact
% 
%              Input :
%             homDBC : Vector of the prescribed DoFs (by global numbering)
%        propContact : structure containing the global numbering of the
%                      canditate contact nodes
%
%             Output :
%        propContact : structure containing the global numbering of the
%                      canditate contact nodes
%
%% Function main body

% Initialize a logical array of zeros(false)
fullyConstrainedNodes = false(propContact.numberOfNodes,1);

% Loop through all the potential contact nodes
for n=1:propContact.numberOfNodes
    
    % Find the coresponding degrees of freedom  
    DOFs = 2*propContact.nodeIDs(n)-1 : 2*propContact.nodeIDs(n);
    
    % Check if both DOFs are part of homogeneous boundary conditions 
    if(min(ismember(DOFs,homDBC)) == true) 
        fullyConstrainedNodes(n) = true; 
    end
end

% Delete the coresponding nodeIDs and calcualte new number of nodes
propContact.nodeIDs(fullyConstrainedNodes) = [];
propContact.numberOfNodes = length(propContact.nodeIDs);

end