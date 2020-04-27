function [contactLength, contactForce, maxContactPressure] = ...
    computePostprocResultantsSignoriniFrictionlessContact2D...
    (mesh, parameters, dHat, lambdaHat, nodeIDs_active)
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
% Computes the length of contact, its reaction force and the maximal
% pressure. Used in unit test for the comparison with analytical values
% according to Hertz contact theory. Cannot detect an element that has all
% three nodes active.
%
%              Input :
%               mesh : Nodes and elements in the mesh
%         parameters : Problem specific technical parameters
%               dHat : Nodal displacement field
%          lambdaHat : Nodal Lagrange Multipliers field
%     nodeIDs_active : The IDs of the active (contact) nodes
%
%             Output :
%      contactLength : Length of contact surface
%       contactForce : The total reaction force on the contact 
% maxContactPressure : The maximum compressive pressure on the contact
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all active contact nodes
% ->
%    1i. Find the elements indices to which nodeI belongs to
%
%   1ii. Loop over remaining to the right active contact nodes
%   ->
%        1ii.1. Find the elements indices to which the nodeJ belongs to
%
%        1ii.2. Loop over all elements containing nodeJ
%        ->
%               1iii.2i. Find the common elements to which both nodes belong
%
%              1iii.2ii. Check if both nodeI and nodeJ share only one element in common and if yes compute all postprocessing resultants at the contact element level
%        <-
%   <-
% <-
%
% 2. Get the maximum contact pressure from finite element edges
%
% 3. Compute the total contact force by means of the Lagrange Multipliers solution field
%
%% Function main body

%% 0. Read input

% Initialize output arrays
contactPressure = [];
contactLength = 0;

%% 1. Loop over all active contact nodes
for i = 1:length(nodeIDs_active)
    %% 1i. Find the elements indices to which nodeI belongs to
    [indexI, ~] = find(nodeIDs_active(i) == mesh.elements);
    
    %% 1ii. Loop over remaining to the right active contact nodes
    for j = i + 1:length(nodeIDs_active)
        %% 1ii.1. Find the elements indices to which the nodeJ belongs to
        [indexJ, ~] = find(nodeIDs_active(j) == mesh.elements);
        
        %% 1ii.2. Loop over all elements containing nodeJ
        for k = 1:length(indexJ)
            %% 1iii.2i. Find the common elements to which both nodes belong
            commonElmnts = find(indexJ(k) == indexI);
            
            %% 1iii.2ii. Check if both nodeI and nodeJ share only one element in common and if yes compute all postprocessing resultants at the contact element level
            if length(commonElmnts) == 1
                % Compute the nodal coordinates of the active nodes on the finite element edge
                nodeI = mesh.nodes(nodeIDs_active(i), 1:2)';
                nodeJ = mesh.nodes(nodeIDs_active(j), 1:2)';
                
                % Compute the displacement field of the active nodes on the finite element edge
                dI = [dHat(2*nodeIDs_active(i) - 1, 1)
                      dHat(2*nodeIDs_active(i), 1)];
                dJ = [dHat(2*nodeIDs_active(j) - 1, 1)
                      dHat(2*nodeIDs_active(j), 1)];

                % Compute the displaced nodal coordinates of the active nodes on the finite element edge
                nodeI_disp = nodeI + dI;
                nodeJ_disp = nodeJ + dJ;
                  
                % Compute the length of the contact element
                contactLength_el = norm(nodeJ_disp - nodeI_disp);
                
                % Compute the contact pressure acting on each node using the Langrange Multipliers solution
                contactPressure_el = (lambdaHat(i, 1) + lambdaHat(j, 1))/2/parameters.t/contactLength_el; 
                contactPressure = [contactPressure; contactPressure_el];

                % Update the contact length
                contactLength = contactLength + contactLength_el;
            end
        end
    end
end

%% 2. Get the maximum contact pressure from finite element edges
maxContactPressure = max(contactPressure);

%% 3. Compute the total contact force by means of the Lagrange Multipliers solution field
contactForce = sum(lambdaHat);

end
