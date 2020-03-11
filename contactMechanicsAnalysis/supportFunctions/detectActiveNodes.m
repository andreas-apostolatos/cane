function activeNodes = detectActiveNodes(mesh,propContact,displacement,segments)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
% Date : 07.03.2020
%
%% Function documentation
%
% Detect the active nodes for the current segment. Loops over all segments
% to check whether the displacement of the node put in an active position
% according to the segment
%
%          Input :
%   contactNodes : structure containing the global numbering of contact
%                  canditate-nodes coordinates of the candidate nodes
%   displacement : Vector of the FULL displacement field
%       segments : data stucture containing informations about the
%                  rigid wall segments (normal vector, parallel vector,
%                  position)
%      
%         Output :
%    active_node : Vector containing global indices of active contact nodes
%
%% Function main body

% initialize variables
activeNodes=[];
k=1;

% set tolerance for contact
tolerance = sqrt(eps);

% loop over displacement_exp vector
% loop over the number of contact nodes segments
for m=1:segments.number
    % loop over every node in that segment
    for n=1:propContact.numberOfNodes
        
        % get displacement of the node of interest
        DOF = 2*propContact.nodeIDs(n)-1 : 2*propContact.nodeIDs(n);
        u = displacement(DOF);
        
        % get node of interest - P
        P = mesh.nodes(propContact.nodeIDs(n),1:2);
        
        % add displacement to the node of interest
        R = P + u';
        
        % get the start and the end point of each segment - A and B
        A = segments.points(1,:,m);
        B = segments.points(2,:,m);
        
        % projection of point on the segment - Rs
        lambda = ((A-B)*(R-A)')/((A-B)*(B-A)');
        Rs = (1-lambda)*A+lambda*B;
        
        % compute distance to the segment - normal*vector
        gap = segments.normals(m,:)*(R-Rs)';
        
        % conditions for geometry (penetration)
        isCnd1 = gap < tolerance;
        isCnd2 = lambda > tolerance;
        isCnd3 = lambda <= 1;
        
        % if all conditions hold
        if (isCnd1 && isCnd2 && isCnd3)
            activeNodes = [activeNodes,k];
        end
        
        % update counter
        k=k+1;
    end
end

end