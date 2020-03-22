function inactive_nodes = detectInactiveNodes...
    (nDOFs,mesh,propContact,displacement_exp,segments)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
% Date : 04.02.2020
%
%% Function documentation
%
% Detect the inactive nodes for the current segment. Loops over all segments
% to check whether the displacement of the node put in an active position
% according to the segment, or whether the Lagrange multipliers are valid
%
%              Input :
%              nDOFs : total number of DOFs in the system
%       contactNodes : structure containing the global numbering of contact
%                      canditate-nodes coordinates of the candidate nodes
%   displacement_exp : Vector of the FULL displacement and Lagrange multipliers
%           segments : data stucture containing informations about the
%                      rigid wall segments (normal vector, parallel vector,
%                      position)
%      
%             Output :
%      inactive_DOFs : The resulting vector containing index of the 
%                      restricted vector of mesh.boundaryNodes containing 
%                      the global indices
%
%% Function main body

% initialize variables
inactive_nodes = zeros(1, segments.number*propContact.numberOfNodes);
k=1;
l=1;

% set tolerance for contact
tolerance = sqrt(eps);

% loop over displacement_exp vector
% loop over the number of contact nodes segments
for m=1:segments.number
    % loop over every contact node in that segment
    for n=1:propContact.numberOfNodes
        
        % get displacement of the node of interest
        DOF = 2*propContact.nodeIDs(n)-1 : 2*propContact.nodeIDs(n);
        u = displacement_exp(DOF);
        
        % get node of interest - P
        P = mesh.nodes(propContact.nodeIDs(n),1:2);
        
        % add displacement to the node of interest
        R = P + u';
        
        % get the start and the end point of each segment - A and B
        A = segments.points(1,:,m);
        B = segments.points(2,:,m);
        
        % projection of point on the segment - Rs
        alpha = ((A-B)*(R-A)')/((A-B)*(B-A)');
        Rs = (1-alpha)*A+alpha*B;
        
        % compute distance to the segment - normal*vector
        gap = segments.normals(m,:)*(R-Rs)';
        
        % conditions for geometry (non-penetration)
        isCnd1 = gap > tolerance;
        isCnd2 = alpha < tolerance;
        isCnd3 = alpha >= 1;
        
        % condition for Lagrange multipliers (non-compressive)
        isCnd4 = displacement_exp(nDOFs+k) > tolerance;
        
        % if any of the conditions hold then the node is inactive
        if (isCnd1 || isCnd2 || isCnd3 || isCnd4)
            inactive_nodes(l) = nDOFs+k;
            l=l+1;
        end
        
        % update counter
        k=k+1;
    end
end

% keep only non-zero entries of the inactive_nodes
inactive_nodes = inactive_nodes(1:l-1);

end