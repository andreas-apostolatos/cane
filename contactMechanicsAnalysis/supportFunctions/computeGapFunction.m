function  propContact = computeGapFunction(mesh,propContact,segments)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% Expands a data structure contactNodes. Gap function is the distance
% between nodes of the structure and the segments points (walls)
%
%              Input :
%       contactNodes :
%           .indices : global numbering of the contact canditate-nodes
%         .positions : coordinates of the candidate nodes
%
%           segments : data stucture containing informations about the
%                      rigid wall segments (normal vector, parallel vector,
%                      position)
%            .points : points(:,:,i) is a list of 2x2 matrices containing
%                      end points of the segment(s)
%
%             Output :
%   contactNodes.gap : Normal distance of every node to each segment
%
%% %% Function main body

% loop through the number of segments
for m=1:segments.number
    % loop through the contact nodes
    for n=1:size(propContact.nodeIDs,1)
        
        % get node of interest - P
        P = mesh.nodes(propContact.nodeIDs(n),1:2);
        
        % get the start and the end point of each segment - A and B
        A = segments.points(1,:,m);
        B = segments.points(2,:,m);
        
        % projection point on the segment - Ps
        lambda = ((A-B)*(P-A)')/((A-B)*(B-A)');
        Ps = (1-lambda)*A+lambda*B;
        
        % compute distance to the segment - normal*vector
        propContact.gap(n,m) = segments.normals(m,:)*(P-Ps)';
    end
end

end