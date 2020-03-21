function  propContact = computeGapFunction(mesh,propContact,segments)
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
% Expands a data structure propContact. Gap function is the distance
% between nodes of the structure and the segments points (walls)
%
%             Input :
%              mesh : Elements and nodes of the mesh
%       propContact :       .nodeIds : global numbering of contact nodes
%                     .numberOfNodes : number of nodes 
%          segments :        .points : a list of 2x2 matrices containing
%                                      end points of the segment(s)
%                            .number : total number of segments
%                           .normals : normal vector of each segment
%            Output :
%  contactNodes.gap : Normal distance of every node to each segment
%
%% %% Function main body

% loop through the number of segments
for m=1:segments.number
    % loop through the contact nodes
    for n=1:size(propContact.nodeIDs,1)
        
        % get node of interest - R = P
        R = mesh.nodes(propContact.nodeIDs(n),1:2);
        
        % get the start and the end point of each segment - A and B
        A = segments.points(1,:,m);
        B = segments.points(2,:,m);
        
        % projection point on the segment - Ps
        alpha = ((A-B)*(R-A)')/((A-B)*(B-A)');
        Rs = (1-alpha)*A+alpha*B;
        
        % compute distance to the segment - normal*vector P-Ps
        propContact.gap(n,m) = segments.normals(m,:)*(R-Rs)';
    end
end

end