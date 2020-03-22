function segments = buildSegmentsData(segments)
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
% Returns the updated variable of the rigid segments with the unit normal
% vectors
% 
%           Input :
% contactSegments : Containts the coordinates of the vertices of the 
%                   straight segments which form the rigid body's boundary:
%                       .points : Array containing the coordinates of the 
%                                 vertices of each line segment X0 - X1 in 
%                                 the form [x0, y0 ; x1 , y1]
%
%          Output :
% contactSegments : Updated data structure containing:
%                   .numSegments : Number of straight segments
%                       .normals : Array of the coordinates of each outward 
%                                  unit normal vector
%
% Function layout :
%
% 0. Read input
%
% 1. loop over all the segments
% ->
%    1i. Compute size of the segment along x- and y- Cartesian coordinates
%   1ii. Compute the not normalized outward normal to the segment vector
% <-
%
%% Function main body

%% 0. Read input

% assign the number of segments
segments.number = size(segments.points,3);

% prealocate vector of normals and directors
segments.normals = zeros(segments.number,2);

%% 1. loop over all the segments
for iSeg = 1:segments.number    
    %% 1i. Compute size of the segment along x- and y- Cartesian coordinates
    dx = segments.points(2,1,iSeg) - segments.points(1,1,iSeg);
    dy = segments.points(2,2,iSeg) - segments.points(1,2,iSeg);
    
    %% 1ii. Compute the not normalized outward normal to the segment vector
    n_tilde = [-dy dx];
    segments.normals(iSeg,:) = n_tilde/norm(n_tilde);
end

end