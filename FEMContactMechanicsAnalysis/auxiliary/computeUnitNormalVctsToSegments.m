function segments = computeUnitNormalVctsToSegments(segments)
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
% Returns the updated variable of the rigid segments with the unit normal
% vectors
% 
%           Input :
% contactSegments : Containts the coordinates of the vertices of the 
%                   straight segments which form the rigid body's boundary:
%                  .numSegments : Number of straight segments
%                       .points : Array containing the coordinates of the 
%                                 vertices of each line segment X0 - X1 in 
%                                 the form [x0, y0 ; x1 , y1]
%
%          Output :
% contactSegments : Updated data structure containing:
%                       .normals : Array of the coordinates of each outward 
%                                  unit normal vector
%
% Function layout :
%
% 0. Read input
%
% 1. loop over all the rigid segments
% ->
%    1i. Compute size of the segment along x- and y- Cartesian coordinates
%   1ii. Compute the not normalized outward normal to the segment vector
% <-
%
%% Function main body

%% 0. Read input

% prealocate vector of normals and directors
segments.normals = zeros(segments.numSegments,2);

%% 1. loop over all the segments
for iSeg = 1:segments.numSegments    
    %% 1i. Compute size of the segment along x- and y- Cartesian coordinates
    dx = segments.points(iSeg, 3) - segments.points(iSeg ,1);
    dy = segments.points(iSeg, 4) - segments.points(iSeg, 2);
    
    %% 1ii. Compute the not normalized outward normal to the segment vector
    segments.normals(iSeg,:) = [-dy dx]/norm([-dy dx]);
end

end