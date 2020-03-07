function segments = buildSegmentsData(segments)
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
%   Expands a data structure segments. The segments data structure contains
%   all information required for next steps of finding the active/inactive 
%   nodes
% 
%              Input :
%    segments.points : Matrix with the coordinates of two wall determining
%                      points in form of [(x0,y0) ; (x1,y1)]. Length in
%                      direction 3 depends on the number of walls
%
%             Output :
%           segments : data stucture containing informations about segments
%                      points and normal vector
%
%% Function main body

% assign the number of segments
segments.number = size(segments.points,3);
% prealocate vector of normals and directors
segments.normals = zeros(segments.number,2);

% loop over the number of segments
for m=1:segments.number
    
    % change of distance in x and y direction
    DX = segments.points(2,1,m)- segments.points(1,1,m);
    DY = segments.points(2,2,m)- segments.points(1,2,m);
    
    % segment normal vector
    normal = [-DY, DX];

    % assign  segment normalized normal vector
    segments.normals(m,:) = normal/norm(normal);
end

end