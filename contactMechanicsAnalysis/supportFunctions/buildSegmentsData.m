function segments = buildSegmentsData(segments)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
% Date : 03.02.2020
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
%           segments : data stucture containing informations about the
%                      rigid wall segments (normal vector, parallel vector,
%                      position)
%
%% Function main body

% create a temporary array to initialize segmetns struct
segments.number = size(segments.points,3);
tmp=zeros(segments.number,2);
segments.normals = tmp;
segments.directors = tmp;
segments.constants = tmp(:,1);

% loop over the number of segments
for i=1:size(segments.points,3)
    % change of distance in x and y direction
    DX = segments.points(2,1,i)- segments.points(1,1,i);
    DY = segments.points(2,2,i)- segments.points(1,2,i);
    % normalized segment director
    director = [DX DY];
    director = director/norm(director);
    % normalized segment normal
    normal = [-DY DX];
    normal = normal/norm(normal);
    % assign variables
    segments.directors(i,:)=director;
    segments.normals(i,:)=normal;
    % dot product between normal and fist point of the segment - what is?
    segments.constants(i)=-dot(normal,segments.points(1,:,i));
end

end