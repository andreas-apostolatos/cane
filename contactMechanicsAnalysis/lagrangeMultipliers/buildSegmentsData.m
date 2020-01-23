function segments = buildSegmentsData(segments)
%% Function documentation
%
%BUILDSEGMENTSDATA Build a data structure segments
%   The segments data structure contains all information required for next
%   steps aka finding the violating nodes
% 
%              Input :
%             points : points(i,:,:) is a list of 2x2 matrices containing
%                      extremities points of the segment(s)
%
%             Output :
%           segments : data stucture containing informations about the
%                      rigid wall segments (normal vector, parallel vector,
%                      position)
%
%% Function main body
% create a temporary array to initialize segmetns struct
tmp=zeros(size(segments.points,3),2);
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

 