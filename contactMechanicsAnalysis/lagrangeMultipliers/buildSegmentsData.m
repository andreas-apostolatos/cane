function segments = buildSegmentsData(segmentPoints)
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
%%
tmp=zeros(size(segmentPoints,3),2);
segments=struct('normals',tmp,'directors',tmp,'constants',tmp(:,1),'number',0);
for i=1:size(segmentPoints,3)
    % change of distance in x and y direction
    DX = segmentPoints(2,1,i)- segmentPoints(1,1,i);
    DY = segmentPoints(2,2,i)- segmentPoints(1,2,i);
    % normalized segment director
    director = [DX DY];
    director = director/norm(director);
    % normalized segment normal
    normal = [-DY DX];
    normal = normal/norm(normal);
    % assign variables
    segments.directors(i,:)=director;
    segments.normals(i,:)=normal;
    segments.constants(i)=-dot(normal,segmentPoints(1,:,i));
end
segments.number=length(segments.constants); % delete this, not used
end

 