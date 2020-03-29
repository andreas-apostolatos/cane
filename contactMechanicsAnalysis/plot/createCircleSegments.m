function segments = createCircleSegments(center,radius,startAngle,endAngle,nSegments)
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
%   Creates a circular set of segments. If segments are orderd from left to
%   right, the normal vector points upwards.
% 
%              Input :
%             center : vector in form of [x0,y0]
%
%             Output :
%           segments : data stucture containing informations about segments
%                      points and normal vector
%
%% Function main body

x_translation = 0;
y_translation = 0;

% Initialize coordinates
coordinates = zeros(nSegments+1,2);

% Create points on the circle
angles = linspace(startAngle,endAngle,nSegments+1); 
coordinates(:,1) = radius * cos(angles) + center(1); 
coordinates(:,2) = radius * sin(angles) + center(2);

% Create segments from points
for i=1:nSegments
    
    % get rotated coordinates
    x0 = coordinates(i,1) + x_translation;
    y0 = coordinates(i,2) + y_translation;
    x1 = coordinates(i+1,1) + x_translation;
    y1 = coordinates(i+1,2) + y_translation;
    
    % add each segment to .points
    segments.points(:,:,i) = [x0, y0; x1,y1];

end