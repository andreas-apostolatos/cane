function contactSegments = createCircleSegments ...
    (center, radius, startAngle, endAngle, numSegments)
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
%   Creates a circular set of segments. If segments are orderd from left to
%   right, the normal vector points upwards.
% 
%              Input :
%             center : vector in form of [x0,y0]
%
%             Output :
%    contactSegments : data stucture containing informations about segments
%                      points and normal vector
%
%% Function main body

x_translation = 0;
y_translation = 0;

contactSegments.numSegments = numSegments;

% Initialize coordinates
coordinates = zeros(numSegments+1,2);

% Create points on the circle
angles = linspace(startAngle,endAngle,numSegments+1); 
coordinates(:,1) = radius * cos(angles) + center(1); 
coordinates(:,2) = radius * sin(angles) + center(2);

% Create segments from points
for iSeg = 1:numSegments
    
    % get rotated coordinates
    x0 = coordinates(iSeg,1) + x_translation;
    y0 = coordinates(iSeg,2) + y_translation;
    x1 = coordinates(iSeg+1,1) + x_translation;
    y1 = coordinates(iSeg+1,2) + y_translation;
    
    % add each segment to .points
    contactSegments.points(iSeg,:) = [x0 y0 x1 y1];
end

end
