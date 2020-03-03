function segments = createCircleSegments(xCenter,yCenter,radius,nSegments)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function main body

% Initialize coordinates
coordinates = zeros(nSegments,2);

% Create points on the circle
angles = linspace(3*pi/4, pi/4, nSegments); 
coordinates(:,1) = radius * cos(angles) + xCenter; 
coordinates(:,2) = radius * sin(angles) + yCenter;

% Create segments from points
for i=1:nSegments-1
    
    x0 = coordinates(i,1);
    y0 = coordinates(i,2);
    x1 = coordinates(i+1,1);
    y1 = coordinates(i+1,2);
    
    segments.points(:,:,i) = [x0, y0; x1,y1];

end