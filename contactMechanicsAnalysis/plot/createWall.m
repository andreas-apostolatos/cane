function points = createWall(wall)
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
% Returns rotated coordinates used in the plot_segments function to
% visualize the orientation of a line segment [(x0,y0) ; (x1,y1)]
%
%   Input : 
%    wall : Matrix with the coordinates of two wall determining points
%
%  Output :
%  points : Matrix with the coordinates of two wall determining points
%           ofset by a normal * scaling factor
%
%% Function main body

% Initialize the output arrays for the rigid wall marking lines
x0(1,2)=0;
x1(1,2)=0;
% scaling factor to determine distance between lines
scalingFactor = 0.05;

% change of distance in x and y direction
DX = wall(2,1)-wall(1,1);
DY = wall(2,2)-wall(1,2);
% normalized segment normal
normal = [-DY,DX];
normal = normal/norm(normal);
% Coordinates for the rigid wall marker lines
x0 = wall(1,:)- normal*scalingFactor;
x1 = wall(2,:)- normal*scalingFactor; 

% output result
points = [x0;x1];

end