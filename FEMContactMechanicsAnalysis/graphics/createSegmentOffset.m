function segmentOffset = createSegmentOffset(segment, normal)
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
% Returns the coordinates of the vertices of a straight segment offsetted
% in the opposite of the normal direction of the provided straight segment,
% which helps in determining in which side of the segment the the rigid 
% body lies.
%
%         Input : 
%       segment : Row vector containing the coordinates of the end vertices 
%                 of a straight line segment
%        normal : Unit normal vector to the straight segment
%
%        Output :
% segmentOffset : Matrix with the coordinates of two wall determining 
%                 points offset by a normal * scaling factor
%
% Function layout :
%
% 1. Assign a scaling factor determining the offset distance
%
% 2. Compute the vertices of the offsetted segment
%
% 3. Return the segment vertices in an array
%
%% Function main body

%% 1. Assign a scaling factor determining the offset distance
scalingFactor = 0.05;

%% 2. Compute the vertices of the offsetted segment
x0 = segment(1,1:2)- normal*scalingFactor;
x1 = segment(1,3:4)- normal*scalingFactor; 

%% 3. Return the segment vertices in an array
segmentOffset = [x0
                 x1];

end