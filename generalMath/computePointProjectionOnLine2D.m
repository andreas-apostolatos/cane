function [point_proj, lambda] = computePointProjectionOnLine2D...
    (point, vertexA, vertexB)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the projected coordinates of a point in 2D space and the
% corresponding line parameter corresponding to its projection on a line
% which is defined by its two vertices. If the line parameter lies in the
% interval [0, 1] the node has a projection on the line segment
% vertexA - vertexB otherwise it's projection lies outside these vertices.
%
% Date: 07.03.2020
%
%            Input :
%            point : Array of 2 coordinates defining the 2D point to be 
%                    projected on line
% vertexA, vertexB : Vertices defining the straight line
%
%           Output :
%       point_proj : Coordinates of the projected point onto the line in 2D
%           lambda : Line parameter
%
% Function layout :
%
% 1. Compute the line parameter in the interval [0, 1]
%
% 2. Compute the coordinates of the projected node on the line defined by vertices vertexA and vertexB
%
%% Function main body

%% 1. Compute the line parameter in the interval [0, 1]
lambda = ((vertexA - vertexB)*(vertexA - point)')/norm(vertexA - vertexB)^2;

%% 2. Compute the coordinates of the projected node on the line defined by vertices vertexA and vertexB
point_proj = (1 - lambda)*vertexA + lambda*vertexB;

end

