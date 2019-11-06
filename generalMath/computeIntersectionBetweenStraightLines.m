function [x,haveIntersection] = computeIntersectionBetweenStraightLines(x1,x2,x3,x4)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the intersection point between the two straight lines defined by
% the given points x1-x2 and x3-x4 in the 2-Dimensional space
%
%            Input :
%            x1,x2 : The end-points of the first line
%            x3,x4 : The end-points of the second line
%
%           Output :
%                x : The intersection point
% haveIntersection : Boolean on whether the two line actually intersect
%
%% Function main body

% Assign a tolerance value
eps = 1e-14;

% Initialize the flag to true
haveIntersection = 1;

% Form matrix J
J = [x2(1,1)-x1(1,1) x3(1,1)-x4(1,1)
     x2(1,2)-x1(1,2) x3(1,2)-x4(1,2)];

if abs(det(J)) >= eps
    % Form right-hand side
    r = [x3(1,1)-x1(1,1)
         x3(1,2)-x1(1,2)];
    
    % Solve the system for the parameters
    lambda = J\r;
    
    % Compute the locations for the intersection points
    xIntersection1 = (1-lambda(1))*x1 + lambda(1)*x2;
    xIntersection2 = (1-lambda(2))*x3 + lambda(2)*x4;
    
    % Condition for existence of an intersection point
    condition = norm(xIntersection1-xIntersection2)>eps||lambda(1)>1||lambda(1)<0||lambda(2)>1||lambda(2)<0;
    
    if condition
       x = 'The lines do not intersect';
       haveIntersection = 0;
    else
        x = xIntersection1;
    end
else
    x = 'The lines do not intersect';
    haveIntersection = 0;
end

