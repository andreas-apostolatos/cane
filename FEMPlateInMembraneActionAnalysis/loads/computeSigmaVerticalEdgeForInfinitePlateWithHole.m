function F = computeSigmaVerticalEdgeForInfinitePlateWithHole(x,y,z,t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the analytical stress resultant sigma_yy and sigma_yx in the 
% Cartesian space for the plain stress benchmark case of an infinite plate 
% with a circular hole subject to tensional loading in +-infty.
%
%   Input : 
%   x,y,z : The Cartesian coordinates on where to evaluate the stress
%           component sigma_xy
%       t : The time instance where to compute the load
%
%  Output :
%       F : The magnitude of the nessecary force at (x,y)
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the curvilinear coordinates of (x,y)
%
% 2. Compute the stress resultants at (theta1,theta2)
%
%% Function main body

%% Compute the position of the point in the physical space

%% 0. Read input

% The load amplitude
T = 10;

% The radius of the quarter of a cirle hole
holeRadius = 1;

% Initialize the output array
F = zeros(3,1);

%% 1. Compute the curvilinear coordinates of (x,y)

% theta1 = radius
radius = sqrt(x^2+y^2);

% theta2 = angle
theta = 2*pi()-atan(y/(-x));

%% 2. Compute the stress resultants at (theta1,theta2)

% Assign some constant parameters
a = 2*theta;
rr = holeRadius/radius;
t = T/2;

% Compute the stress components over the curvilinear basis
sigmaTheta1 = t*(1-rr^2) + t*(1-4*rr^2+3*rr^4)*cos(a);
sigmaTheta2 = t*(1+rr^2) - t*(1+3*rr^4)*cos(a);
sigmaTheta1Theta2 = -t*(1+2*rr^2-3*rr^4)*sin(a);

% Compute the stress component sigma_yx over the Cartesian space
F(1,1) = -(sigmaTheta1*sin(theta)*cos(theta) - sigmaTheta2*sin(theta)*cos(theta) - sigmaTheta1Theta2*((sin(theta))^2-(cos(theta))^2));

% Compute the stress component sigma_yy over the Cartesian space
F(2,1) = sigmaTheta1*(sin(theta))^2 + sigmaTheta2*(cos(theta))^2 + 2*sigmaTheta1Theta2*sin(theta)*cos(theta);

end
