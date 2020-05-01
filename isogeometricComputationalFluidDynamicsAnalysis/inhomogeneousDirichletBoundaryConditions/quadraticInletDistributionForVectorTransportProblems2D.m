function scalar = ...
    quadraticInletDistributionForVectorTransportProblems2D ...
    (parameters, x, y, z, t)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the concetration of the vectorial quantity at the given location 
% y on a vertical (y-directed) line which is defined by a lower and an 
% upper value w.r.t. its y-coordinate. The scalar concetration distribution 
% along this line is assumed to be quadratic reaching its maximum value at 
% the middle point of this line.
%
%      Input :
% parameters : Flow parameters (dummy variable for this function)
%      x,y,z : Cartesian coordinates of the node where the prescribed 
%              (inhomogeneous) Dirichlet boundary conditions are applied
%          t : The time instance where to compute the prescribed Dirichlet 
%              boundary conditions
%
%     Output :
%     scalar : The concetration of the scalar quantity at the given
%              location    
%
%% Function main body

% Maximum concetration of the scalar in the center of the application line
Smax = 6e1;

% Height of the channel
Height = 10;

% Y-component of the lower edge
ylo = -Height/2;
% ylo = Height+1;

% Y-component of the upper edge
yup = Height/2;
% yup = Height+2;

% Start time of the simulation
Tstart = 10;

% end time of the simulation
Tend = 0;

% Compute the scalar value at the given location
if y>=ylo && y<=yup
    scalar = -Smax*(yup-y)*(ylo-y)/(yup-ylo)^2;%*cos(5*pi*t/(Tend-Tstart));
else
    scalar = -Smax*(yup-y)*(ylo-y)/(yup-ylo)^2;%*cos(5*pi*t/(Tend-Tstart));
end

end

