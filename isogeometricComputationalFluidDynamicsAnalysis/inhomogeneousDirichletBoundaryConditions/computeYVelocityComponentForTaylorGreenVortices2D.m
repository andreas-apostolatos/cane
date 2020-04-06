function v = computeYVelocityComponentForTaylorGreenVortices2D ...
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
% Returns the y-component of the velocity field for the Taylor-Green vortex 
% benchmark
%
%      Input :
%   parameters : Flow parameters,
%                   .nue : Dynamic viscosity
%      x,y,z : Cartesian coordinates of the node where the prescribed 
%              (inhomogeneous) Dirichlet boundary conditions are applied
%          t : The time instance where to compute the prescribed Dirichlet 
%              boundary conditions
%
%     Output :
%          v : The y-component of the velocity field    
%
%% Function main body

% Get the dynamic viscosity
nue = parameters.nue;

% Compute the y-component of the velocity field
v = sin(x)*cos(y)*exp(-2*t*nue);

end