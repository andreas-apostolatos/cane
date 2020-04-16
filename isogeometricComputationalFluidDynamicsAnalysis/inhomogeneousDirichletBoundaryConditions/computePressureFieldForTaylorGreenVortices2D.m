function p = computePressureFieldForTaylorGreenVortices2D ...
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
% Returns the pressure field for the Taylor-Green vortex benchmark
%
%        Input :
%   parameters : Flow parameters,
%                   .nue : Dynamic viscosity
%        x,y,z : Cartesian coordinates of the node where the prescribed 
%                (inhomogeneous) Dirichlet boundary conditions are applied
%            t : The time instance where to compute the prescribed 
%                Dirichlet boundary conditions
%
%     Output :
%          p : The pressure field    
%
%% Function main body

% Get the dynamic viscosity
nue = parameters.nue;

% Compute the pressure field
p = -.25*(cos(2*x) + cos(2*y))*exp(-4*t*nue);

end