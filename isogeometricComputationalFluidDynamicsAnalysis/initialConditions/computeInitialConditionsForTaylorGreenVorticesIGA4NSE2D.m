function [up, upDot, upDDot, noTimeStep, propFldDynamics] = ...
    computeInitialConditionsForTaylorGreenVorticesIGA4NSE2D ... 
    (BSplinePatch, numDOFs, propFldDynamics, caseName, pathToOutput, ...
    tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the vector with the initial conditions for the case of the Taylor
% -Green vortices benchmark. Note that for the isogeometric setting only
% the bilinear elements can be used sue to the non-interpolatory nature of
% the basis functions for the high-order elements.
%
%           Input :
%    BSplinePatch : The B-Spline patch over which to compute the initial
%                   conditions corresponding to the Taylor-Green vortices 
%                   flow in 2D
%         numDOFs : Number of DOFs
% propFldDynamics : Structure containing information on the time
%                   integration of the fluid dynamics problem (dummy)
%        caseName : Name of the case (dummy)
%    pathToOutput : Path to where the results are written (dummy)
%             tab : Tabulation for writting into the command window (dummy)
%          outMsg : Enabled information on the command window if it is set
%                   as 'outputEnabled' (dummy)
%
%          Output :
%              up : The vector of DoFs containing the initial conditions 
%                   for the velocity and pressure field
%           upDot : The vector of DoFs containing the initial conditions 
%                   for the acceleration and the pressure rate field
%          upDDot : Dummny variable needed only for computational 
%                   structural dynamics
%      noTimeStep : Number of time step to start from
% propFldDynamics : Updated structure containing information on the time
%                   integration of the fluid dynamics problem (dummy)
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the DOFs and assign the initial conditions
%   
%% Function main body

%% 0. Read input

% Get the B-Spline patch
% Check the NURBS geometry input
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
end

% Get the dynamic viscosity of the problem
nue = BSplinePatch.parameters.nue;

% Get the number of Control Points in the xi-direction
numCPs_xi = length(BSplinePatch.CP(:, 1, 1));

% Initialize output arrays
up = zeros(numDOFs, 1);
upDot = zeros(numDOFs, 1);
upDDot = 'undefined';

% Number of time step
noTimeStep = 0;

%% 1. Loop over all the DOFs and assign the initial conditions
for iDOFs = 1:numDOFs
    % Get the numbering of the corresponding Control Point
    p = ceil(iDOFs/3);
    j = ceil(p/numCPs_xi);
    i = p - (j - 1)*numCPs_xi;
    dir = iDOFs - ((j - 1)*numCPs_xi + i - 1)*3;
    
    % Get the Cartesian coordinates of the Control Point
    x = BSplinePatch.CP(i, j, 1);
    y = BSplinePatch.CP(i, j, 2);
    
    % Assign the initial condition according to the nature of the DoF (i.e. 
    % if it velocity in x or y-direction or pressure)
    if dir == 1 % velocity and acceleration in x-direction
        up(iDOFs) = -cos(x)*sin(y);
        upDot(iDOFs) = 2*nue*cos(x)*sin(y);
    elseif dir == 2 % velocity and acceleration in y-direction
        up(iDOFs) = sin(x)*cos(y);
        upDot(iDOFs) = -2*nue*sin(x)*cos(y);
    elseif dir == 3 % pressure and its rate
        up(iDOFs) = -.25*(cos(2*x) + cos(2*y));
        upDot(iDOFs) = nue*(cos(2*x) + cos(2*y));
    end
end

end