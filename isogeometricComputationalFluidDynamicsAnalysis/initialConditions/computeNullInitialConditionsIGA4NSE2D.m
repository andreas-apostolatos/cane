function [up, upDot, upDDot, noTimeStep, propFldDynamics] = ...
    computeNullInitialConditionsIGA4NSE2D ...
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
% Returns vectors of zeros for the application of null initial conditions
% for an isogeometric Navier-Stokes flow in 2D.
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
%% Function main body

noTimeStep = 0;
up = zeros(numDOFs,1);
upDot = zeros(numDOFs,1);
upDDot = 'undefined';

end