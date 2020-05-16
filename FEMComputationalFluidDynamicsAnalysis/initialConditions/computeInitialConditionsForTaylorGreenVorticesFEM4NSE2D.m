function [up, upDot, upDDot, noTimeStep] = ...
    computeInitialConditionsForTaylorGreenVorticesFEM4NSE2D ... 
    (propAnalysis, fldMsh, DOF4Output, parameters, fldDynamics, ...
    VTKResultFile, caseName, pathToFile)
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
% Returns the vector with the initial conditions (at t=0) for the case of
% the Taylor-Green vortices benchmark in 2D.
%
%             Input :
%      propAnalysis : Structure containing general information about the 
%                     analysis,
%                        .type : The analysis type
%            fldMsh : Nodes and elements of the fluid mesh
%        DOF4Output : Array containing the arrangment of the DOFs for 
%                     printing them out
%        parameters : Flow parameters (null variable for this function)
%       fldDynamics : Transient analysis parameters for the CFD simulation
%     VTKResultFile : The name of the result file in the output folder
%                     where to get the initial conditions for the transient 
%                     simulation
%          caseName : The name of the case
%        pathToFile : Path to the file where to get the initial conditions 
%                     from
%
%            Output :
%                up : The vector of DoFs containing the initial conditions 
%                     for the velocity and pressure field
%             upDot : The vector of DoFs containing the initial conditions 
%                     for the acceleration and the pressure rate field
%            upDDot : Dummny variable needed only for computational 
%                     structural dynamics
%       numTimeStep : The number of the time step when the simulation
%                     starts
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the DOFs and assign the initial conditions
%   
%% Function main body

%% 0. Read input

% Get the dynamic viscosity of the problem
nue = parameters.nue;

% Get the total number of nodes and DOFs
numNodes = size(fldMsh.nodes, 1);
numDOFs = numNodes*propAnalysis.noFields;

% Initialize output arrays
up = zeros(numDOFs, 1);
upDot = zeros(numDOFs, 1);
upDDot = 'undefined';

% Number of time step
noTimeStep = 0;

%% 1. Loop over all the Nodes and assign initial conditions to each DOF
counter = 1;
for iNodes = 1:numNodes
    % Compute node coordinates
    x = fldMsh.nodes(iNodes, 2);
    y = fldMsh.nodes(iNodes, 3);
    
    % Assign the initial conditions at t = 0s;
    % velocity and acceleration in x-direction
    up(counter) = -cos(x)*sin(y);
    upDot(counter) = 2*nue*cos(x)*sin(y);
    
    % velocity and acceleration in y-direction
    up(counter+1) = sin(x)*cos(y);
    upDot(counter+1) = -2*nue*sin(x)*cos(y);
    
    % pressure and its rate
    up(counter+2) = -0.25*(cos(2*x) + cos(2*y));
    upDot(counter+2) = nue*(cos(2*x) + cos(2*y));
    
    % Update DOF counter
    counter = counter+3;
end

end
