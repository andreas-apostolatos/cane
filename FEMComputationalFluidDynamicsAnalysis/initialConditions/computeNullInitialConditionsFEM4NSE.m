function [up, upDot, upDDot, numTimeStep] = ...
    computeNullInitialConditionsFEM4NSE...
    (propAnalysis, fldMsh, DOF4Output, parameters, fldDynamics, VTKResultFile, ...
    caseName, pathToFile)
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
% for a classical finite element discretization for the Navier-Stokes flow 
% in 2D.
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
% 1. Compute null initial conditions
%
%% Function main body

%% 0. Read input

% Get starting time of the simulation
numTimeStep = 0;

% Number of DOFs per node
if strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
    noDOFsNode = 3;
elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
    noDOFsNode = 4;
end

%% 1. Compute null initial conditions

% Number of nodes
noNodes = length(fldMsh.nodes(:,1));

% Number of degees of freedom
noDOFs = noDOFsNode*noNodes;

% Initialize output arrays
up = zeros(noDOFs, 1);
upDot = zeros(noDOFs, 1);
upDDot = 'undefined';

end
