function [u, uDot, uDDot, numTimeStep, propTransientAnalysis] = ...
    computeInitCndsFEMPlateInMembraneAction...
    (analysis, strMsh, DOF4Output, propParameters, propTransientAnalysis, ...
    VTKResultFile, caseName, pathToOutput)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
%
%                 Input :
%              analysis : Structure containing general information on the 
%                         analysis,
%                               .type : Analysis type
%                strMsh : Nodes and elements in the mesh
%            DOF4Output : Arrangement of the DOFs for the output (dummy 
%                         variable for this function)
%        propParameters : Problem specific technical parameters
% propTransientAnalysis : Structure containing information on the time 
%                         integration,
%                           .timeDependence : Steady-state or transient 
%                                             analysis
%                                   .method : The time integration scheme
%                                       .T0 : The start time of the 
%                                             simulation
%                                     .TEnd : The end time of the 
%                                             simulation
%                                 .noTSteps : The number of the time steps
%         VTKResultFile : VTK result file to be read for restart (dummy for
%                         this function)
%              caseName : The name of the case in the inputGiD case folder
%          pathToOutput : Define the path to where to write out the results
%                outMsg : On outputting information
%
%                Output :
%                     u : The nodal displacement field
%                  uDot : The velocity field
%                 uDDot : The acceleration field
%           numTimeStep : The ID of the starting time step
% propTransientAnalysis : Updated structure containing information on the
%                         time integration
%
% Function layout :
%
% 0. Read input
%
% 1. Find the prescribed and the free DOFs of the system
%
% 2. Solve the transient problem
%
% 3. Appendix
%
%% Function main body
numNodes = length(strMsh.nodes(:, 1));
numDOFs = 2*numNodes;
u = zeros(numDOFs, 1);
uDot = zeros(numDOFs, 1);
uDDot = zeros(numDOFs, 1);

% Modify structure of transient analysis properties
% propTransientAnalysis.dt = propTransientAnalysis.noTimeSteps;
%propTransientAnalysis.dt = 0.1;
%propTransientAnalysis.dt =nTSteps
Tint = propTransientAnalysis.TEnd - propTransientAnalysis.T0;
propTransientAnalysis.dt = propTransientAnalysis.noTimeSteps;
propTransientAnalysis.noTimeSteps = Tint / propTransientAnalysis.dt;

numTimeStep = propTransientAnalysis.T0/propTransientAnalysis.dt;
% if(noTimeStep == 0)
%     noTimeStep = 1;

end
