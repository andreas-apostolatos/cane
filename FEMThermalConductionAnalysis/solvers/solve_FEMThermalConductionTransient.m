function [THistory, WComplete, minElSize] = ...
    solve_FEMThermalConductionTransient...
    (strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateInhomDOFs, propNBC,computeLoadVct, ...
    propParameters, computeBodyForceVct, propAnalysis, computeInitCnds, ...
    computeProblemMatricesSteadyState, propNLinearAnalysis, propIDBC, ...
    propHeatDynamics, solve_LinearSystem, solve_FEMSystem, propGaussInt, ...
    propVTK, caseName, outMsg)
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
% Returns the temperature field and the minimum element area size
% corresponding to the solution of the 2D heat transfer analysis problem
% considering a linear transient analysis.
%
%              Input :
%           analysis : Information on the analysis type
%                         .type : The analysis type
%             strMsh : Nodes and elements in the mesh
%            homDOFs : The global numbering of the nodes where homogeneous
%                      Dirichlet boundary conditions are applied 
%          inhomDOFs : The global numbering of the nodes where 
%                      inhomogeneous Dirichlet boundary conditions are 
%                      applied
%    valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                      Dirichlet boundary conditions are applied
%                NBC : On the Neumann Dirichlet boundary conditions  
%                       .nodes : The nodes where Neumann boundary 
%                                   conditions are applied
%                    .loadType : The type of the load for each Neumann node
%                   .fctHandle : The function handle for each Neumann node
%                                for the computation of the load vector 
%                                (these functions are under the folder 
%                                load)
%      computeLoadVct : Function handle to the load vector computation
%          parameters : Problem specific technical parameters
% propNLinearAnalysis : Properties of the nonlinear scheme    
%                       .scheme : The employed nonlinear scheme
%                    .tolerance : The residual tolerance
%                      .maxIter : The maximum number of the nonlinear 
%                                 iterations
%     propStrDynamics : .timeDependence : Steady-state or transient analysis
%                           .scheme : The time integration scheme
%                               .T0 : The start time of the simulation
%                             .TEnd : The end time of the simulation
%                          .nTSteps : The number of the time steps
%  solve_LinearSystem : Function handle to the linear equation system
%                       solver
%            gaussInt : Number of Gauss Points for the integral evaluation
%                        .domainNoGP : Number of Gauss Points for the
%                                      evaluation of the domain integrals
%                      .boundaryNoGP : Number of Gauss Points for the
%                                      evaluation of the boundary integrals
%            caseName : The name of the case in the inputGiD case folder
%              outMsg : On outputting information
%
%              Output :
%            THistory : The history of the temperature field
%           WComplete : The history of the energy rate field
%           minElSize : The minimum element area size in the mesh
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
if strcmp(outMsg, 'outputEnabled')
    fprintf('______________________________________________________________\n');
    fprintf('##############################################################\n');
    fprintf('Computation of the displacement field for a thermal conduction\n');
    fprintf('problem has been initiated\n');
    fprintf('______________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Title for the output file
title = 'linear transient 2D heat transfer analysis';

% Number of nodes in the mesh
noNodes = length(strMsh.nodes(:, 1));

% Number of DOFs in the mesh
numDOFs = noNodes;

% GLobal DOF numbering
DOFNumbering = 1:numDOFs;

% Path where to write the output
pathToOutput = '../../outputVTK/FEMThermalConductionAnalysis/';

% Write output to VTK
if propVTK.isOutput
    propVTK.writeOutputToFile = @writeOutputFEMThermalConductionAnalysisToVTK;
else
    propVTK.writeOutputToFile = 'undefined';
end

% Define tabulation for outputting on the command window
tab = '\t'; 

% Assign dummy variables
computeConstantProblemMatrices = 'undefined';
propALE = 'undefined';
computeUpdatedMesh = 'undefined';

% Title for the output file
title = 'linear transient 2D heat transfer analysis';

% Get the DOF numbering for each component of the temperature field
DOF4Output = 1:numDOFs;

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs, inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs, prescribedDoFs)) = [];

%% 2. Solve the transient problem
[THistory, WComplete, minElSize] = solve_FEMTransientAnalysis ...
    (propAnalysis, strMsh, DOFNumbering, freeDOFs, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, updateInhomDOFs, propALE, computeInitCnds, ...
    computeBodyForceVct, propNBC, computeLoadVct, propParameters, ...
    solve_FEMSystem, computeConstantProblemMatrices, ...
    @computeMassMtxFEMThermalConductionAnalysis, ...
    computeProblemMatricesSteadyState, computeUpdatedMesh, ...
    solve_LinearSystem, propHeatDynamics, propNLinearAnalysis, ...
    propIDBC, propGaussInt, propVTK, caseName, pathToOutput, title, ...
    DOF4Output, tab, outMsg);

%% 3. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('\nThermal conduction analysis took %.2d seconds \n\n', computationalTime);
    fprintf('________________________Linear Analysis Ended______________________\n');
    fprintf('####################################################################\n\n\n');
end

end
