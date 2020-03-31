function [dHistory, minElSize] = ...
    solve_FEMPlateInMembraneActionNLinearTransient...
    (propAnalysis, strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, ...
    computeLoadVct, parameters, computeBodyForces, propNLinearAnalysis, ...
    propStrDynamics, solve_LinearSystem, propGaussInt, propOutput, ...
    caseName, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the displacement field and the minimum element area size
% corresponding to the solution of the plate in membrane action problem
% considering a geometrically nonlinear transient analysis.
%
%              Input :
%       propAnalysis : Structure containing general information about the 
%                      analysis,
%                         .type : The analysis type
%             strMsh : Nodes and elements in the mesh
%            homDOFs : The global numbering of the nodes where homogeneous
%                      Dirichlet boundary conditions are applied 
%          inhomDOFs : The global numbering of the nodes where 
%                      inhomogeneous Dirichlet boundary conditions are 
%                      applied
%    valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                      Dirichlet boundary conditions are applied
%            propNBC : Structure containing properties on the Neumann 
%                      boundary conditions,
%                       .nodes : The nodes where Neumann boundary 
%                                   conditions are applied
%                    .loadType : The type of the load for each Neumann node
%                   .fctHandle : The function handle for each Neumann node
%                                for the computation of the load vector 
%                                (these functions are under the folder 
%                                load)
%      computeLoadVct : Function handle to the load vector computation
%          parameters : Problem specific technical parameters
%   computeBodyForces : Function handle to the computation of the body 
%                       force vector
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
%        propGaussInt : Number of Gauss Points for the integral evaluation
%                        .domainNoGP : Number of Gauss Points for the
%                                      evaluation of the domain integrals
%                      .boundaryNoGP : Number of Gauss Points for the
%                                      evaluation of the boundary integrals
%          propOutput : Structure containing information on writting the
%                       results for postprocessing,
%                                 .isOutput : Flag on whether the results 
%                                             to be written out
%                        .writeOutputToFile : Function handle to the
%                                             writting out of the results
%                                             in a VTK format
%            caseName : The name of the case in the inputGiD case folder
%              outMsg : On outputting information
%
%              Output :
%            dHistory : The history of the displacement field
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
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________\n');
    fprintf('###################################################################\n');
    fprintf('Computation of the displacement field for a geometrically nonlinear\n');
    fprintf('plate in transient membrane action problem has been initiated\n');
    fprintf('___________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Number of nodes in the mesh
noNodes = length(strMsh.nodes(:, 1));

% Number of DOFs in the mesh
nDOFs = 2*noNodes;

% GLobal DOF numbering
DOFNumbering = 1:nDOFs;

% Path where to write the output
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

% Define tabulation for outputting on the command window
tab = '\t'; 

% Assign dummy variables
computeConstantMatrices = 'undefined';
propALE = 'undefined';
computeUpdatedMesh = 'undefined';

% Title for the output file
title = 'Geometrically nonlinear transient plane stress analysis';

% Get the DOF numbering for each component of the displacement field and
% the pressure seperately
DOF4Output = [1:2:nDOFs - 1
              2:2:nDOFs];

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs, inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs, prescribedDoFs)) = [];

%% 2. Solve the transient problem
[dHistory,minElSize] = solve_FEMTransientAnalysis...
    (propAnalysis, strMsh, DOFNumbering, freeDOFs, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, propALE, @computeInitCndsFEMPlateInMembraneAction, ...
    computeBodyForces, propNBC, computeLoadVct, parameters, ...
    @solve_FEMNLinearSystem, computeConstantMatrices, ...
    @computeMassMtxFEMPlateInMembraneAction, ...
    @computeTangentStiffMtxResVctFEMPlateInMembraneAction, ...
    computeUpdatedMesh, solve_LinearSystem, propStrDynamics, ...
    propNLinearAnalysis, propGaussInt, propOutput, caseName, pathToOutput, ...
    title, DOF4Output, tab, outMsg);

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('\nNonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('________________________Linear Analysis Ended______________________\n');
    fprintf('####################################################################\n\n\n');
end

end
