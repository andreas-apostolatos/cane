function [up, FComplete, isConverged, minElSize] = ...
    solve_FEMVMSStabSteadyStateNSE ...
    (fldMsh, up, homDOFs, inhomDOFs, valuesInhomDOFs, uMeshALE, ...
    propParameters, computeBodyForces, propAnalysis, solve_LinearSystem, ...
    propFldDynamics, propNLinearAnalysis, numIterStep, propGaussInt, ...
    propOutput, caseName, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Solve the transient nonlinear Navier-Stokes problem in 2D using the 
% Bossak time integration method for the temporal discretization and the 
% classical triangular basis functions for the spatial discretization.
%
%               Input :
%              fldMsh : Nodes and elements for the fluid mesh
%                  up : Initial conditions
%             homDOFs : The global numbering of the DOFs where homogeneous
%                       Dirichlet boundary conditions are applied
%           inhomDOFs : The global numbering of the DOFs where
%                       inhomogeneous Dirichlet boundary conditions are
%                       applied
%     valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                       Dirichlet boundary conditions are applied
%            uMeshALE : mesh velocity related to the Arbitrary
%                       Lagrangian-Eulerian method
%      propParameters : Flow parameters
%   computeBodyForces : Function handle to the computation of the body
%                       force vector
%        propAnalysis : Structure containing information about the analysis
%                            .type : The analysis type
%  solve_LinearSystem : Function handle to the solver for the linear 
%                       equation system
%     propFldDynamics : On the transient analysis :
%                             .method : The time integration method
%                          .alphaBeta : (parameter for the Bossak method)
%                              .gamma : (parameter for the Bossak method)
%                             .TStart : Start time of the simulation
%                               .TEnd : End time of the simulation
%                                 .nT : Number of time steps
%                                 .dt : Time step
% propNLinearAnalysis : On the nonlinear analysis :
%                               .method : The nonlinear solution method
%                                  .eps : The residual tolerance
%                              .maxIter : The maximum number of nonlinear
%                                         iterations
%         numIterStep : Time step number which is useful for writting out 
%                       the results if the analysis is called in a loop,
%                       like in an optimization loop
%        propGaussInt : Structure containing information on the Gaussian
%                       quadrature,
%                             .type : 'default', 'user'
%                       .domainNoGP : Number of Gauss Points for the domain
%                                     integration
%                     .boundaryNoGP : Number of Gauss Points for the
%                                     boundary integration
%          propOutput : Structure containing information on writting the
%                       results for postprocessing,
%                                 .isOutput : Flag on whether the results 
%                                             to be written out
%                        .writeOutputToFile : Function handle to the
%                                             writting out of the results
%                            .VTKResultFile : Specifies the name of the
%                                             VTK result file from which
%                                             the simulation to be
%                                             restarted. If it is
%                                             specified as 'undefined' the 
%                                             simulation starts from time 
%                                             TStart
%            caseName : String defining the case name
%              outMsg : On printing information during analysis in the
%                       command window
%
%              Output :
%                  up : The solution field in terms of the velocity and the 
%                       pressure field
%           FComplete : The complete force vector
%         isConverged : Flag on the convegence of the nonlinear system
%           minElSize : The minimum element area size over the isogeometric
%                       mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Find the prescribed and the free DOFs of the system
%
% 2. Solve the mesh motion problem and update the mesh node locations and velocities
%
% 3. Solve the steady-state nonlinear Navier-Stokes stabilized finite element equation system
%
% 4. Write out the results into a VTK file
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________\n');
    fprintf('###################################################################\n');
    fprintf('Nonlinear steady-state analysis for the FEM-based 2D incompressible\n');
    fprintf('Navier-Stokes flow problem has been initiated. \n\n');
    if strcmp(propNLinearAnalysis.method,'Newton')
        fprintf('Nonlinear method : Newton method \n');
    end
    fprintf('Residual tolerance: %d \n',propNLinearAnalysis.eps);
    fprintf('Maximum number of nonlinear iterations = %d \n',propNLinearAnalysis.maxIter);
    if propFldDynamics.isAdaptive
        fprintf('Adaptive time step chosen \n');
    else
        fprintf('Constant time step: %d (seconds) \n',propFldDynamics.dt);
    end
    fprintf('___________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMComputationalFluidDynamicsAnalysis/';

% Dummy variables
uSaved = 'undefined';
uDotSaved = 'undefined';
uDDotSaved = 'undefined';
uDot = 'undefined';
uDDot = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
precompStiffMtx = 'undefined';
precomResVct = 'undefined';
nodesSaved = 'undefined';
t = 'undefined';

% Define tabulation for the output to the command window
tab = '\t';

% Compute the number of nodes in the fluid mesh
numNodes = length(fldMsh.nodes(:, 1));

% Compute the number of degrees of freedom
numDOFs = 3*numNodes;

% Assign a sequential numbering to the system DOFs
DOFNumbering = 1:numDOFs;

% Computation of the force vector
F = zeros(numDOFs, 1);

% Title for the VTK files
title = 'Steady-state stabilized finite element formulation for the 2D incopmpressible Navier Stokes equations';

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs, inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs, prescribedDoFs)) = [];

%% 2. Solve the steady-state nonlinear Navier-Stokes stabilized finite element equation system
[up, FComplete, isConverged, minElSize] = solve_FEMNLinearSystem ...
    (propAnalysis, uSaved, uDotSaved, uDDotSaved, fldMsh, F, ...
    computeBodyForces, propParameters, up, uDot, uDDot, massMtx, ...
    dampMtx, precompStiffMtx, precomResVct, ...
    @computeFEMVMSStabMtxAndVct4NLinear4NSE, DOFNumbering, ...
    freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, uMeshALE, ...
    solve_LinearSystem, propFldDynamics, t, propNLinearAnalysis, ...
    propGaussInt, tab, outMsg);

%% 3. Write out the results into file
if isfield(propOutput, 'isOutput')
    if isa(propOutput.isOutput, 'logical')
        if propOutput.isOutput
            if isfield(propOutput, 'writeOutputToFile')
                if isa(propOutput.writeOutputToFile, 'function_handle')
                    if strcmp(outMsg,'outputEnabled')
                        fprintf('>> Writting out the results to "%s"\n\n',strcat(pathToOutput, caseName, '/'));
                    end
                    DOF4Output = [1:3:numDOFs - 2
                                  2:3:numDOFs - 1
                                  3:3:numDOFs];
                    propOutput.writeOutputToFile(propAnalysis, propNLinearAnalysis, ...
                        propFldDynamics, fldMsh, propParameters, up, uDot, uDDot, DOF4Output, ...
                        caseName, pathToOutput, title, numIterStep);
                else
                    error('Variable propVTK.writeOutputToFile should define a function handle');
                end
            else
                error('Structure propVTK should define variable writeOutputToFile');
            end
        end
    else
        error('Structure propVTK should define boolean isOutput');
    end
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Steady-state nonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('_______________Steady-State Nonlinear Analysis Ended_______________\n');
    fprintf('###################################################################\n\n\n');
end

end