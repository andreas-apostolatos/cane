function [uHistory_A, uHistory_B, resVctHistory_A, resVctHistory_B, ...
    minElSize_A, minElSize_B] = solve_FEMCoupledTransientAnalysis ...
    (propAnalysis_A, propAnalysis_B, mesh_A, mesh_B, DOFNumbering_A, DOFNumbering_B,...
    freeDOFs_A, freeDOFs_B, homDOFs_A, homDOFs_B, inhomDOFs_A, inhomDOFs_B, ...
    valuesInhomDOFs_A, valuesInhomDOFs_B, updateInhomDOFs_A, updateInhomDOFs_B, ...
    propALE_A, propALE_B, computeInitCnds_A, computeInitCnds_B, ...
    computeBodyForceVct_A, computeBodyForceVct_B, propNBC_A, propNBC_B, ...
    computeLoadVct_A, computeLoadVct_B, propParameters_A, propParameters_B, ...
    solve_FEMCoupledSystem, solve_FEMEquationSystem_A, solve_FEMEquationSystem_B, ...
    computeConstMtx_A, computeConstMtx_B, computeMassMtx_A, computeMassMtx_B, ...
    computeMtxSteadyState_A, computeMtxSteadyState_B, ...
    computeUpdatedMesh_A, computeUpdatedMesh_B, solve_LinearSystem_A, ...
    solve_LinearSystem_B, propTransientAnalysis_A, propTransientAnalysis_B, ...
    propNLinearAnalysis_A, propNLinearAnalysis_B, propCoupling, propIDBC_A, propIDBC_B, ...
    propPostProc_A, propPostProc_B, propGaussInt_A, propGaussInt_B, propOutput_A, ...
    propOutput_B, caseName_A, caseName_B, pathToOutput_A, pathToOutput_B, title_A, ...
    title_B, DOF4Output_A, DOF4Output_B, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the solution of a coupled transient problem together with the
% residual history assuming mathing time spans and time discretizations.
%
%                         Input :
% propAnalysis_A,propAnalysis_B : Structure containing general information 
%                                 on the analysis A and B,
%                                       .type : Analysis type
%                 mesh_A,mesh_B : Structures containing information on the 
%                                 meshes corresponding to analyses A and B,
%                             .initialNodes : The initial nodes of the mesh
%                                    .nodes : The updated nodes of the mesh
%                                 .elements : The elements of the mesh
% DOFNumbering_A,DOFNumbering_B : DOF numbering for the analyses A and B
%         freeDOFs_A,freeDOFs_B : Vectors containing the free DOFs of the
%                                 analyses A and B
%           homDOFs_A,homDOFs_B : Vectors containing the global numbering 
%                                 of the DOFs where homogeneneous Dirichlet
%                                 boundary conditions are applied
%       inhomDOFs_A,inhomDOFs_B : Vectors containing the global numbering 
%                                 of the DOFs where inhomogeneneous 
%                                 Dirichlet boundary conditions are applied
%             valuesInhomDOFs_A : Prescribed values of the DOFs on the 
%                                 inhomogeneous Dirichlet boundary of
%                                 analysis A
%             valuesInhomDOFs_B : Prescribed values of the DOFs on the 
%                                 inhomogeneous Dirichlet boundary of
%                                 analysis B
%             updateInhomDOFs_A : Function handle to the update of the
%                                 inhomogeneous Dirichlet boundary
%                                 conditions for analysis A
%             updateInhomDOFs_B : Function handle to the update of the
%                                 inhomogeneous Dirichlet boundary
%                                 conditions for analysis B
%           propALE_A,propALE_B : Structures containing information about 
%                                 the ALE motion of the fluid along ALE
%                                 interfaces of analyses A and B,
%                                   .nodes : IDs of the fluid nodes on the 
%                                            ALE boundary
%                               .fctHandle : Function handles to the
%                                            computation of the ALE motion
%                                            at each node on the ALE
%                                            boundary
%                                            boundary
%                                  .isFree : Vector of flags indicating
%                                            whether the fluid motion is
%                                            dictated by the ALE motion at
%                                            each node on the ALE boundary
%             computeInitCnds_A : Function handle to the computation of the
%                                 initial conditions for analysis A
%             computeInitCnds_B : Function handle to the computation of the
%                                 initial conditions for analysis B
%         computeBodyForceVct_A : Function handle to the computation of the
%                                 body force vector for analysis A
%         computeBodyForceVct_B : Function handle to the computation of the
%                                 body force vector for analysis B
%           propNBC_A,propNBC_B : Structures containing information on the
%                                 Neumann boundary conditions for analyses
%                                 A and B
%              computeLoadVct_A : Function handle to the computation of the
%                                 boundary load vector for analysis A
%              computeLoadVct_B : Function handle to the computation of the
%                                 boundary load vector for analysis B
%              propParameters_A : Structure containing information on the
%                                 parameters of analysis A
%              propParameters_B : Structure containing information on the
%                                 parameters of analysis B
%        solve_FEMCoupledSystem : Function handle to the solution algorithm
%                                 for the coupled system
%     solve_FEMEquationSystem_A : solve_FEMLinearSystem or 
%                                 solver_FEMNLinearSystem
%     solve_FEMEquationSystem_B : solve_FEMLinearSystem or 
%                                 solver_FEMNLinearSystem
%             computeConstMtx_A : Function handle to the computation of the
%                                 constant part of the left-hand side
%                                 matrix for analysis A
%             computeConstMtx_B : Function handle to the computation of the
%                                 constant part of the left-hand side
%                                 matrix for analysis B
%              computeMassMtx_A : Function handle to the computation of the
%                                 mass matrix of system A
%              computeMassMtx_B : Function handle to the computation of the
%                                 mass matrix of system B
%       computeMtxSteadyState_A : Function handle to the computation of the
%                                 steady-state left-hand side matrix of
%                                 system A
%       computeMtxSteadyState_B : Function handle to the computation of the
%                                 steady-state left-hand side matrix of
%                                 system B
%          computeUpdatedMesh_A : Function handle to the update of the mesh
%                                 for system A
%          computeUpdatedMesh_B : Function handle to the update of the mesh
%                                 for system B
%          solve_LinearSystem_A : Function handle to the linear equation
%                                 system solver for system A
%          solve_LinearSystem_B : Function handle to the linear equation
%                                 system solver for system B
%       propTransientAnalysis_A : Structure containing information on the 
%                                 structural dynamics of system A
%                              .method : The time integration method
%                           .computeProblemMtrcsTransient : Function handle
%                                                           to the
%                                                           computation of
%                                                           the transient
%                                                           problem 
%                                                           matrices
%                                      .computeUpdatedVct : Function handle
%                                                           to the update
%                                                           of the solution 
%                                                           and its time
%                                                           derivatives
%                                                 .TStart : Start time of 
%                                                           the simulation
%                                                   .TEnd : End time of the 
%                                                           simulation
%                                            .noTimeSteps : Number of time 
%                                                           steps
%                                                     .dt : Time step size
%       propTransientAnalysis_B : Structure containing information on the 
%                                 structural dynamics of system B
%                              .method : The time integration method
%                           .computeProblemMtrcsTransient : Function handle
%                                                           to the
%                                                           computation of
%                                                           the transient
%                                                           problem 
%                                                           matrices
%                                      .computeUpdatedVct : Function handle
%                                                           to the update
%                                                           of the solution 
%                                                           and its time
%                                                           derivatives
%                                                 .TStart : Start time of 
%                                                           the simulation
%                                                   .TEnd : End time of the 
%                                                           simulation
%                                            .noTimeSteps : Number of time 
%                                                           steps
%                                                     .dt : Time step size
%         propNLinearAnalysis_A : Structure containing information on the 
%                                 nonlinear analysis for system A
%                                .method : The nonlinear solution scheme
%                                   .eps : The residual tolerance
%                               .maxIter : The maximum number of nonlinear
%                                          iterations
%         propNLinearAnalysis_B : Structure containing information on the 
%                                 nonlinear analysis for system B
%                                .method : The nonlinear solution scheme
%                                   .eps : The residual tolerance
%                               .maxIter : The maximum number of nonlinear
%                                          iterations
%                  propCoupling : Structure containing information on the 
%                                 coupling algorithm,
%                               .relaxation : Under-relaxation factor
%                                      .tol : Relative residual tolerance
%                                             on the displacement field
%         propIDBC_A,propIDBC_B : Structures containing information on the
%                                 update of the inhomogeneous Dirichlet
%                                 boundary conditions for systems A and B
% propPostProc_A,propPostProc_B : Structures containing information on the
%                                 computation of postprocessing resultants 
%                                 for analyses A and B,
%                               .nameDomain : Name of the domain onto which
%                                             postprocessing resultants are
%                                             to be computed
%                              .nodesDomain : IDs of the nodes which are on
%                                             the selected domain
%                          .computePostProc : Function handle to the
%                                             computation of the desirable 
%                                             resultant
% propGaussInt_A,propGaussInt_B : Structure containing information on the
%                                 numerical quadrature for analyses A and B
%                              .type : 'default', 'user'
%                        .domainNoGP : Number of Gauss Points for the 
%                                      domain integration
%                      .boundaryNoGP : Number of Gauss Points for the
%                                      boundary integration
%     propOutput_A,propOutput_B : Structure containing information on
%                                 writting the results of solutions to the 
%                                 A and B systems in files,
%                                  .isOutput : Flag on whether the results 
%                                              to be written out
%                         .writeOutputToFile : Function handle to the
%                                              writting out of the results
%                                              in a VTK format
%                             .VTKResultFile : Specifies the name of the
%                                              VTK result file from which
%                                              the simulation to be
%                                              restarted. If it is
%                                              specified as 'undefined'
%                                              the simulation starts from
%                                              time TStart
%         caseName_A,caseName_B : Name of the cases for analyses A and B
% pathToOutput_A,pathToOutput_B : Paths to where the output files are
%                                 stored
%               title_A,title_B : Titles of the analyses in the
%                                 corresponding output files
%     DOF4Output_A,DOF4Output_B : Arrays containing the re-organisation of
%                                 the DOFs suitable for the output files
%                           tab : Tabulation for the output messages on the
%                                 command window
%                        outMsg : Enables outputting information on the
%                                 command window if it is chosen as
%                                 'outputEnabled'
%
%                        Output :
%         uHistory_A,uHistory_B : Arrays containing the solution of systems 
%                                 A and B for all time steps or 'undefined' 
%                                 if the results are written out to file
%               resVctHistory_A : Array containing the history of the full
%                                 right-hand side residual vector for
%                                 system A or 'undefined' if the results 
%                                 are written out to file
%               resVctHistory_B : Array containing the history of the full
%                                 right-hand side residual vector for
%                                 system B or 'undefined' if the results 
%                                 are written out to file
%       minElSize_A,minElSize_B : The minimum element edge lengths for the
%                                 meshes corresponding to analyses A and B
%
% Function layout :
%
% 0. Read input
%
% 1. Get the initial values for the discrete solution vector and its first and second order rate
%
% 2. Initialize the simulation time
%
% 3. Write out the initial conditions to a file
%
% 4. Compute the mass matrices of the systems
%
% 5. Compute part of the residual vector and the stiffness matrix which stays constant throughout the transient analysis
%
% 6. Compute the damping matrices of the systems
%
% 7. Loop over all the time instances of the simulation
% ->
%    7i. Update the simulation time
%
%   7ii. Update the time counter
%
%  7iii. Preamble of the time stepping iterations
%
%   7iv. Update the values of the inhomogeneous Dirichlet boundary conditions
%
%    7v. Compute the load vector at the current time step
%
%   7vi. Save the discrete primary field and its first and second time derivatives
%
%  7vii. Solve the coupled system
%
% 7viii. Write out the results into a VTK file or save them into an output file
% <-
%
%% Function main body
if ~isfield(propTransientAnalysis_A, 'method') || ~isfield(propTransientAnalysis_B, 'method')
    error('Undefined time integration method propStrDynamics.method');
end
if ~isfield(propTransientAnalysis_A, 'T0') || ~isfield(propTransientAnalysis_B, 'T0')
    error('Undefined start time propStrDynamics.TStart');
end
if ~isfield(propTransientAnalysis_A, 'TEnd') || ~isfield(propTransientAnalysis_B, 'TEnd')
    error('Undefined end time propStrDynamics.TEnd');
end
if ~isfield(propTransientAnalysis_A, 'noTimeSteps') || ~isfield(propTransientAnalysis_B, 'noTimeSteps')
    error('Undefined number of time steps propStrDynamics.noTimeSteps');
end
if ~isfield(propTransientAnalysis_A, 'dt') || ~isfield(propTransientAnalysis_B, 'dt')
    error('Undefined time step size propStrDynamics.dt');
end
if isfield(propTransientAnalysis_A, 'damping')
    if ~isfield(propTransientAnalysis_A.damping, 'method')
        error('Undefined damping method propStrDynamics.damping.method');
    end
    if ~isfield(propTransientAnalysis_A.damping, 'computeDampMtx')
        error('Undefined damping method propStrDynamics.damping.computeDampMtx');
    else
        if isempty(propTransientAnalysis_A.damping.computeDampMtx)
            error('String propStrDynamics.damping.computeDampMtx does not define a function');
        end
    end
end
if isfield(propTransientAnalysis_B, 'damping')
    if ~isfield(propTransientAnalysis_B.damping, 'method')
        error('Undefined damping method propStrDynamics.damping.method');
    end
    if ~isfield(propTransientAnalysis_B.damping, 'computeDampMtx')
        error('Undefined damping method propStrDynamics.damping.computeDampMtx');
    else
        if isempty(propTransientAnalysis_B.damping.computeDampMtx)
            error('String propStrDynamics.damping.computeDampMtx does not define a function');
        end
    end
end
if strcmp(outMsg, 'outputEnabled')
    fprintf(strcat(tab, 'Transient analysis properties for problem %s\n'), propAnalysis_A.type);
    fprintf(strcat(tab, '--------------------------------------------\n\n'));
    fprintf(strcat(tab, 'Time integration method: %s \n'), propTransientAnalysis_A.method);
    fprintf(strcat(tab, 'Start time = %d \n'), propTransientAnalysis_A.T0);
    fprintf(strcat(tab, 'End time = %d \n'), propTransientAnalysis_A.TEnd);
    fprintf(strcat(tab, 'No. time steps = %d \n'), propTransientAnalysis_A.noTimeSteps);
    fprintf(strcat(tab, 'Time step size = %d \n'), propTransientAnalysis_A.dt);
    if isfield(propTransientAnalysis_A, 'damping')
        fprintf(strcat(tab, 'Damping method: %s \n'), propTransientAnalysis_A.damping.method);
    else
        fprintf(strcat(tab, 'No damping is selected\n'));
    end
    fprintf('\n');
    fprintf(strcat(tab, 'Transient analysis properties for problem %s\n'), propAnalysis_B.type);
    fprintf(strcat(tab, '--------------------------------------------\n\n'));
    fprintf(strcat(tab, 'Time integration method: %s \n'), propTransientAnalysis_B.method);
    fprintf(strcat(tab, 'Start time = %d \n'), propTransientAnalysis_B.T0);
    fprintf(strcat(tab, 'End time = %d \n'), propTransientAnalysis_B.TEnd);
    fprintf(strcat(tab, 'No. time steps = %d \n'), propTransientAnalysis_B.noTimeSteps);
    fprintf(strcat(tab, 'Time step size = %d \n'), propTransientAnalysis_B.dt);
    if isfield(propTransientAnalysis_B, 'damping')
        fprintf(strcat(tab, 'Damping method: %s \n'), propTransientAnalysis_B.damping.method);
    else
        fprintf(strcat(tab, 'No damping is selected\n'));
    end
    fprintf('\n');
end
if ~ischar(propALE_A) && ~isempty(propALE_A)
    if ~isa(computeUpdatedMesh_A, 'function_handle')
        error('ALE nodes are found but computeUpdatedMesh is not a function handle for problem A');
    end
end
if ~ischar(propALE_B) && ~isempty(propALE_B)
    if ~isa(computeUpdatedMesh_B, 'function_handle')
        error('ALE nodes are found but computeUpdatedMesh is not a function handle for problem B');
    end
end

%% 0. Read input

% Create a common time span and time discretization
if ~strcmp(propTransientAnalysis_A.timeDependence, propTransientAnalysis_B.timeDependence)
    error('Time dependencies for the problems A and B do not match');
end
if propTransientAnalysis_A.T0 ~= propTransientAnalysis_B.T0
    error('Start times for problems A and B do not match');
end
if propTransientAnalysis_A.TEnd ~= propTransientAnalysis_B.TEnd
    error('End times for the problems A and B do not match');
end
if propTransientAnalysis_A.noTimeSteps ~= propTransientAnalysis_B.noTimeSteps
    error('Number of time steps for problems A and B do not match');
end
propTransientAnalysis.T0 = propTransientAnalysis_A.T0;
propTransientAnalysis.TEnd = propTransientAnalysis_A.TEnd;
propTransientAnalysis.noTimeSteps = propTransientAnalysis_A.noTimeSteps;
propTransientAnalysis.dt = propTransientAnalysis_A.dt;

% Compute the number of nodes and DOFs

% BVP A
% -----

numNodes_A = length(mesh_A.nodes(:,1));
if strcmp(propAnalysis_A.type, 'NAVIER_STOKES_2D')
    numDOFs_A = 3*numNodes_A;
elseif strcmp(propAnalysis_A.type, 'NAVIER_STOKES_3D')
    numDOFs_A = 4*numNodes_A;
elseif strcmp(propAnalysis_A.type, 'planeStress') || strcmp(propAnalysis_A.type, 'planeStrain')
    numDOFs_A = 2*numNodes_A;
elseif strcmp(propAnalysisA.type, 'THERMAL_CONDUCTION_2D') || strcmp(propAnalysis_A.type, 'SDOF')
    numDOFs_A = numNodes_A;
else
    error('Select a valid analysis type in analysis_A.type');
end

% BVP B
% -----

numNodes_B = length(mesh_B.nodes(:,1));
if strcmp(propAnalysis_B.type, 'NAVIER_STOKES_2D')
    numDOFs_B = 3*numNodes_B;
elseif strcmp(propAnalysis_B.type, 'NAVIER_STOKES_3D')
    numDOFs_B = 4*numNodes_B;
elseif strcmp(propAnalysis_B.type, 'planeStress') || strcmp(propAnalysis_B.type, 'planeStrain')
    numDOFs_B = 2*numNodes_B;
elseif strcmp(propAnalysis_B.type, 'THERMAL_CONDUCTION_2D') || strcmp(propAnalysis_B.type, 'SDOF')
    numDOFs_B = numNodes_B;
else
    error('Select a valid analysis type in analysis_B.type');
end

% Check the function handle for the computation of the updated time
% derivatives of the field
if ~isa(propTransientAnalysis_A.computeUpdatedVct, 'function_handle') || ...
        ~isa(propTransientAnalysis_B.computeUpdatedVct, 'function_handle')
    error('Function handle propTransientAnalysis.computeUpdatedVct undefined');
end

% Initialize the output arrays according to whether the values are saved in
% the workspace or not

% BVP A
% -----

isWriteOutVariables_A = false;
if isfield(propOutput_A, 'isOutput')
    if isa(propOutput_A.isOutput, 'logical')
        if propOutput_A.isOutput
            if isfield(propOutput_A, 'writeOutputToFile')
                if isa(propOutput_A.writeOutputToFile, 'function_handle')
                    isWriteOutVariables_A = true;
                else
                    error('Variable propOutput_A.writeOutputToFile should define a function handle');
                end
            else
                error('Structure propOutput_A should define variable writeOutputToFile');
            end
        end
    else
        error('Structure propOutput_A.isOutput should be a boolean');
    end
else
    error('Structure propOutput_A should define boolean isOutput');
end
if isWriteOutVariables_A
    uHistory_A = 'undefined';
    resVctHistory_A = 'undefined';
else
    uHistory_A = zeros(numDOFs_A, propTransientAnalysis.noTimeSteps + 1);
    resVctHistory_A = zeros(numDOFs_A, propTransientAnalysis.noTimeSteps + 1);
end

% Initialize interface force vectors
FInterface_A = zeros(numDOFs_A, 1);
FInterface_B = zeros(numDOFs_B, 1);

% BVP B
% -----

isWriteOutVariables_B = false;
if isfield(propOutput_B, 'isOutput')
    if isa(propOutput_B.isOutput, 'logical')
        if propOutput_B.isOutput
            if isfield(propOutput_B, 'writeOutputToFile')
                if isa(propOutput_B.writeOutputToFile, 'function_handle')
                    isWriteOutVariables_B = true;
                else
                    error('Variable propOutput_B.writeOutputToFile should define a function handle');
                end
            else
                error('Structure propOutput_B should define variable writeOutputToFile');
            end
        end
    else
        error('Structure propOutput_B.isOutput should be a boolean');
    end
else
    error('Structure propOutput_B should define boolean isOutput');
end
if isWriteOutVariables_B
    uHistory_B = 'undefined';
    resVctHistory_B = 'undefined';
else
    uHistory_B = zeros(numDOFs_B, propTransientAnalysis.noTimeSteps + 1);
    resVctHistory_B = zeros(numDOFs_B, propTransientAnalysis.noTimeSteps + 1);
end

%% 1. Get the initial values for the discrete solution vector and its first and second order rate

% BVP A
% -----

if isfield(propOutput_A, 'VTKResultFile')
    if ischar(propOutput_A.VTKResultFile)
        [u_A, uDot_A, uDDot_A, numTimeStep_A] = computeInitCnds_A ...
            (propAnalysis_A, mesh_A, DOF4Output_A, propParameters_A, ...
            propTransientAnalysis_A, propOutput_A.VTKResultFile, ...
            caseName_A, pathToOutput_A);
    else
        error('Variable propOutput_A.VTKResultFile should be a string');
    end
else
    error('Structure propOutput_A should define string VTKResultFile');
end

% BVP B
% -----

if isfield(propOutput_B, 'VTKResultFile')
    if ischar(propOutput_B.VTKResultFile)
        [u_B, uDot_B, uDDot_B, numTimeStep_B] = computeInitCnds_B ...
            (propAnalysis_B, mesh_B, DOF4Output_B, propParameters_B, ...
            propTransientAnalysis_B, propOutput_B.VTKResultFile, ...
            caseName_B, pathToOutput_B);
    else
        error('Variable propOutput_B.VTKResultFile should be a string');
    end
else
    error('Structure propOutput_B should define string VTKResultFile');
end

% Check if the time step counters match
if numTimeStep_A ~= numTimeStep_B
    error('Start time step IDs are not matching');
else
    numTimeStep = numTimeStep_A; 
end

%% 2. Initialize the simulation time
t = propTransientAnalysis.dt*numTimeStep;

%% 3. Write out the initial conditions to a file

% Advance time step counter
numTimeStep = numTimeStep + 1;

% BVP A
% -----

if isWriteOutVariables_A
    propOutput_A.writeOutputToFile ...
        (propAnalysis_A, propNLinearAnalysis_A, propTransientAnalysis_A, ...
        mesh_A, propParameters_A, u_A, uDot_A, uDDot_A, DOF4Output_A, ...
        caseName_A, pathToOutput_A, title_A, numTimeStep);
else
    uHistory_A(:,numTimeStep) = u_A;
end

% BVP B
% -----

if isWriteOutVariables_B
    propOutput_B.writeOutputToFile ...
        (propAnalysis_B, propNLinearAnalysis_B, propTransientAnalysis_B, ...
        mesh_B, propParameters_B, u_B, uDot_B, uDDot_B, DOF4Output_B, ...
        caseName_B, pathToOutput_B, title_B, numTimeStep);
else
    uHistory_B(:,numTimeStep) = u_B;
end

%% 4. Compute the mass matrices of the systems

% BVP A
% -----

if isa(computeMassMtx_A, 'function_handle')
    if strcmp(outMsg, 'outputEnabled')
        fprintf(strcat(tab, 'Computing the mass matrix of %s\n'), ...
            propAnalysis_A.type);
        fprintf(strcat(tab, '-------------------------------------\n\n'));
    end
    massMtx_A = computeMassMtx_A ...
        (propAnalysis_A, mesh_A, propParameters_A, propGaussInt_A);
else
    error('Variable computeMassMtx_A is not defining a function handle as expected');
end

% BVP B
% -----

if isa(computeMassMtx_A, 'function_handle')
    if strcmp(outMsg, 'outputEnabled')
        fprintf(strcat(tab, 'Computing the mass matrix of %s\n'), ...
            propAnalysis_B.type);
        fprintf(strcat(tab, '-------------------------------------\n\n'));
    end
    massMtx_B = computeMassMtx_B ...
        (propAnalysis_B, mesh_B, propParameters_B, propGaussInt_B);
else
    error('Variable computeMassMtx_B is not defining a function handle as expected');
end

%% 5. Compute part of the residual vector and the stiffness matrix which stays constant throughout the transient analysis

% BVP A
% -----

if ~ischar(computeConstMtx_A)
    [precompStiffMtx_A, precomResVct_A] = computeConstantProblemMatrices();
else
    precompStiffMtx_A = 'undefined';
    precomResVct_A = 'undefined';
end

% BVP B
% -----

if ~ischar(computeConstMtx_B)
    [precompStiffMtx_B, precomResVct_B] = computeConstantProblemMatrices();
else
    precompStiffMtx_B = 'undefined';
    precomResVct_B = 'undefined';
end

%% 6. Compute the damping matrices of the systems

% BVP A
% -----

if isfield(propTransientAnalysis_A, 'damping')
    if strcmp(outMsg, 'outputEnabled')
        fprintf(strcat(tab, 'Computing the damping matrix of %s\n'), ...
            propAnalysis_A.type);
        fprintf(strcat(tab, '-----------------------------------------\n\n'));
    end 
    dampMtx_A = propTransientAnalysis_A.damping.computeDampMtx ...
        (propAnalysis_A, mesh_A, propParameters_A, propGaussInt_A);
else
    dampMtx_A = 'undefined';
end

% BVP B
% -----

if isfield(propTransientAnalysis_B, 'damping')
    if strcmp(outMsg, 'outputEnabled')
        fprintf(strcat(tab, 'Computing the damping matrix of %s\n'), ...
            propAnalysis_B.type);
        fprintf(strcat(tab, '-----------------------------------------\n\n'));
    end 
    dampMtx_B = propTransientAnalysis_A.damping.computeDampMtx ...
        (propAnalysis_B, mesh_B, propParameters_B, propGaussInt_B);
else
    dampMtx_B = 'undefined';
end

%% 7. Loop over all the time instances of the simulation
if strcmp(outMsg, 'outputEnabled')
    fprintf('\tLooping over all the time steps \n\t------------------------------- \n\n');
end
while t < propTransientAnalysis.TEnd && numTimeStep <= propTransientAnalysis.noTimeSteps
    %% 7i. Update the simulation time
    t = t + propTransientAnalysis.dt;
    
	%% 7ii. Update the time counter
    numTimeStep = numTimeStep + 1;
    
    %% 7iii. Preamble of the time stepping iterations
    if strcmp(outMsg, 'outputEnabled')
        msgTS = sprintf(strcat(tab, '\tTime step %d/%d at real time %d seconds with dt=%d \n \n'), ...
            numTimeStep - 1, propTransientAnalysis.noTimeSteps, t, ...
            propTransientAnalysis.dt);
        fprintf(msgTS);
    end
    
    %% 7iv. Update the values of the inhomogeneous Dirichlet boundary conditions
    
    % BVP A
    % -----
    
    if isa(updateInhomDOFs_A, 'function_handle')
        valuesInhomDOFs_A = updateInhomDOFs_A ...
            (mesh_A, propParameters_A, t);
    end
    
    % BVP B
    % -----
    
    if isa(updateInhomDOFs_B, 'function_handle')
        valuesInhomDOFs_B = updateInhomDOFs_B ...
            (mesh_B, propParameters_B, t);
    end

    %% 7v. Compute the load vector at the current time step
    
    % BVP A
    % -----
    
    if isa(computeLoadVct_A, 'function_handle')
        F_A = computeLoadVct_A ... 
            (mesh_A, propAnalysis_A, propNBC_A, t, propGaussInt_A, '');
    else
        F_A = zeros(numDOFs_A, 1);
    end
    
    % BVP B
    % -----
    
    if isa(computeLoadVct_B, 'function_handle')
        F_B = computeLoadVct_B ... 
            (mesh_B, propAnalysis_B, propNBC_B, t, propGaussInt_B, '');
    else
        F_B = zeros(numDOFs_B, 1);
    end
        
    %% 7vi. Save the discrete primary field and its first and second time derivatives
    
    % BVP A
    % -----
    
    uSaved_A = u_A;
    uDotSaved_A = uDot_A;
    uDDotSaved_A = uDDot_A;
    
    % BVP B
    % -----
    
    uSaved_B = u_B;
    uDotSaved_B = uDot_B;
    uDDotSaved_B = uDDot_B;
    
    %% 7vii. Solve the coupled system
    [u_A, u_B, uDot_A, uDot_B, uDDot_A, uDDot_B, resVct_A, resVct_B, ...
        mesh_A, mesh_B, FInterface_A, FInterface_B, minElSize_A, minElSize_B] = ...
        solve_FEMCoupledSystem ...
        (propAnalysis_A, propAnalysis_B, uSaved_A, uSaved_B, ...
        uDotSaved_A, uDotSaved_B, uDDotSaved_A, uDDotSaved_B, ...
        mesh_A, mesh_B, F_A, F_B, FInterface_A, FInterface_B, ...
        computeBodyForceVct_A, computeBodyForceVct_B, propParameters_A, ...
        propParameters_B, u_A, u_B, uDot_A, uDot_B, uDDot_A, uDDot_B, ...
        massMtx_A, massMtx_B, dampMtx_A, dampMtx_B, precompStiffMtx_A, ...
        precompStiffMtx_B, precomResVct_A, precomResVct_B, ...
        computeMtxSteadyState_A, computeMtxSteadyState_B, ...
        solve_FEMEquationSystem_A, solve_FEMEquationSystem_B, ...
        DOFNumbering_A, DOFNumbering_B, freeDOFs_A, freeDOFs_B, homDOFs_A, ...
        homDOFs_B, inhomDOFs_A, inhomDOFs_B, valuesInhomDOFs_A, valuesInhomDOFs_B, ...
        propALE_A, propALE_B, computeUpdatedMesh_A, computeUpdatedMesh_B, ...
        solve_LinearSystem_A, solve_LinearSystem_B, ...
        propTransientAnalysis_A, propTransientAnalysis_B, t, ...
        propNLinearAnalysis_A, propNLinearAnalysis_B, propCoupling, ...
        propPostProc_A, propPostProc_B, propGaussInt_A, propGaussInt_B, ...
        strcat(tab,'\t'), outMsg);
    
    %% 7viii. Write out the results into a VTK file or save them into an output file
    
    % BVP A
    % -----
    
    if isWriteOutVariables_A
        propOutput_A.writeOutputToFile ...
            (propAnalysis_A, propNLinearAnalysis_A, propTransientAnalysis_A, ...
            mesh_A, propParameters_A, u_A, uDot_A, uDDot_A, DOF4Output_A, caseName_A, ...
            pathToOutput_A, title_A, numTimeStep);
    else
        uHistory_A(:, numTimeStep) = u_A;
        resVctHistory_A(:, numTimeStep) = resVct_A;
    end
    
    % BVP B
    % -----
    
    if isWriteOutVariables_B
        propOutput_B.writeOutputToFile ...
            (propAnalysis_B, propNLinearAnalysis_B, propTransientAnalysis_B, ...
            mesh_B, propParameters_B, u_B, uDot_B, uDDot_B, DOF4Output_B, caseName_B, ...
            pathToOutput_B, title_B, numTimeStep);
    else
        uHistory_B(:, numTimeStep) = u_B;
        resVctHistory_B(:, numTimeStep) = resVct_B;
    end
end

end