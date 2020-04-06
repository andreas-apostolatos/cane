function [uHistory,minElSize] = solve_FEMTransientAnalysis ...
    (propAnalysis, msh, DOFNumbering, freeDOFs, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, updateInhomDOFs, propALE, computeInitCnds, ...
    computeBodyForceVct, propNBC, computeLoadVct, parameters, ...
    solve_FEMEquationSystem, computeConstantProblemMatrices, ...
    computeMassMtx, computeProblemMatricesSteadyState, computeUpdatedMesh, ...
    solve_LinearSystem, propTransientAnalysis, propNLinearAnalysis, ...
    propIDBC, propGaussInt, propOutput, caseName, pathToOutput, title, ...
    DOF4Output, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the history of the primary variable and the load vector of a
% transient analysis over a field which it discrezed with a classical
% finite element basis.
%
%                            Input :
%                     propAnalysis : Structure containing general 
%                                    information about the analysis,
%                                       .type : The analysis type
%                              msh : Nodes and elements for the FE mesh
%                     DOFNumbering : The global numbering of the DOFs
%                         freeDOFs : The global numbering of the DOFs which 
%                                    are uncontrained
%                          homDOFs : The global numbering of the DOFs over
%                                    which homogeneous Dirichlet boundary
%                                    conditions are applied
%                        inhomDOFs : The global numbering of the DOFs over
%                                    which inhomogeneous Dirichlet boundary
%                                    conditions are applied
%                  valuesInhomDOFs : Prescribed values on the nodes where 
%                                    inhomogeneous boundary conditions are 
%                                    applied
%                  updateInhomDOFs : Function handle to the update of the 
%                                    inhomogeneous Dirichlet boundary 
%                                    conditions
%                          propALE : Structure containing information on 
%                                    the nodes along the ALE boundary,
%                                        .nodes : The sequence of the nodal 
%                                                 coordinates on the ALE 
%                                                 boundary
%                                    .fcthandle : Function handle to the 
%                                                 computation of the ALE 
%                                                 motion
%                  computeInitCnds : Function handle to the initial 
%                                    boundary conditions computation
%              computeBodyForceVct : Function handle to the computation of 
%                                    the body force vector [bx by]'
%          computeBodyForceLoadVct : Function handle to the computation of
%                                    the body force vector
%                          propNBC : Structure containing information on 
%                                    the Neumann boundary conditions,
%                                        .nodes : The nodes where Neumann 
%                                                 boundary conditions are 
%                                                 applied
%                                     .loadType : The type of the load for 
%                                                 each Neumann node
%                                    .fctHandle : The function handle for 
%                                                 each Neumann node for the 
%                                                 computation of the load 
%                                                 vector (these functions 
%                                                 are under the folder 
%                                                 load)
%                   computeLoadVct : Function handle to the computation of
%                                    the boundary load vector
%                       parameters : Flow parameters
%          solve_FEMEquationSystem : Function handle to the solution of the
%                                    equation system resulting out of the
%                                    classical finite element 
%                                    discretization
%   computeConstantProblemMatrices : Handle to function which returns the
%                                    parts of the problem matrices which 
%                                    are time invariant
%                   computeMassMtx : Function handle to the computation of
%                                    the mass matrix of the system
% computeProblemMatricesSteadyState : Function handle to the computation of 
%                                    the problem matrices nessecary for 
%                                    solving the equation system (linear or 
%                                    nonlinear)
%               computeUpdatedMesh : Function handle to the computation of
%                                    the updated mesh location and 
%                                    velocities at the nodes
%               solve_LinearSystem : Function handle to the solver for the 
%                                    linear equation system
%             propIDBC : Structure containing information on the 
%                        inhomogeneous Dirichlet boundary conditions
%            propTransientAnalysis : On the transient analysis :
%                                 .method : The time integration method
%                              .alphaBeta : (parameter for the Bossak 
%                                           scheme)
%                                  .gamma : (parameter for the Bossak 
%                                           scheme)
%                                 .TStart : Start time of the simulation
%                                   .TEnd : End time of the simulation
%                                     .nT : Number of time steps
%                                     .dt : Time step
%              propNLinearAnalysis : On the nonlinear analysis :
%                                 .method : The nonlinear solution scheme
%                                    .eps : The residual tolerance
%                                .maxIter : The maximum number of nonlinear
%                                           iterations
%                         propIDBC : Structure containing information on 
%                                    the inhomogeneous Dirichlet boundary 
%                                    conditions
%                     propGaussInt : On the spatial integration
%                                         .type : 'default', 'user'
%                                   .domainNoGP : Number of Gauss Points 
%                                                 for the domain 
%                                                 integration
%                                 .boundaryNoGP : Number of Gauss Points 
%                                                 for the boundary 
%                                                 integration
%                         .postProcComponent : Which component to plot in 
%                                              the postprocessing
%          propOutput : Structure containing information on writting the
%                       results for postprocessing,
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
%                         caseName : String defining the case name
%                     pathToOutput : Path to where to write out the output 
%                                    file
%                            title : The title to the VTK file
%                       DOF4Output : Array containing the arrangment of the 
%                                    DOFs for printing them out
%                              tab : Tabulation for the output in the 
%                                    command window
%                           outMsg : On printing information during 
%                                    analysis in the command window
%                           Output :
%                         uHistory : The history of the primary field
%                                    throughout the transient analysis
%                        minElSize : The minimum element edge size in the
%                                    mesh
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
% 4. Compute the mass matrix of the problem
%
% 5. Compute part of the residual vector and the stiffness matrix which stays constant throughout the transient analysis
%
% 6. Compute the damping matrix of the problem
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
%    7v. Solve the mesh motion problem and update the mesh node locations and velocities
%
%   7vi. Compute the load vector at the current time step
%
%  7vii. Save the discrete primary field and its first and second time derivatives
%
% 7viii. Solve the equation system
%
%   7ix. Update the time derivatives of the field
%
%    7x. Write out the results into a VTK file or save them into an output variable
% <-
%
%% Function main body
if ~isfield(propTransientAnalysis, 'method')
    error('Undefined time integration method propStrDynamics.method');
end
if ~isfield(propTransientAnalysis, 'T0')
    error('Undefined start time propStrDynamics.TStart');
end
if ~isfield(propTransientAnalysis, 'TEnd')
    error('Undefined end time propStrDynamics.TEnd');
end
if ~isfield(propTransientAnalysis, 'noTimeSteps')
    error('Undefined number of time steps propStrDynamics.noTimeSteps');
end
if ~isfield(propTransientAnalysis, 'dt')
    error('Undefined time step size propStrDynamics.dt');
end
if isfield(propTransientAnalysis, 'damping')
    if ~isfield(propTransientAnalysis.damping, 'method')
        error('Undefined damping method propStrDynamics.damping.method');
    end
    if ~isfield(propTransientAnalysis.damping, 'computeDampMtx')
        error('Undefined damping method propStrDynamics.damping.computeDampMtx');
    else
        if isempty(propTransientAnalysis.damping.computeDampMtx)
            error('String propStrDynamics.damping.computeDampMtx does not define a function');
        end
    end
end
if ~ischar(propALE) && ~isempty(propALE)
    if ~isa(computeUpdatedMesh, 'function_handle')
        error('ALE nodes are found but computeUpdatedMesh is not a function handle');
    end
end
if strcmp(outMsg, 'outputEnabled')
    fprintf(strcat(tab, 'Transient analysis properties\n'));
    fprintf(strcat(tab, '-----------------------------\n\n'));
    fprintf(strcat(tab, 'Time integration method: %s \n'), propTransientAnalysis.method);
    fprintf(strcat(tab, 'Start time = %d \n'), propTransientAnalysis.T0);
    fprintf(strcat(tab, 'End time = %d \n'), propTransientAnalysis.TEnd);
    fprintf(strcat(tab, 'No. time steps = %d \n'), propTransientAnalysis.noTimeSteps);
    fprintf(strcat(tab, 'Time step size = %d \n'), propTransientAnalysis.dt);
    if isfield(propTransientAnalysis, 'damping')
        fprintf(strcat(tab, 'Damping method: %s \n'), propTransientAnalysis.damping.method);
    else
        fprintf(strcat(tab, 'No damping is selected\n'));
    end
    fprintf('\n');
end

%% 0. Read input

% Compute the number of nodes in the fluid mesh
noNodes = length(msh.nodes(:,1));

% Compute the number of degrees of freedom
if strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
    noDOFs = 3*noNodes;
elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
    noDOFs = 4*noNodes;
elseif strcmp(propAnalysis.type, 'planeStress') || strcmp(propAnalysis.type, 'planeStrain')
    noDOFs = 2*noNodes;
else
    error('Select a valid analysis type in analysis.type');
end

% Initialize the scalar concetration history
FHistory = sparse(noDOFs,propTransientAnalysis.noTimeSteps + 1);

% Initialize the output array for the state variable
if isfield(propOutput, 'isOutput')
    if isa(propOutput.isOutput, 'logical')
        if propOutput.isOutput
            if isfield(propOutput, 'writeOutputToFile')
                if isa(propOutput.writeOutputToFile, 'function_handle')
                    uHistory = 'undefined';
                else
                    error('Variable propVTK.writeOutputToFile should define a function handle');
                end
            else
                error('Structure propVTK should define variable writeOutputToFile');
            end
        else
            uHistory = zeros(noDOFs,propTransientAnalysis.noTimeSteps + 1); 
        end
    else
        error('Structure propVTK.isOutput should be a boolean');
    end
else
    error('Structure propVTK should define boolean isOutput');
end

%% 1. Get the initial values for the discrete solution vector and its first and second order rate
if isfield(propOutput, 'VTKResultFile')
    if ischar(propOutput.VTKResultFile)
        [u, uDot, uDDot, numTimeStep] = computeInitCnds...
            (propAnalysis, msh, DOF4Output, parameters, propTransientAnalysis, ...
            propOutput.VTKResultFile, caseName, pathToOutput);
    else
        error('Variable propVTK.VTKResultFile should be a string');
    end
else
    error('Structure propVTK should define string VTKResultFile');
end

%% 2. Initialize the simulation time
t = propTransientAnalysis.dt*numTimeStep;

%% 3. Write out the initial conditions to a file
if propOutput.isOutput
    if isa(propOutput.writeOutputToFile, 'function_handle')
        propOutput.writeOutputToFile(propAnalysis, propNLinearAnalysis, ...
            propTransientAnalysis, msh, parameters, u, uDot, uDDot, ...
            DOF4Output, caseName, pathToOutput, title, numTimeStep);
    else
        uHistory(:,1) = u;
    end
end

%% 4. Compute the mass matrix of the problem
if isa(computeMassMtx, 'function_handle')
    if strcmp(outMsg, 'outputEnabled')
        fprintf(strcat(tab, 'Computing the mass matrix of the system\n'));
        fprintf(strcat(tab, '---------------------------------------\n\n'));
    end
    massMtx = computeMassMtx(propAnalysis, msh, parameters, propGaussInt);
else
    error('Variable computeMassMtx is not defining a function handle as expected');
end

%% 5. Compute part of the residual vector and the stiffness matrix which stays constant throughout the transient analysis
if ~ischar(computeConstantProblemMatrices)
    KConstant = computeConstantProblemMatrices();
end

%% 6. Compute the damping matrix of the problem
if isfield(propTransientAnalysis, 'damping')
    if strcmp(outMsg, 'outputEnabled')
        fprintf(strcat(tab, 'Computing the damping matrix of the system\n'));
        fprintf(strcat(tab, '------------------------------------------\n\n'));
    end 
    dampMtx = propTransientAnalysis.damping.computeDampMtx...
        (propAnalysis, msh, parameters, propGaussInt);
else
    dampMtx = 'undefined';
end

%% 7. Loop over all the time instances of the simulation
if strcmp(outMsg, 'outputEnabled')
    fprintf('\tLooping over all the time steps \n\t------------------------------- \n\n');
end
while t < propTransientAnalysis.TEnd && numTimeStep < propTransientAnalysis.noTimeSteps
    %% 7i. Update the simulation time
    t = t + propTransientAnalysis.dt;
    
	%% 7ii. Update the time counter
    numTimeStep = numTimeStep + 1;
    
    %% 7iii. Preamble of the time stepping iterations
    if strcmp(outMsg, 'outputEnabled')
        if ~ischar(propNLinearAnalysis)
            msgTS = sprintf(strcat(tab, '\tTime step %d/%d at real time %d seconds with dt=%d and maxNoNRIter=%d \n \n'), ...
                numTimeStep, propTransientAnalysis.noTimeSteps, t, propTransientAnalysis.dt, propNLinearAnalysis.maxIter);
        else
            msgTS = sprintf(strcat(tab, '\tTime step %d/%d at real time %d seconds with dt=%d \n \n'), ...
                numTimeStep, propTransientAnalysis.noTimeSteps, t, propTransientAnalysis.dt);
        end
        fprintf(msgTS);
    end
    
    %% 7iv. Update the values of the inhomogeneous Dirichlet boundary conditions
    if isa(updateInhomDOFs, 'function_handle')
        valuesInhomDOFs = updateInhomDOFs(valuesInhomDOFs, propIDBC);
    end
    
    %% 7v. Solve the mesh motion problem and update the mesh node locations and velocities
    nodesSaved = msh.nodes;
    if ~ischar(propALE) && ~isempty(propALE)
        [msh, uMeshALE, inhomDOFs, valuesInhomDOFs] = ...
            computeUpdatedMesh...
            (msh, homDOFs, inhomDOFs, valuesInhomDOFs, nodesSaved, ...
            propALE, solve_LinearSystem, propTransientAnalysis, t);
    else
        uMeshALE = 'undefined';
    end
    
    %% 7vi. Compute the load vector at the current time step
    if ~ischar(computeLoadVct)
        FHistory(:, numTimeStep) = computeLoadVct ... 
            (msh, propAnalysis, propNBC, t, propGaussInt, '');
    elseif strcmp(computeLoadVct,'undefined')
        FHistory(:, numTimeStep) = zeros(noDOFs, 1);
    end
        
    %% 7vii. Save the discrete primary field and its first and second time derivatives
    uSaved = u;
    uDotSaved = uDot;
    uDDotSaved = uDDot;
    
    %% 7viii. Solve the equation system
    [u, ~, isConverged, minElSize] = solve_FEMEquationSystem...
        (propAnalysis, uSaved, uDotSaved, uDDotSaved, msh, ...
        FHistory(:, numTimeStep), computeBodyForceVct, parameters, u, ...
        uDot, uDDot, massMtx, dampMtx, computeProblemMatricesSteadyState, ...
        DOFNumbering, freeDOFs, homDOFs, inhomDOFs, valuesInhomDOFs, ...
        uMeshALE, solve_LinearSystem, propTransientAnalysis, t, ...
        propNLinearAnalysis, propGaussInt, strcat(tab,'\t'), outMsg);
    
    %% 7ix. Update the time derivatives of the field
    if isConverged
        if isa(propTransientAnalysis.computeUpdatedVct, 'function_handle')
            [uDot, uDDot] = propTransientAnalysis.computeUpdatedVct ...
                (u, uSaved, uDotSaved, uDDotSaved, propTransientAnalysis);
        else
            error('Function handle propTransientAnalysis.computeUpdatedVct undefined');
        end
    end
    
    %% 7x. Write out the results into a VTK file or save them into an output variable
    if isa(propOutput.writeOutputToFile, 'function_handle')
        propOutput.writeOutputToFile...
            (propAnalysis, propNLinearAnalysis, propTransientAnalysis, ...
            msh, parameters, u, uDot, uDDot, DOF4Output, caseName, ...
            pathToOutput, title, numTimeStep);
    else
        uHistory(:,numTimeStep + 1) = u;
    end
end

end