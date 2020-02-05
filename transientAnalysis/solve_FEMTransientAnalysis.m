function [uHistory,minElSize] = solve_FEMTransientAnalysis...
    (analysis,msh,DOFNumbering,freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,...
    nodesALE,computeInitCnds,VTKResultFile,computeBodyForceVct,NBC,...
    computeLoadVct,parameters,solve_FEMEquationSystem,...
    computeConstantProblemMatrices,computeMassMtx,...
    computeProblemMatricesSteadyState,computeUpdatedMesh,...
    solve_LinearSystem,propTransientAnalysis,propNLinearAnalysis,...
    propGaussInt,caseName,pathToOutput,title,DOF4Output,writeOutputToVTK,...
    tab,outMsg)
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
%                         analysis : Information on the analysis type
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
%                         nodesALE : The nodes on the ALE boundary:
%                                         .nodes : The sequence of the 
%                                                  nodal coordinates on the 
%                                                  ALE boundary
%                                     .fcthandle : Function handle to the 
%                                                  computation of the ALE 
%                                                  motion
%                  computeInitCnds : Function handle to the initial 
%                                    boundary conditions computation
%                    VTKResultFile : The name of the result file in the 
%                                    output folder where to get the initial 
%                                    conditions for the transient 
%                                    simulation
%              computeBodyForceVct : Function handle to the computation of 
%                                    the body force vector [bx by]'
%          computeBodyForceLoadVct : Function handle to the computation of
%                                    the body force vector
%                              NBC : On the Neumann boundary conditions    
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
%                         caseName : String defining the case name
%                     pathToOutput : Path to where to write out the output 
%                                    file
%                            title : The title to the VTK file
%                       DOF4Output : Array containing the arrangment of the 
%                                    DOFs for printing them out
%                 writeOutputToVTK : Function handle to write out the 
%                                    results from each time step to a VTK 
%                                    file
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
%   7iv. Solve the mesh motion problem and update the mesh node locations and velocities
%
%    7v. Compute the load vector at the current time step
%
%   7vi. Save the discrete primary field and its first and second time derivatives
%
%  7vii. Solve the equation system
%
% 7viii. Update the time derivatives of the field
%
%   7ix. Write out the results into a VTK file or save them into an output variable
% <-
% 
%
%% Function main body
if ~isfield(propTransientAnalysis,'method')
    error('Undefined time integration method propStrDynamics.method');
end
if ~isfield(propTransientAnalysis,'T0')
    error('Undefined start time propStrDynamics.TStart');
end
if ~isfield(propTransientAnalysis,'TEnd')
    error('Undefined end time propStrDynamics.TEnd');
end
if ~isfield(propTransientAnalysis,'noTimeSteps')
    error('Undefined number of time steps propStrDynamics.noTimeSteps');
end
if ~isfield(propTransientAnalysis,'dt')
    error('Undefined time step size propStrDynamics.dt');
end
if isfield(propTransientAnalysis,'damping')
    if ~isfield(propTransientAnalysis.damping,'method')
        error('Undefined damping method propStrDynamics.damping.method');
    end
    if ~isfield(propTransientAnalysis.damping,'computeDampMtx')
        error('Undefined damping method propStrDynamics.damping.computeDampMtx');
    else
        if isempty(propTransientAnalysis.damping.computeDampMtx)
            error('String propStrDynamics.damping.computeDampMtx does not define a function');
        end
    end
end
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'Transient analysis properties\n'));
    fprintf(strcat(tab,'-----------------------------\n\n'));
    fprintf(strcat(tab,'Time integration method: %s \n'),propTransientAnalysis.method);
    fprintf(strcat(tab,'Start time = %d \n'),propTransientAnalysis.T0);
    fprintf(strcat(tab,'End time = %d \n'),propTransientAnalysis.TEnd);
    fprintf(strcat(tab,'No. time steps = %d \n'),propTransientAnalysis.noTimeSteps);
    fprintf(strcat(tab,'Time step size = %d \n'),propTransientAnalysis.dt);
    if isfield(propTransientAnalysis,'damping')
        fprintf(strcat(tab,'Damping method: %s \n'),propTransientAnalysis.damping.method);
    else
        fprintf(strcat(tab,'No damping is selected\n'));
    end
    fprintf('\n');
end

%% 0. Read input

% Compute the number of nodes in the fluid mesh
noNodes = length(msh.nodes(:,1));

% Compute the number of degrees of freedom
if strcmp(analysis.type,'NAVIER_STOKES_2D')
    noDOFs = 3*noNodes;
elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
    noDOFs = 4*noNodes;
elseif strcmp(analysis.type,'planeStress') || strcmp(analysis.type,'planeStrain')
    noDOFs = 2*noNodes;
else
    error('Select a valid analysis type in analysis.type');
end

% Initialize the scalar concetration history
FHistory = sparse(noDOFs,propTransientAnalysis.noTimeSteps + 1);

% Initialize the output array for the state variable
if isa(writeOutputToVTK,'function_handle')
    uHistory = 'undefined';
else
    uHistory = zeros(noDOFs,propTransientAnalysis.noTimeSteps + 1); 
end

%% 1. Get the initial values for the discrete solution vector and its first and second order rate
[u,uDot,uDDot,noTimeStep] = computeInitCnds...
    (analysis,msh,DOF4Output,parameters,propTransientAnalysis,...
    VTKResultFile,caseName,pathToOutput);

%% 2. Initialize the simulation time
t = propTransientAnalysis.dt*noTimeStep;

%% 3. Write out the initial conditions to a file
if isa(writeOutputToVTK,'function_handle')
    writeOutputToVTK(analysis,propNLinearAnalysis,propTransientAnalysis,...
        msh,parameters,u,uDot,uDDot,DOF4Output,caseName,pathToOutput,...
        title,noTimeStep);
else
    uHistory(:,1) = u;
end

%% 4. Compute the mass matrix of the problem
if isa(computeMassMtx,'function_handle')
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(tab,'Computing the mass matrix of the system\n'));
        fprintf(strcat(tab,'---------------------------------------\n\n'));
    end
    massMtx = computeMassMtx(analysis,msh,parameters,propGaussInt);
else
    error('Variable computeMassMtx is not defining a function handle as expected');
end

%% 5. Compute part of the residual vector and the stiffness matrix which stays constant throughout the transient analysis
if ~ischar(computeConstantProblemMatrices)
    KConstant = computeConstantProblemMatrices();
end

%% 6. Compute the damping matrix of the problem
if isfield(propTransientAnalysis,'damping')
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(tab,'Computing the damping matrix of the system\n'));
        fprintf(strcat(tab,'------------------------------------------\n\n'));
    end 
    dampMtx = propTransientAnalysis.damping.computeDampMtx...
        (analysis,msh,parameters,propGaussInt);
else
    dampMtx = 'undefined';
end

%% 7. Loop over all the time instances of the simulation
if strcmp(outMsg,'outputEnabled')
    fprintf('\tLooping over all the time steps \n\t------------------------------- \n\n');
end
while t < propTransientAnalysis.TEnd && noTimeStep < propTransientAnalysis.noTimeSteps
    %% 7i. Update the simulation time
    t = t + propTransientAnalysis.dt;
    
	%% 7ii. Update the time counter
    noTimeStep = noTimeStep + 1;
    
    %% 7iii. Preamble of the time stepping iterations
    if strcmp(outMsg,'outputEnabled')
        if ~ischar(propNLinearAnalysis)
            msgTS = sprintf(strcat(tab,'\tTime step %d/%d at real time %d seconds with dt=%d and maxNoNRIter=%d \n \n'),...
                noTimeStep,propTransientAnalysis.noTimeSteps,t,propTransientAnalysis.dt,propNLinearAnalysis.maxIter);
        else
            msgTS = sprintf(strcat(tab,'\tTime step %d/%d at real time %d seconds with dt=%d \n \n'),...
                noTimeStep,propTransientAnalysis.noTimeSteps,t,propTransientAnalysis.dt);
        end
        fprintf(msgTS);
    end
    
    %% 7iv. Solve the mesh motion problem and update the mesh node locations and velocities
    if ~ischar(nodesALE)
        [msh,uMeshALE,inhomDOFs,valuesInhomDOFs] = ...
            computeUpdatedMesh(msh,homDOFs,inhomDOFs,valuesInhomDOFs,nodesALE,...
            solve_LinearSystem,propTransientAnalysis,t);
    elseif strcmp(nodesALE,'undefined')
        uMeshALE = 'undefined';
    end
    
    %% 7v. Compute the load vector at the current time step
    if ~ischar(computeLoadVct)
        FHistory(:,noTimeStep) = computeLoadVct(msh,analysis,NBC,t,propGaussInt,'');
    elseif strcmp(computeLoadVct,'undefined')
        FHistory(:,noTimeStep) = zeros(noDOFs,1);  
    end
        
    %% 7vi. Save the discrete primary field and its first and second time derivatives
    uSaved = u;
    uDotSaved = uDot;
    uDDotSaved = uDDot;
    
    %% 7vii. Solve the equation system
    [u,~,hasConverged,minElSize] = solve_FEMEquationSystem...
        (analysis,uSaved,uDotSaved,uDDotSaved,msh,FHistory(:,noTimeStep),...
        computeBodyForceVct,parameters,u,uDot,uDDot,massMtx,dampMtx,...
        computeProblemMatricesSteadyState,DOFNumbering,freeDOFs,homDOFs,...
        inhomDOFs,valuesInhomDOFs,uMeshALE,solve_LinearSystem,...
        propTransientAnalysis,t,propNLinearAnalysis,propGaussInt,...
        strcat(tab,'\t'),outMsg);
    
    %% 7viii. Update the time derivatives of the field
    if hasConverged
        if isa(propTransientAnalysis.computeUpdatedVct,'function_handle')
            [uDot,uDDot] = propTransientAnalysis.computeUpdatedVct ...
                (u,uSaved,uDotSaved,uDDotSaved,propTransientAnalysis);
        else
            error('Function handle propTransientAnalysis.computeUpdatedVct undefined');
        end
    end
    
    %% 7ix. Write out the results into a VTK file or save them into an output variable
    if isa(writeOutputToVTK,'function_handle')
        writeOutputToVTK(analysis,propNLinearAnalysis,propTransientAnalysis,...
            msh,parameters,u,uDot,uDDot,DOF4Output,caseName,pathToOutput,...
            title,noTimeStep);
    else
        uHistory(:,noTimeStep + 1) = u;
    end
end

end
