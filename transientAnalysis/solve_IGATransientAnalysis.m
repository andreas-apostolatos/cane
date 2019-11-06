function [uHistory,resHistory,BSplinePatches,propCoupling,minElAreaSize] = ...
    solve_IGATransientAnalysis...
    (analysis,BSplinePatches,connections,freeDOFs,homDOFs,inhomDOFs,...
    valuesInhomDOFs,masterDOFs,slaveDOFs,nodesALE,computeInitCnds,...
    solve_IGAEquationSystem,computeConstantMatrices,computeMassMtx,...
    computeProblemMatricesSteadyState,computeUpdatedMesh,...
    computeUpdatedGeometry,solve_LinearSystem,propCoupling,...
    propTransientAnalysis,propNLinearAnalysis,propPostproc,caseName,...
    pathToOutput,title,propOutput,isReferenceUpdated,...
    propEmpireCoSimulation,tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the history of the primary variable corresponding to a transient 
% analysis over a field which is discrezed with isogeometric analysis.
%
%                          Input :
%                       analysis : Information on the analysis type
%                 BSplinePatches : Cell array of B-Spline patches each of 
%                                  which containing:
%                                   .p,q : The polynomial orders of the 
%                                          B-Spline surface in both 
%                                          parametric directions
%                                .Xi,Eta : The knot vectors in both 
%                                          parametric directions
%                                    .CP : The set of control points and 
%                                          weights
%                               .isNURBS : Flag on whether the basis is a 
%                                          NURBS or a B-Spline
%                            .parameters : Technical parameters for the 
%                                          structure
%                               .homDOFs : The global numbering of the DOFs 
%                                          where homogeneous Dirichlet 
%                                          boundary conditions are applied
%                             .inhomDOFs : The global numbering of the DOFs 
%                                          where inhomogeneous Dirichlet 
%                                          boundary conditions are applied
%                       .valuesInhomDOFs : The values on the DOFs 
%                                          corresponding to the application 
%                                          of inhomogeneous Dirichlet 
%                                          boundary conditions
%                            masterDOFs : Master DOFs of the system
%                             slaveDOFs : Slave DOFs of the system
%                                  .NBC : Structure containing information 
%                                         on the application of the Neumann 
%                                         boundary conditions
%                                           .noCnd : Number of Neumann 
%                                                    boundary conditions
%                                 .xiLoadExtension : Cell array {.noCnd} 
%                                                    containing the load 
%                                                    extensions in the xi-
%                                                    direction
%                                .etaLoadExtension : Cell array {.noCnd} 
%                                                    containing the load 
%                                                    extensions in the eta-
%                                                    direction
%                                   .loadAmplitude : Array (1,.noCnd) 
%                                                    containing the load 
%                                                    amplitudes
%                                   .loadDirection : Array (1,.noCnd) 
%                                                    containing the load 
%                                                    directions
%                                  .computeLoadVct : Cell array {.noCnd} 
%                                                    containing the 
%                                                    function name for the 
%                                                    computation of the 
%                                                    load vector
%                                  .isConservative : Array (1,.noCnd) of 
%                                                    flags indicating 
%                                                    whether the load is 
%                                                    conservative or not
%                    connections : Array defining the connection of the
%                                  B-Spline patches as well as the
%                                  extension of the coupling boundaries
%                       freeDOFs : The global numbering of the DOFs which 
%                                  are uncontrained
%                        homDOFs : The global numbering of the DOFs over
%                                  which homogeneous Dirichlet boundary
%                                  conditions are applied
%                      inhomDOFs : The global numbering of the DOFs over
%                                  which inhomogeneous Dirichlet boundary
%                                  conditions are applied
%                valuesInhomDOFs : Prescribed values on the nodes where 
%                                  inhomogeneous boundary conditions are 
%                                  applied
%                       nodesALE : The nodes on the ALE boundary:
%                                       .nodes : The sequence of the nodal 
%                                                coordinates on the ALE 
%                                                boundary
%                                   .fcthandle : Function handle to the 
%                                                computation of the ALE 
%                                                motion
%                computeInitCnds : Function handle to the initial boundary 
%                                  conditions computation
%                            NBC : On the Neumann boundary conditions    
%                                      .nodes : The nodes where Neumann 
%                                               boundary conditions are 
%                                               applied
%                                   .loadType : The type of the load for 
%                                               each Neumann node
%                                  .fctHandle : The function handle for 
%                                               each Neumann node for the 
%                                               computation of the load 
%                                               vector (these functions are 
%                                               under the folder load)
%        solve_IGAEquationSystem : Function handle to the solution of the
%                                  equation system resulting out of the
%                                  isogeometric discretization
%        computeConstantMatrices : Function handle which returns the parts 
%                                  of the problem matrices which are linear 
%                                  with respect to the solution field
%                 computeMassMtx : Function handle to the computation of
%                                  the mass matrix of the system
% computeProblemMtrcsSteadyState : Function handle to the computation of 
%                                  the problem matrices corresponding to
%                                  the steady-state problem
%             computeUpdatedMesh : Function handle to the computation of
%                                  the updated mesh location and velocities 
%                                  at the nodes
%         computeUpdatedGeometry : Function handle to the computation of
%                                  the updated geometry
%                   propCoupling : Coupling properties for a multipatch
%                                  geometry
%             solve_LinearSystem : Function handle to the linear equation
%                                  system solver
%          propTransientAnalysis : On the transient analysis :
%                                        .method : The time integration method
%                                     .alphaBeta : (parameter for the 
%                                                   Bossak method)
%                                         .gamma : (parameter for the 
%                                                   Bossak method)
%                                        .TStart : Start time of the 
%                                                  simulation
%                                          .TEnd : End time of the 
%                                                  simulation
%                                   .noTimeSteps : Number of time steps
%                                            .dt : Time step
%                                       .damping : On the damping
%                                                       .method : Damping
%                                                                 method
%                  .computeProblemMtrcsTransient : Function handle to the
%                                                  computation of the
%                                                  problem matrices
%                                                  according to the defined
%                                                  time integration method
%                                                  given the system matrix 
%                                                  of the steady-state
%                                                  system as well as the
%                                                  mass matrix of the
%                                                  problem
%                             .computeUpdatedVct : Function handle to the 
%                                                  computation of the 
%                                                  updated values of the 
%                                                  discretized 
%                                  field and its time derivatives nessecary 
%                                  for marching in time
%            propNLinearAnalysis : On the nonlinear analysis :
%                                 .method : The nonlinear solution method
%                                    .eps : The residual tolerance
%                                .maxIter : The maximum number of nonlinear
%                                           iterations
%                   propPostproc :            .resultant : Array of strings 
%                                                          containing the 
%                                                          name of each 
%                                                          resultant to be 
%                                                          written out
%                                      .computeResultant : Array of strings 
%                                                          representing the 
%                                                          function handle 
%                                                          for the 
%                                                          computation of 
%                                                          the desirable 
%                                                          resultant
%                       caseName : String defining the case name
%                   pathToOutput : Path to where to write out the output 
%                                  file
%                          title : The title to the GiD file
%                     propOutput : Properties for writting out the results
%                              .writeOutput : Function handle to outputting 
%                                             the results in a folder
%                           .writeFrequency : Frequency for writting out
%                                             the results
%                            .saveFrequency : Frequency for saving temorary
%                                             results at the time step
%                                             level
%             isReferenceUpdated : Flag on whether the reference geometry 
%                                  is updated
%         propEmpireCoSimulation : Properties for the co-simulation with 
%                                  EMPIRE
%                            .isCoSimulation : Flag on whether 
%                                              co-simulation with EMPIRE is 
%                                              assumed
%                          .isInterfaceLayer : Flag on whether the matlab 
%                                              client is used as an 
%                                              interface layer
%             .isGaussPointProvidedDirichlet : Flag on whether the Gauss
%                                              Point data along the
%                                              Dirichlet boundaries are
%                                              provided
%             .isGaussPointProvidedInterface : Flag on whether the Gauss
%                                              Point data along the
%                                              interface boundaries are
%                                              provided
%                          .propIntDirichlet : Integration along the
%                                              Dirichlet boundaries
%                                              properties,
%                                                .type : 'default', 'user'
%                                               .noGPs : Number of Gauss
%                                                        Points
%                          .propIntInterface : Integration along the
%                                              interface boundaries
%                                              properties,
%                                                .type : 'default', 'user'
%                                               .noGPs : Number of Gauss
%                                                        Points
%                              .strMatlabXml : Name of the xml file for 
%                                              connection to Empire
%                            tab : Tabulation for the output messages
%                         outMsg : On printing information during analysis 
%                                  in the command window
%
%                         Output :
%                       uHistory : The history of the primary field
%                                  throughout the transient analysis
%                     resHistory : Structure array containing the
%                                  convergence history for each time step
%                 BSplinePatches : The updated array of B-Spline patches
%                   propCoupling : The updated structure of the coupling 
%                                  properties
%                      minElSize : The minimum element edge size in the
%                                  mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Initialize connection to EMPIRE
%
% 2. Send the multipatch isogeometric surface to EMPIRE
%
% 3. Get and write the initial values for the discrete solution vector and its first and second order rate, write out the initial conditions and initialize output array
%
% 4. Initialize the simulation time
%
% 5. Compute matrix which stays constant throughout the transient simulation
%
% 6. Compute the mass matrix of the problem
%
% 7. Compute the damping matrix of the problem
%
% 8. Perform a static analysis for the structure to come into equilibrium with its internal forces before starting the time loop
%
% 9. Update the displaced Control Point coordinates in case restart option is enabled
%
% 10. Loop over all the time instances of the simulation
% ->
%     10i. Update the simulation time
%
%    10ii. Preamble of the time stepping iterations
%
%   10iii. Solve the mesh motion problem and update the mesh node locations and velocities
%
%    10iv. Save the discrete primary field and its first and second time derivatives
%
%     10v. Solve the equation system
%
%    10vi. Update the time derivatives of the field
%
%   10vii. Save and write out the results into a GiD file
%
%  10viii. Save temporary data at the current time step
% <-
%
% 11. Erase temporary data
%
% 12. Finalize connection to EMPIRE
%
% 13. Check if all output variables are assigned a value
%
%% Function main body
if ~isfield(propTransientAnalysis,'method')
    error('Undefined time integration method propStrDynamics.method');
end
if ~isfield(propTransientAnalysis,'TStart')
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
    fprintf(strcat(tab,'Start time = %d \n'),propTransientAnalysis.TStart);
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

% Dummy arrays for this function
plot_IGANLinear = 'undefined';
graph = 'undefinded';

% Initialize restart flag
isRestart = false;

% Empty arrays for this function
masterDOFs = [];
slaveDOFs = [];

% Get number of DOFs of the system
noDOFs = length(freeDOFs) + length(homDOFs) + length(inhomDOFs);

% Get the number of weak Dirichlet boundary conditions
noWeakDBCCnd = 0;
for iPatches = 1:length(BSplinePatches)
    if isfield(BSplinePatches{iPatches},'weakDBC')
        if isfield(BSplinePatches{iPatches}.weakDBC,'noCnd')
            noWeakDBCCnd = noWeakDBCCnd + ...
                BSplinePatches{iPatches}.weakDBC.noCnd;
        end
    end
end

% Initialize a counters
counterRes = 1;
counterPrim = 1;

% Initialize output arrays
uHistory = zeros(noDOFs,propTransientAnalysis.noTimeSteps);
if ~ischar(propNLinearAnalysis)
    resHistory = zeros(propNLinearAnalysis.maxIter,propTransientAnalysis.noTimeSteps);
else
    resHistory = zeros(length('undefined'),propTransientAnalysis.noTimeSteps);
end

%% 1. Initialize connection to EMPIRE
if propEmpireCoSimulation.isCoSimulation
    initializeConnectionWithEmpire(propEmpireCoSimulation.strMatlabXml,tab,outMsg);
end

%% 2. Send the multipatch isogeometric surface to EMPIRE
if propEmpireCoSimulation.isCoSimulation
    sendBSplineGeometryToEmpire(BSplinePatches,connections,propEmpireCoSimulation,tab,outMsg);
    if propEmpireCoSimulation.isInterfaceLayer
        if strcmp(outMsg,'outputEnabled')
            fprintf(strcat(tab,'Interface layer option chosen, no analysis is performed\n'));
        end
        EMPIRE_API_Disconnect();
        return;
    end
end

%% 3. Get and write the initial values for the discrete solution vector and its first and second order rate, write out the initial conditions and initialize output array
[u,uDot,uDDot,noTimeStep,~] = computeInitCnds...
    (BSplinePatches,noDOFs,propTransientAnalysis,caseName,pathToOutput,tab,outMsg);
if noTimeStep ~= 0
    isRestart = true;
    uHistory(:,1:length(u(1,:))) = u;
    counterPrim = length(u(1,:)) + 1;
    u = u(:,end);
else
    uHistory(:,counterPrim) = u;
    counterPrim = counterPrim + 1;
    if ~ischar(propOutput)
        if isfield(propOutput,'writeOutput')
            if isa(propOutput.writeOutput,'function_handle')
                propOutput.writeOutput(analysis,BSplinePatches,u,...
                    uDot,uDDot,propNLinearAnalysis,propTransientAnalysis,...
                    propPostproc,caseName,pathToOutput,title,noTimeStep);
            end
        end
    end
end

%% 4. Initialize the simulation time
t = propTransientAnalysis.TStart + propTransientAnalysis.dt*noTimeStep;

%% 5. Compute matrix which stays constant throughout the transient simulation
if isa(computeConstantMatrices,'function_handle')
    constMtx = computeConstantMatrices(BSplinePatches,connections,noDOFs,propCoupling);
else
    constMtx = 'undefined';
end

%% 6. Compute the mass matrix of the problem
if isa(computeMassMtx,'function_handle')
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(tab,'Computing the mass matrix of the system\n'));
        fprintf(strcat(tab,'---------------------------------------\n\n'));
    end
    massMtx = computeMassMtx(BSplinePatches,noDOFs);
else
    error('Variable computeMassMtx is not defining a function handle as expected');
end

%% 7. Compute the damping matrix of the problem
if isfield(propTransientAnalysis,'damping')
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(tab,'Computing the damping matrix of the system\n'));
        fprintf(strcat(tab,'------------------------------------------\n\n'));
    end 
    dampMtx = propTransientAnalysis.damping.computeDampMtx...
        (constMtx,computeProblemMatricesSteadyState,...
        massMtx,noDOFs,BSplinePatches,connections,...
        propCoupling,propTransientAnalysis.TStart,propTransientAnalysis,...
        noWeakDBCCnd,isReferenceUpdated,tab,outMsg);
else
    dampMtx = 'undefined';
end

%% 8. Perform a static analysis for the structure to come into equilibrium with its internal forces before starting the time loop
propTransientAnalysis.isStaticStep = true;
if ~isRestart    
    if strcmp(outMsg,'outputEnabled')
        fprintf(strcat(tab,'Performing a static step to bring the structure in equilibrium with its internal forces \n'));
        fprintf(strcat(tab,'--------------------------------------------------------------------------------------- \n\n'));
    end
    [u,~,resHistory(:,counterRes),hasConverged,FComplete,rankD,condK,...
        minEig,BSplinePatches,propCoupling,minElAreaSize] = ...
        solve_IGAEquationSystem...
        (analysis,u,uDot,uDDot,BSplinePatches,connections,...
        u,uDot,uDDot,constMtx,massMtx,dampMtx,computeProblemMatricesSteadyState,...
        computeUpdatedGeometry,freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,...
        masterDOFs,slaveDOFs,solve_LinearSystem,t,propCoupling,...
        propTransientAnalysis,propNLinearAnalysis,plot_IGANLinear,...
        isReferenceUpdated,propEmpireCoSimulation.isCoSimulation,strcat(tab),graph,outMsg);
    counterRes = counterRes + 1;
    if ~hasConverged
        warning(strcat(tab,'\t','Nonlinear analysis did not converge'));
    end
    uHistory(:,counterPrim) = u;
    counterPrim = counterPrim + 1;
    if ~ischar(propOutput)
        if isfield(propOutput,'writeOutput')
            if isa(propOutput.writeOutput,'function_handle') && noTimeStep == 0  
                propOutput.writeOutput(analysis,BSplinePatches,u,...
                    uDot,uDDot,propNLinearAnalysis,propTransientAnalysis,...
                    propPostproc,caseName,pathToOutput,title,noTimeStep + 1);
            end
        end
    end
    if isfield(propOutput,'saveFrequency')
        save(['dataTemporary_' caseName]);
    end
end
propTransientAnalysis.isStaticStep = false;

%% 9. Update the displaced Control Point coordinates in case restart option is enabled
if isRestart
    BSplinePatches = computeUpdatedGeometry(BSplinePatches,u);
end

%% 10. Loop over all the time instances of the simulation
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'Looping over all the time steps \n'));
    fprintf(strcat(tab,'------------------------------- \n\n'));
end
for iTimeSteps = noTimeStep + 1:propTransientAnalysis.noTimeSteps
    %% 10i. Update the simulation time
    t = t + propTransientAnalysis.dt;
    
    %% 10ii. Preamble of the time stepping iterations
    if strcmp(outMsg,'outputEnabled')
        if ~ischar(propNLinearAnalysis)
            msgTS = sprintf(strcat(tab,'Time step %d/%d at real time %d seconds with dt=%d and maxNoNRIter=%d \n \n'),...
                iTimeSteps,propTransientAnalysis.noTimeSteps,t,propTransientAnalysis.dt,propNLinearAnalysis.maxIter);
        else
            msgTS = sprintf(strcat(tab,'Time step %d/%d at real time %d seconds with dt=%d\n\n'),...
                iTimeSteps,propTransientAnalysis.noTimeSteps,t,propTransientAnalysis.dt);
        end
        fprintf(msgTS);
    end
    
    %% 10iii. Solve the mesh motion problem and update the mesh node locations and velocities
    if ~ischar(nodesALE)
        [BSplinePatches,uMeshALE,inhomDOFs,valuesInhomDOFs] = ...
            computeUpdatedMesh(BSplinePatches,homDOFs,inhomDOFs,...
            valuesInhomDOFs,nodesALE,propTransientAnalysis,t);
    elseif strcmp(nodesALE,'undefined')
        uMeshALE = 'undefined';
    end
    
    %% 10iv. Save the discrete primary field and its first and second time derivatives
    uSaved = u;
    uDotSaved = uDot;
    uDDotSaved = uDDot;
    
    %% 10v. Solve the equation system
    [u,~,resHistory(:,counterRes),hasConverged,FComplete,rankD,...
        condK,minEig,BSplinePatches,propCoupling,minElAreaSize] = ...
        solve_IGAEquationSystem...
        (analysis,uSaved,uDotSaved,uDDotSaved,BSplinePatches,connections,...
        u,uDot,uDDot,constMtx,massMtx,dampMtx,computeProblemMatricesSteadyState,...
        computeUpdatedGeometry,freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,...
        masterDOFs,slaveDOFs,solve_LinearSystem,t,propCoupling,...
        propTransientAnalysis,propNLinearAnalysis,plot_IGANLinear,...
        isReferenceUpdated,propEmpireCoSimulation.isCoSimulation,strcat(tab,'\t'),...
        graph,outMsg);
    counterRes = counterRes + 1;
    
    %% 10vi. Update the time derivatives of the field
    if hasConverged
        if isa(propTransientAnalysis.computeUpdatedVct,'function_handle')
            [uDot,uDDot] = propTransientAnalysis.computeUpdatedVct ...
                (u,uSaved,uDotSaved,uDDotSaved,propTransientAnalysis);
        else
            error('Function handle propTransientAnalysis.computeUpdatedVct undefined');
        end
    end
    
    %% 10vii. Save and write out the results into a GiD file
    uHistory(:,counterPrim) = u;
    counterPrim = counterPrim + 1;
    if ~ischar(propOutput)
        if isfield(propOutput,'writeOutput')
            if isa(propOutput.writeOutput,'function_handle')
                if mod((iTimeSteps - noTimeStep - 1),propOutput.writeFrequency) == 0
                    propOutput.writeOutput(analysis,BSplinePatches,u,...
                        uDot,uDDot,propNLinearAnalysis,propTransientAnalysis,...
                        propPostproc,caseName,pathToOutput,title,iTimeSteps + 1);
                end
            end
        end
    end
    
    %% 10viii. Save temporary data at the current time step
    if isfield(propOutput,'saveFrequency')
        if mod((iTimeSteps - noTimeStep - 1),propOutput.saveFrequency) == 0
            save(['dataTemporary_' caseName]);
        end
    end
end

%% 11. Erase temporary data
if isfield(propOutput,'saveFrequency')
    delete(['dataTemporary_' caseName '.mat']);
end

%% 12. Finalize connection to EMPIRE
if propEmpireCoSimulation.isCoSimulation
    EMPIRE_API_Disconnect();
end

%% 13. Check if all output variables are assigned a value
if ~exist('minElAreaSize','var')
    minElAreaSize = 'undefined';
end

end
