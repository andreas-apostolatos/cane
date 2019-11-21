function [up,FComplete,hasConverged,minElSize] = solve_FEMVMSStabSteadyStateNSE2D...
    (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,nodesALE,parameters,...
    computeBodyForces,analysis,computeInitCnds,VTKResultFile,...
    solve_LinearSystem,propFldDynamics,propNLinearAnalysis,gaussInt,...
    caseName,outMsg)
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
%             homDOFs : The global numbering of the DOFs where homogeneous
%                       Dirichlet boundary conditions are applied
%           inhomDOFs : The global numbering of the DOFs where
%                       inhomogeneous Dirichlet boundary conditions are
%                       applied
%     valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                       Dirichlet boundary conditions are applied
%            nodesALE : The nodes on the ALE boundary:
%                           .nodes : The sequence of the nodal coordinates
%                                    on the ALE boundary
%                       .fcthandle : Function handle to the computation of
%                                    the ALE motion
%          parameters : Flow parameters
%   computeBodyForces : Function handle to the computation of the body
%                       force vector
%            analysis : .type : The analysis type
%     computeInitCnds : Function handle to the initial boundary conditions 
%                       computation
%       VTKResultFile : The name of the result file in the output folder
%                        where to get the initial conditions for the
%                        transient simulation
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
%            gaussInt : On the Gauss Point integration
%                             .type : 'default', 'user'
%                       .domainNoGP : Number of Gauss Points for the domain
%                                     integration
%                     .boundaryNoGP : Number of Gauss Points for the
%                                     boundary integration
%            caseName : String defining the case name
%              outMsg : On printing information during analysis in the
%                       command window
%
%              Output :
%                  up : The solution field in terms of the velocity and the 
%                       pressure field
%           FComplete : The complete force vector
%        hasConverged : Flag on the convegence of the nonlinear system
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

    % start measuring computational time
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
uMeshALE = 'undefined';
massMtx = 'undefined';
dampMtx = 'undefined';
t = 'undefined';
noTimeStep = 0;

% Define tabulation for the output to the command window
tab = '\t';

% Compute the number of nodes in the fluid mesh
noNodes = length(fldMsh.nodes(:,1));

% Compute the number of degrees of freedom
noDOFs = 3*noNodes;

% Assign a sequential numbering to the system DOFs
DOFNumbering = 1:noDOFs;

% Get the DOF numbering for each component of the displacement field and
% the pressure seperately
DOF4Output = [1:3:noDOFs-2; 2:3:noDOFs-1; 3:3:noDOFs];

% Prediction
up = zeros(noDOFs,1);

% Computation of the force vector
F = zeros(noDOFs,1);

% Title for the VTK files
title = 'Steady-state stabilized finite element formulation for the 2D incopmpressible Navier Stokes equations';

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs,inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs,prescribedDoFs)) = [];

%% 2. Solve the mesh motion problem and update the mesh node locations and velocities
if ~ischar(nodesALE)
    [fldMsh,uMeshALE,inhomDOFs,valuesInhomDOFs] = ...
        computeUpdatedMeshAndVelocitiesPseudoStrALE2D...
        (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,nodesALE,...
        solve_LinearSystem,propFldDynamics,t);
elseif strcmp(nodesALE,'undefined')
    uMeshALE = 'undefined';
end

%% 3. Solve the steady-state nonlinear Navier-Stokes stabilized finite element equation system
[up,FComplete,hasConverged,minElSize] = solve_FEMNLinearSystem...
    (analysis,uSaved,uDotSaved,uDDotSaved,fldMsh,F,computeBodyForces,...
    parameters,up,uDot,uDDot,massMtx,dampMtx,...
    @computeFEMVMSStabMtxAndVct4SteadyStateNLinear4NSE2D,DOFNumbering,...
    freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,uMeshALE,solve_LinearSystem,...
    propFldDynamics,t,propNLinearAnalysis,gaussInt,tab,outMsg);

%% 4. Write out the results into a VTK file
writeOutputFEMIncompressibleFlowToVTK(analysis,propNLinearAnalysis,...
    propFldDynamics,fldMsh,parameters,up,uDot,uDDot,DOF4Output,caseName,...
    pathToOutput,title,noTimeStep);

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Steady-state nonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('_______________Steady-State Nonlinear Analysis Ended_______________\n');
    fprintf('###################################################################\n\n\n');
end

end
