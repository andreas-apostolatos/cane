function [upHistory,minElSize] = solve_FEMVMSStabTransientNSEBossakTI2D...
    (fldMsh,homDOFs,inhomDOFs,valuesInhomDOFs,nodesALE,parameters,...
    computeBodyForces,analysis,computeInitCnds,VTKResultFile,...
    solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
    writeOutputToFile,gaussInt,caseName,outMsg)
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
% Bossak time integration scheme for the temporal discretization and the 
% classical triangular basis functions for the spatial discretization.
%
% Unconditional stability of the Bossak time integration scheme is ensured 
% if the following relations hold :
% 
% - alphaBeta <= .5 
% - beta >= gamma/2 >= .25
% - alphaBeta + gamma >= .25
%
%                Input :
%               fldMsh : Nodes and elements for the fluid mesh
%              homDOFs : The global numbering of the DOFs where homogeneous
%                        Dirichlet boundary conditions are applied
%            inhomDOFs : The global numbering of the DOFs where
%                        inhomogeneous Dirichlet boundary conditions are
%                        applied
%      valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                        Dirichlet boundary conditions are applied
%             nodesALE : The nodes on the ALE boundary:
%                            .nodes : The sequence of the nodal coordinates
%                                     on the ALE boundary
%                        .fcthandle : Function handle to the computation of
%                                     the ALE motion
%           parameters : Flow parameters
%    computeBodyForces : Function handle to the computation of the body
%                        force vector
%             analysis : .type : The analysis type
%      computeInitCnds : Function handle to the initial boundary conditions 
%                        computation
%        VTKResultFile : The name of the result file in the output folder
%                         where to get the initial conditions for the
%                         transient simulation
%   solve_LinearSystem : Function handle to the solver for the linear 
%                        equation system
%      propFldDynamics : On the transient analysis :
%                              .method : The time integration method
%                           .alphaBeta : (parameter for the Bossak scheme)
%                               .gamma : (parameter for the Bossak scheme)
%                              .TStart : Start time of the simulation
%                                .TEnd : End time of the simulation
%                                  .nT : Number of time steps
%                                  .dt : Time step
%  propNLinearAnalysis : On the nonlinear analysis :
%                                .method : The nonlinear solution scheme
%                                   .eps : The residual tolerance
%                               .maxIter : The maximum number of nonlinear
%                                          iterations
%    writeOutputToFile : Function handle to writting the output into a file
%             gaussInt : On the Gauss Point integration
%                              .type : 'default', 'user'
%                        .domainNoGP : Number of Gauss Points for the domain
%                                      integration
%                      .boundaryNoGP : Number of Gauss Points for the
%                                      boundary integration
%             caseName : String defining the case name
%               outMsg : On printing information during analysis in the
%                        command window
%
%               Output :
%            upHistory : The history of the velocity and pressure field
%                        throughout the transient analysis
%            minElSize : The minimum element area size over the isogeometric
%                        mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Find the prescribed and the free DOFs of the system
%
% 2. Solve the transient nonlinear Navier-Stokes stabilized finite element equation system using the Bossak time integration scheme
%
% 3. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________\n');
    fprintf('###################################################################\n');
    fprintf('Nonlinear transient analysis for the FEM-based 2D incompressible\n');
    fprintf('Navier-Stokes flow problem using the Bossak time integration scheme \n');
    fprintf('has been initiated. \n\n');
    if strcmp(propNLinearAnalysis.method,'Newton')
        fprintf('Nonlinear scheme : Newton method \n');
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

% Find the dimensionality of the problem
if strcmp(analysis.type,'NAVIER_STOKES_3D')
    isAnalysis3D = true;
    noDOFsNode = 4;
elseif strcmp(analysis.type,'NAVIER_STOKES_2D')
    isAnalysis3D = false;
    noDOFsNode = 3;
else
    error('Error in the dimensionality of the analysis');
end
% Output data to a VTK format
pathToOutput = '../outputVTK/FEMComputationalFluidDynamicsAnalysis/';

% Dummy variables
computeConstantMatrices = 'undefined';
NBC = 'undefined';
computeLoadVct = 'undefined';

% Define tabulation for the output in the command window
tab = '\t';

% Compute the number of nodes in the fluid mesh
noNodes = length(fldMsh.nodes(:,1));

% Compute the number of degrees of freedom
noDOFs = noDOFsNode*noNodes;

% Assign a sequential numbering to the system DOFs
DOFNumbering = 1:noDOFs;

% Get the DOF numbering for each component of the displacement field and
% the pressure seperately
if isAnalysis3D
    DOF4Output = [1:4:noDOFs-3; 2:4:noDOFs-2; 3:4:noDOFs-1; 4:4:noDOFs];
else
    DOF4Output = [1:3:noDOFs-2; 2:3:noDOFs-1; 3:3:noDOFs];
end

% Title for the VTK files
title = 'Stabilized finite element formulation for the 2D incopmpressible Navier Stokes equations';

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs,inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs,prescribedDoFs)) = [];

%% 2. Solve the transient nonlinear Navier-Stokes stabilized finite element equation system using the Bossak time integration scheme
[upHistory,minElSize] = solve_FEMTransientAnalysis...
    (analysis,fldMsh,DOFNumbering,freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,...
    nodesALE,computeInitCnds,VTKResultFile,computeBodyForces,NBC,...
    computeLoadVct,parameters,@solve_FEMNLinearSystem,...
    computeConstantMatrices,@computeMassMtx4FEMVMSStabNSE2D,...
    @computeFEMVMSStabMtxAndVct4NLinear4NSE,...
    @computeUpdatedMeshAndVelocitiesPseudoStrALE2D,solve_LinearSystem,...
    propFldDynamics,propNLinearAnalysis,gaussInt,caseName,pathToOutput,...
    title,DOF4Output,writeOutputToFile,tab,outMsg);

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Transient nonlinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('________________Transient Nonlinear Analysis Ended_________________\n');
    fprintf('###################################################################\n\n\n');
end

end
