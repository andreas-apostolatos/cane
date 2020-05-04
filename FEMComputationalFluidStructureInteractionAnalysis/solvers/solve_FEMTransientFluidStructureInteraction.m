function [upHistory, dHistory, resVctHistoryFld, resVctHistoryStr, ...
    minElSizeFld, minElSizeStr] = ...
    solve_FEMTransientFluidStructureInteraction ...
    (propAnalysisFld, propAnalysisStr, fldMsh, strMsh, homDOFsFld, homDOFsStr, ...
    inhomDOFsFld, inhomDOFsStr, valuesInhomDOFsFld, valuesInhomDOFsStr, ...
    updateInhomDOFsFld, updateInhomDOFsStr, propALE, ...
    propParametersFld, propParametersStr, computeBodyForcesFld, computeBodyForcesStr, ...
    computeInitCndsFld, computeInitCndsStr, solve_FEMSystemFld, solve_FEMSystemStr, ...
    computeMassMtxFld, computeMassMtxStr, computeMtxSteadyStateFld, ...
    computeMtxSteadyStateStr, solve_LinearSystemFld, solve_LinearSystemStr, ...
    propFldDynamics, propStrDynamics, propNLinearAnalysisFld, ...
    propNLinearAnalysisStr, propFSI, propIDBCFld, propIDBCStr, ...
    propPostProcFld, propPostProcStr, propGaussIntFld, propGaussIntStr, ...
    propOutputFld, propOutputStr, caseNameFld, caseNameStr, ...
    outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the solution of the transient coupled fluid-structure interaction
% simulation using the Finite Element Method.
%
%                    Input :
%          propAnalysisFld : Structure containing general information on 
%                            the fluid analysis,
%                               .type : Analysis type
%          propAnalysisStr : Structure containing general information on 
%                            the structural analysis,
%                               .type : Analysis type
%           fldMsh, strMsh : Structures containing information on the fluid
%                            and structural meshes,
%                             .initialNodes : The initial nodes of the mesh
%                                    .nodes : The updated nodes of the mesh
%                                 .elements : The elements of the mesh
%   homDOFsFld, homDOFsStr : The global numbering of the DOFs where
%                            homogeneous Dirichlet boundary conditions are
%                            applied on the fluid and on the structure
%             inhomDOFsFld : The global numbering of the DOFs where
%                            inhomogeneous Dirichlet boundary conditions 
%                            are applied on the fluid side
%             inhomDOFsStr : The global numbering of the DOFs where
%                            inhomogeneous Dirichlet boundary conditions 
%                            are applied on the structural side
%       valuesInhomDOFsFld : The prescribed values of the DOFs where
%                            inhomogeneous Dirichlet boundary conditions 
%                            are applied on the fluid side
%       updateInhomDOFsStr : The prescribed values of the DOFs where
%                            inhomogeneous Dirichlet boundary conditions 
%                            are applied on the structural side
%               propALEFld : Structure containing information about the ALE
%                            motion of the fluid
%        propParametersFld : Structure containing information on the flow 
%                            parameters,
%                                .nue : Dynamic viscosity
%        propParametersStr : Structure containing information on the
%                            material properties of the structure,
%                                  .E : Young's modulus
%                                .nue : Poisson ratio
%                                  .t : Thickness
%     computeBodyForcesFld : Function handle to the compuation of the body
%                            forces for the fluid
%     computeBodyForcesStr : Function handle to the compuation of the body
%                            forces for the structure
%       computeInitCndsFld : Function handle to the compuation of the 
%                            initial conditions for the fluid problem
%       computeInitCndsStr : Function handle to the compuation of the 
%                            initial conditions for the structural problem
%       solve_FEMSystemFld : solve_FEMLinearSystem, solve_FEMNLinearSystem
%       solve_FEMSystemStr : solve_FEMLinearSystem, solve_FEMNLinearSystem
%        computeMassMtxFld : Function handle to the computation of the 
%                            fluid mass matrix
%        computeMassMtxStr : Function handle to the computation of the
%                            structural mass matrix
% computeMtxSteadyStateFld : Function handle to the computation of the
%                            steady-state matrices of the fluid problem
%    solve_LinearSystemFld : Function handle to the linear equation system
%                            solver for the fluid problem
%    solve_LinearSystemStr : Function handle to the linear equation system
%                            solver for the structural problem
%          propFldDynamics : Structure containing information on the 
%                            fluid dynamics,
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
%                                              .alphaBeta : Parameter for 
%                                                           the Bossak 
%                                                           scheme
%                                                  .gamma : Parameter for 
%                                                           the Bossak 
%                                                           scheme
%                                                 .TStart : Start time of 
%                                                           the simulation
%                                                   .TEnd : End time of the 
%                                                           simulation
%                                            .noTimeSteps : Number of time 
%                                                           steps
%                                                     .dt : Time step size
%          propStrDynamics : Structure containing information on the 
%                            structural dynamics,
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
%   propNLinearAnalysisFld : Structure containing information on the 
%                            nonlinear fluid analysis,
%                                .method : The nonlinear solution scheme
%                                   .eps : The residual tolerance
%                               .maxIter : The maximum number of nonlinear
%                                          iterations
%   propNLinearAnalysisStr : Structure containing information on the 
%                            nonlinear structural analysis,
%                                .method : The nonlinear solution scheme
%                                   .eps : The residual tolerance
%                               .maxIter : The maximum number of nonlinear
%                                          iterations
%                  propFSI : Structure containing information on the FSI
%                            coupling algorithm,
%                               .relaxation : Under-relaxation factor
%                                      .tol : Relative residual tolerance
%                                             on the displacement field
%              propIDBCFld : Structure containing information on the
%                            inhomogeneous Dirichlet boundary conditions 
%                            for the fluid problem
%              propIDBCStr : Structure containing information on the
%                            inhomogeneous Dirichlet boundary conditions 
%                            for the structural problem
%          propPostProcFld : Structure containing information on the
%                            computation of postprocessing resultants for
%                            the fluid problem,
%                               .nameDomain : Name of the domain onto which
%                                             postprocessing resultants are
%                                             to be computed
%                              .nodesDomain : IDs of the nodes which are on
%                                             the selected domain
%                          .computePostProc : Function handle to the
%                                             computation of the desirable 
%                                             resultant
%          propPostProcStr : Structure containing information on the
%                            computation of postprocessing resultants for
%                            the structural problem,
%                               .nameDomain : Name of the domain onto which
%                                             postprocessing resultants are
%                                             to be computed
%                              .nodesDomain : IDs of the nodes which are on
%                                             the selected domain
%                          .computePostProc : Function handle to the
%                                             computation of the desirable 
%                                             resultant
%          propGaussIntFld : Structure containing information on the
%                            numerical integration for the fluid problem,
%                              .type : 'default', 'user'
%                        .domainNoGP : Number of Gauss Points for the 
%                                      domain integration
%                      .boundaryNoGP : Number of Gauss Points for the
%                                      boundary integration
%          propGaussIntStr : Structure containing information on the
%                            numerical integration for the structural 
%                            problem,
%                              .type : 'default', 'user'
%                        .domainNoGP : Number of Gauss Points for the 
%                                      domain integration
%                      .boundaryNoGP : Number of Gauss Points for the
%                                      boundary integration
%            propOutputFld : Structure containing information on writting 
%                            the results of the fluid solution in file,
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
%            propOutputStr : Structure containing information on writting 
%                            the results of the fluid solution in file,
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
%             caseNameFld : Name of the case for the fluid problem
%             caseNameStr : Name of the case for the structural problem
%                  outMsg : On printing information during analysis in the
%                           command window
%
%                 Output :
%              upHistory : The solution history of the fluid problem
%               dHistory : The solution history of the structural problem
%       resVctFldHistory : The history of the residual vector for the fluid
%                          problem
%       resVctStrHistory : The history of the residual vector for the
%                          structural problem
%           minElSizeFld : The minimum element area size over the finite
%                         element fluid mesh
%           minElSizeStr : The minimum element area size over the finite
%                         element structural mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Find the prescribed and the free DOFs of the fluid system
%
% 2. Find the prescribed and the free DOFs of the fluid system
%
% 3. Solve the transient coupled system
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('______________________________________________________________________\n');
    fprintf('######################################################################\n');
    fprintf('Computational Fluid-Structure Interaction analysis has been initiated. \n\n');
    if ~isfield(propAnalysisFld, 'type')
        error('Structure propAnalysisFld must define field type');
    else
        if ~ischar(propAnalysisFld.type)
            error('propAnalysisFld.type must be a character array');
        end
    end
    fprintf('Fluid solver: %s \n\n', propAnalysisFld.type);
    if ~isfield(propAnalysisStr, 'type')
        error('Structure propAnalysisStr must define field type');
    else
        if ~ischar(propAnalysisStr.type)
            error('propAnalysisStr.type must be a character array');
        end
    end
    fprintf('Structural solver: %s \n\n', propAnalysisStr.type);
    fprintf('Coupling properties\n');
    fprintf('-------------------\n\n');
    if ~isfield(propFSI, 'relaxation')
        error('propFSI must contain field relaxation')
    end
    fprintf('Under-relaxation factor: %d\n', propFSI.relaxation);
    if ~isfield(propFSI, 'tol')
        error('propFSI must contain field tol')
    end
    fprintf('Coupling Iteration tolerance: %d\n', propFSI.tol);
    if ~isfield(propFSI, 'maxIter')
        error('propFSI must contain field maxIter')
    end
    fprintf('Maximum number of coupling iterations: %d\n\n', propFSI.maxIter);
    fprintf('______________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Check the structural and fluid dynamics properties for potential mismatch
if ~strcmp(propFldDynamics.timeDependence, propStrDynamics.timeDependence)
    error('Time dependencies for the fluid and structural dynamics do not match');
end
if propFldDynamics.T0 ~= propStrDynamics.T0
    error('Start times for the fluid and structural dynamics do not match');
end
if propFldDynamics.TEnd ~= propStrDynamics.TEnd
    error('End times for the fluid and structural dynamics do not match');
end
if propFldDynamics.noTimeSteps ~= propStrDynamics.noTimeSteps
    error('Number of time steps for the fluid and structural dynamics do not match');
end

% Find the dimensionality of the fluid problem
if strcmp(propAnalysisFld.type, 'NAVIER_STOKES_2D')
    numDOFsNodeFld = 3;
else
    error('Error in the dimensionality of the analysis');
end
if strcmp(propAnalysisStr.type, 'planeStress') || ...
        strcmp(propAnalysisStr.type, 'planeStrain')
    numDOFsNodeStr = 2;
elseif strcmp(propAnalysisStr.type, 'SDOF')
    numDOFsNodeStr = 1;
else
    error('Error in the dimensionality of the analysis');
end

% Output data to a VTK format
pathToOutputFld = '../../outputVTK/FEMFluidStructureInteraction/';
pathToOutputStr = '../../outputVTK/FEMFluidStructureInteraction/';

% Dummy variables
propNBCFld = 'undefined';
propNBCStr = 'undefined';
propALEStr = 'undefined';
computeUpdatedMeshStr = 'undefined';
computeLoadVctFld = 'undefined';
computeConstMtxFld = 'undefined';
computeConstMtxStr = 'undefined';

% Function handle to the load vector computation for the structural problem
if strcmp(propAnalysisStr.type, 'planeStress') || ...
        strcmp(propAnalysisStr.type, 'planeStrain')
    computeLoadVctStr = @computeLoadVctFEMPlateInMembraneAction;
else
    computeLoadVctStr = 'undefined';
end

% Function handle to the computation of the updated mesh
computeUpdatedMeshFld = @computeUpdatedMeshAndVelocitiesPseudoStrALE2D;

% Compute the number of nodes in both meshes
numNodesFld = length(fldMsh.nodes(:, 1));
numNodesStr = length(strMsh.nodes(:, 1));

% Compute the number of degrees of freedom
numDOFsFld = numDOFsNodeFld*numNodesFld;
numDOFsStr = numDOFsNodeStr*numNodesStr;

% Assign a sequential numbering to the system DOFs
DOFNumberingFld = 1:numDOFsFld;
DOFNumberingStr = 1:numDOFsStr;

% Get the DOF numbering for each component of the velocity field and the 
% pressure seperately
DOF4OutputFld = [1:3:numDOFsFld - 2
                 2:3:numDOFsFld - 1
                 3:3:numDOFsFld];

% Get the DOF numbering for the displacement field
DOF4OutputStr = [1:2:numDOFsStr - 1
                 2:2:numDOFsStr];

% Title for the VTK files
titleFld = 'Stabilized finite element formulation for the 2D incopmpressible Navier Stokes equations';
titleStr = 'Geometrically nonlinear transient plane stress analysis';

% Define tabulation for the output in the command window
tab = '\t';

%% 1. Find the prescribed and the free DOFs of the fluid system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDOFsFld = mergesorted(homDOFsFld, inhomDOFsFld);
prescribedDOFsFld = unique(prescribedDOFsFld);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFsFld = DOFNumberingFld;
freeDOFsFld(ismember(freeDOFsFld, prescribedDOFsFld)) = [];

%% 2. Find the prescribed and the free DOFs of the fluid system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDOFsStr = mergesorted(homDOFsStr, inhomDOFsStr);
prescribedDOFsStr = unique(prescribedDOFsStr);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFsStr = DOFNumberingStr;
freeDOFsStr(ismember(freeDOFsStr, prescribedDOFsStr)) = [];

%% 3. Solve the transient coupled system
[upHistory, dHistory, resVctHistoryFld, resVctHistoryStr, ...
    minElSizeFld, minElSizeStr] = ...
    solve_FEMCoupledTransientAnalysis ...
    (propAnalysisFld, propAnalysisStr, fldMsh, strMsh, DOFNumberingFld, DOFNumberingStr,...
    freeDOFsFld, freeDOFsStr, homDOFsFld, homDOFsStr, inhomDOFsFld, inhomDOFsStr, ...
    valuesInhomDOFsFld, valuesInhomDOFsStr, updateInhomDOFsFld, updateInhomDOFsStr, ...
    propALE, propALEStr, computeInitCndsFld, computeInitCndsStr, ...
    computeBodyForcesFld, computeBodyForcesStr, propNBCFld, propNBCStr, ...
    computeLoadVctFld, computeLoadVctStr, propParametersFld, propParametersStr, ...
    @solve_FEMFluidStructureInteractionPartitionedGaussSeidel, ...
    solve_FEMSystemFld, solve_FEMSystemStr, computeConstMtxFld, computeConstMtxStr, ...
    computeMassMtxFld, computeMassMtxStr, computeMtxSteadyStateFld, computeMtxSteadyStateStr, ...
    computeUpdatedMeshFld, computeUpdatedMeshStr, solve_LinearSystemFld, ...
    solve_LinearSystemStr, propFldDynamics, propStrDynamics, ...
    propNLinearAnalysisFld, propNLinearAnalysisStr, propFSI, propIDBCFld, propIDBCStr, ...
    propPostProcFld, propPostProcStr, propGaussIntFld, propGaussIntStr, propOutputFld, ...
    propOutputStr, caseNameFld, caseNameStr, pathToOutputFld, pathToOutputStr, titleFld, ...
    titleStr, DOF4OutputFld, DOF4OutputStr, tab, outMsg);

%% 4. Appendix
if strcmp(outMsg, 'outputEnabled')
    computationalTime = toc;
    fprintf('Fluid-Structure Interaction analysis took %.2d seconds \n\n',computationalTime);
    fprintf('______________Fluid-Structure Interaction Analysis Ended______________\n');
    fprintf('######################################################################\n\n\n');
end

end