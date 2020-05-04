function [up, d, upDot, dDot, upDDot, dDDot, resVctFld, resVctStr, ...
    fldMsh, strMsh, FFSIFld, FFSIStr, minElSizeFld, minElSizeStr] = ...
    solve_FEMFluidStructureInteractionPartitionedGaussSeidel ...
    (propAnalysisFld, propAnalysisStr, upSaved, dSaved, ...
    upDotSaved, dDotSaved, upDDotSaved, dDDotSaved, ...
    fldMsh, strMsh, FFld, FStr, FFSIFld, FFSIStr, computeBodyForceVctFld, ...
    computeBodyForceVctStr, propParametersFld, propParametersStr, up, d, ...
    upDot, dDot, upDDot, dDDot, massMtxFld, massMtxStr, dampMtxFld, ...
    dampMtxStr, precompStiffMtxFld, precompStiffMtxStr, precomResVctFld, ...
    precomResVctStr, computeMtxSteadyStateFld, computeMtxSteadyStateStr, ...
    solve_FEMEquationSystemFld, solve_FEMEquationSystemStr, ...
    DOFNumberingFld, DOFNumberingStr, freeDOFsFld, freeDOFsStr, homDOFsFld, ...
    homDOFsStr, inhomDOFsFld, inhomDOFsStr, valuesInhomDOFsFld, valuesInhomDOFsStr, ...
    propALE, propALEStr, computeUpdatedMeshFld, computeUpdatedMeshStr, ...
    solve_LinearSystemFld, solve_LinearSystemStr, propFldDynamics, propStrDynamics, ...
    t, propNLinearAnalysisFld, propNLinearAnalysisStr, propFSI, ...
    propPostProcFld, propPostProcStr, propGaussIntFld, propGaussIntStr, ...
    tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the converged solution for the coupled fluid-structure
% interaction iterations using the partitioned Gauss-Seidel approach.
%
%                     Input :
%             fldMsh,strMsh : Structures containing the fluid and the
%                             structural meshes,
%                           .initialNodes : Initial nodes in the finite
%                                           element mesh
%                                  .nodes : Current nodes in the finite
%                                           element mesh
%                               .elements : Elements in the finite element
%                                           mesh
%                 FFld,FStr : Externally applied load vectors for the fluid
%                             and the structural problems
%           FFSIFld,FFSIStr : Interface forces vectors from the previous
%                             Gauss-Seidel coupling iterations
%    computeBodyForceVctFld : Function handle to the computation of the 
%                             body forces for the structural problem
%    computeBodyForceVctStr : Function handle to the computation of the 
%                             body forces for the fluid problem
%         propParametersFld : Structure containing the flow parameters,
%                                  .nue : Dynamic viscosity
%         propParametersStr : Structure containing the material properties 
%                             of the structure,
%                                    .k : Spring stiffness (SDOF)
%                                    .m : Spring mass (SDOF)
%                                    .E : Young's modulus (MDOF)
%                                  .nue : Poisson ratio (MDOF)
%                                  .rho : Structural density (MDOF)
%                                    .t : Structural thickness (MDOF)
%                        up : Initial guess for the velocity/pressure field
%                         d : Initial guess for the displacement field
%                     upDot : Fluid acceleration field of the previous time
%                             step
%                      dDot : Structural velocity field of the previous 
%                             time step
%                    upDDot : Dummy variable
%                     dDDot : Structural accelration field of the previous
%                             time step
%                massMtxFld : Fluid mass matrix (dummy for the VMS-
%                             stabilized solver because the mass matrix is 
%                             nonlinear to the primal solution field due to 
%                             VMS stabilizaation)
%                massMtxStr : Structural mass matrix
%                dampMtxFld : Fluid damping matrix (dummy variable)
%                dampMtxStr : Structural damping matrix
%        precompStiffMtxFld : Constant part of the fluid tangent stiffness 
%                             matrix
%        precompStiffMtxStr : Constant part of the structural tangent 
%                             stiffness matrix
%           precomResVctFld : Constant part of the fluid residual vector
%           precomResVctStr : Constant part of the structural residual 
%                             vector
%  computeMtxSteadyStateFld : Function handle to the computation of the
%                             steady-state tangent stiffness matrix and
%                             residual vector for the fluid problem
%  computeMtxSteadyStateStr : Function handle to the computation of the
%                             steady-state tangent stiffness matrix and
%                             residual vector for the structural problem
% solve_FEMEquationSystemFld : solve_FEMLinearSystem or
%                              solver_FEMNLinearSystem
% solve_FEMEquationSystemStr : solve_FEMLinearSystem or
%                              solver_FEMNLinearSystem
%            DOFNumberingFld : Global numbering of the DOFs for the fluid
%                              problem
%            DOFNumberingStr : Global numbering of the DOFs for the
%                              structural problem
%    freeDOFsFld,freeDOFsStr : Vectors containing the free DOFs of the
%                              fluid and the structural problems
%      homDOFsFld,homDOFsStr : Vectors containing the global numbering of
%                              the DOFs where homogeneneous Dirichlet
%                              boundary conditions are applied
%  inhomDOFsFld,inhomDOFsStr : Vectors containing the global numbering of
%                              the DOFs where inhomogeneneous Dirichlet
%                              boundary conditions are applied
%         valuesInhomDOFsFld : Prescribed values of the DOFs on the 
%                              inhomogeneous Dirichlet boundary of the
%                              fluid
%         valuesInhomDOFsStr : Prescribed values of the DOFs on the 
%                              inhomogeneous Dirichlet boundary of the
%                              structure
%                    propALE : Structure containing information about the 
%                              ALE motion of the fluid along the FSI
%                              interface,
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
%                 propALEStr : Dummy variable for this function
%      computeUpdatedMeshFld : Function handle to the computation of the
%                              updated mesh corresponding to the fluid ALE
%                              formulation
%      computeUpdatedMeshStr : Dummy variable for this function
%      solve_LinearSystemFld : Function handle to the linear equation
%                              system solver for the fluid problem
%      solve_LinearSystemStr : Function handle to the linear equation
%                              system solver for the structural problem
%            propFldDynamics : Structure containing information on the 
%                              fluid dynamics,
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
%            propStrDynamics : Structure containing information on the 
%                              structural dynamics,
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
%                          t : Time instance
%     propNLinearAnalysisFld : Structure containing information on the 
%                              nonlinear fluid analysis,
%                                .method : The nonlinear solution scheme
%                                   .eps : The residual tolerance
%                               .maxIter : The maximum number of nonlinear
%                                          iterations
%     propNLinearAnalysisStr : Structure containing information on the 
%                              nonlinear structural analysis,
%                                .method : The nonlinear solution scheme
%                                   .eps : The residual tolerance
%                               .maxIter : The maximum number of nonlinear
%                                          iterations
%                    propFSI : Structure containing information on the FSI
%                              coupling algorithm,
%                               .relaxation : Under-relaxation factor
%                                      .tol : Relative residual tolerance
%                                             on the displacement field
%            propPostProcFld : Structure containing information on the
%                              computation of postprocessing resultants for
%                              the fluid problem,
%                               .nameDomain : Name of the domain onto which
%                                             postprocessing resultants are
%                                             to be computed
%                              .nodesDomain : IDs of the nodes which are on
%                                             the selected domain
%                          .computePostProc : Function handle to the
%                                             computation of the desirable 
%                                             resultant
%            propPostProcStr : Structure containing information on the
%                              computation of postprocessing resultants for
%                              the structural problem,
%                               .nameDomain : Name of the domain onto which
%                                             postprocessing resultants are
%                                             to be computed
%                              .nodesDomain : IDs of the nodes which are on
%                                             the selected domain
%                          .computePostProc : Function handle to the
%                                             computation of the desirable 
%                                             resultant
%            propGaussIntFld : Structure containing information on the
%                              numerical integration for the fluid problem,
%                              .type : 'default', 'user'
%                        .domainNoGP : Number of Gauss Points for the 
%                                      domain integration
%                      .boundaryNoGP : Number of Gauss Points for the
%                                      boundary integration
%            propGaussIntStr : Structure containing information on the
%                              numerical integration for the structural 
%                              problem,
%                              .type : 'default', 'user'
%                        .domainNoGP : Number of Gauss Points for the 
%                                      domain integration
%                      .boundaryNoGP : Number of Gauss Points for the
%                                      boundary integration
%                                tab : Tabulation when outputting
%                                      information on the command window
%                             outMsg : Enables outputting information on
%                                      the command window when chosen as
%                                      'outputEnabled'
%
%                     Output :
%                         up : The velocity/pressure fluid fields of the
%                              coupled problem at the given time step
%                          d : The displacement structural field of the
%                              coupled problem at the given time step
%                      upDot : The fluid acceleration field of the coupled
%                              system at the given time step
%                       dDot : The structural velocity field of the coupled
%                              system at the given time step
%                     upDDot : Variable containing dummy content
%                      dDDot : The structural acceleration field of the
%                              coupled system at the given time step
%                  resVctFld : The complete fluid residual vector of the
%                              coupled system at the given time step
%                  resVctStr : The complete structural residual vector of
%                              the coupled system at the given time step
%                     fldMsh : The updated fluid mesh due to the mesh
%                              motion
%                     strMsh : The updated structural mesh (only makes 
%                              sense for updated Lagrangian formulations)
%                    FFSIFld : The forces acting on the fluid FSI interface
%                              of the coupled system at the given time step
%                    FFSIStr : The forces acting on the structural FSI 
%                              interface of the coupled system at the given 
%                              time step
%               minElSizeFld : The minimum element edge size in the fluid
%                              finite element mesh
%               minElSizeStr : The minimume element edge size in the
%                              structural finite element mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the Gauss-Seidel iterations
% ->
%    1i. Solve the CSD problem
%
%   1ii. Update the displacement of the CSD problem using a relaxation method
%
%  1iii. Compute the relative error based on the displacement field of the CSD problem
%
%   1iv. Print progress message
%
%    1v. Check convergence
%
%   1vi. Update the time derivatives of the CSD solution field
%
%  1vii. Solve the mesh motion problem and update the nodal coordinates and velocities
%
% 1viii. Solve the CFD problem
%
%   1ix. Update the time derivatives of the CFD solution fields
%
%    1x. Compute the forces acting on the FSI interface
%
%   1xi. Save the solution of the CSD problem at the current Gauss-Seidel iteration
%
%  1xii. Update the coupling iteration counter
% <-
%
% 2. Check if the maximum number of Gauss-Seidel iterations has been exceeded
%
%% Function main body

%% 0. Read input

% Check if ALE motion for the fluid problem is defined
isALE = false;
if ~ischar(propALE)
    if isstruct(propALE)
        if ~isfield(propALE, 'nodes')
            error('Structure propALE must defined member variable nodes');
        else
            numNodesALE = length(propALE.nodes(:, 1));
        end
        if ~isfield(propALE, 'fctHandle')
            error('Structure propALE must defined member variable fctHandle');
        else
            numFctHandlesALE = length(propALE.fctHandle);
        end
        if ~isfield(propALE, 'isFree')
            error('Structure propALE must defined member variable isFree');
        else
            numIsFreeALE = length(propALE.isFree(:, 1));
        end
        if numNodesALE ~= numFctHandlesALE || numNodesALE ~= numIsFreeALE || ...
                numFctHandlesALE ~= numIsFreeALE
            error('Arrays propALE.nodes, propALE.fctHandle and propALE.isFree must have the same length');
        end
        isALE = true;
    end
end
if ~isALE
    error('No ALE motion is defined for the fluid problem but fluid-structure interaction analysis is selected');
end

% Initialize dummy arrays
uMeshALEStr = 'undefined';

% Save the old Cartesian coordinates of the nodes in the fluid mesh
nodesSaved = fldMsh.nodes;

% Save arrays containing the global DOF numberings
homDOFsFldSaved = homDOFsFld;
inhomDOFsFldSaved = inhomDOFsFld;
valuesInhomDOFsFldSaved = valuesInhomDOFsFld;
freeDOFsFldSaved = freeDOFsFld;

% Initialize coupling iteration counter
counterCI = 1;

%% 1. Loop over all the Gauss-Seidel iterations
while counterCI <= propFSI.maxIter + 1
    %% 1i. Solve the CSD problem
    [d, resVctStr, ~, minElSizeStr] = ...
        solve_FEMEquationSystemStr ...
        (propAnalysisStr, dSaved, dDotSaved, dDDotSaved, strMsh, FStr + FFSIStr, ...
        computeBodyForceVctStr, propParametersStr, d, dDot, dDDot, ...
        massMtxStr, dampMtxStr, precompStiffMtxStr, precomResVctStr, ...
        computeMtxSteadyStateStr, DOFNumberingStr, freeDOFsStr, ...
        homDOFsStr, inhomDOFsStr, valuesInhomDOFsStr, uMeshALEStr, ...
        solve_LinearSystemStr, propStrDynamics, t, propNLinearAnalysisStr, ...
        propGaussIntStr, strcat(tab,'\t'), '');
    
    %% 1ii. Update the displacement of the CSD problem using a relaxation method
    if counterCI > 1
        d = d*propFSI.relaxation + d_k*(1 - propFSI.relaxation);
    end
    
    %% 1iii. Compute the relative error based on the displacement field of the CSD problem
    if counterCI ~= 1
        if abs(dSaved) > propFSI.tol
            errCI = abs(d - d_k)/abs(dSaved);
        else
            errCI = abs(d - d_k);
        end
    end
    
    %% 1iv. Print progress message
    if counterCI ~= 1
        if strcmp(outMsg, 'outputEnabled')
            msgCI = sprintf(strcat([tab '\t'], '\tCoupling iteration %d with relative error %d\n'), ...
                counterCI - 1, errCI);
            fprintf(msgCI);
        end
    end
    
    %% 1v. Check convergence
    if counterCI ~= 1
        if errCI < propFSI.tol
            if strcmp(outMsg, 'outputEnabled')
                fprintf('\n');
            end
            break;
        end
    end
    
    %% 1vi. Update the time derivatives of the CSD solution field
    [dDot, dDDot] = propStrDynamics.computeUpdatedVct ...
        (d, dSaved, dDotSaved, dDDotSaved, propStrDynamics);
    
    %% 1vii. Solve the mesh motion problem and update the nodal coordinates and velocities
    propALE.propUser.u = 0;
    propALE.propUser.v = d;
    propALE.propUser.w = 0;
    [fldMsh, uMeshALEFld, homDOFsFld, inhomDOFsFld, valuesInhomDOFsFld, ...
        freeDOFsFld] = computeUpdatedMeshFld ...
        (fldMsh, homDOFsFldSaved, inhomDOFsFldSaved, valuesInhomDOFsFldSaved, freeDOFsFldSaved, ...
        nodesSaved, propALE, solve_LinearSystemFld, propFldDynamics, t);
    
    %% 1viii. Solve the CFD problem
    [up, resVctFld, ~, minElSizeFld] = ...
        solve_FEMEquationSystemFld ...
        (propAnalysisFld, upSaved, upDotSaved, upDDotSaved, fldMsh, FFld, ...
        computeBodyForceVctFld, propParametersFld, up, upDot, upDDot, ...
        massMtxFld, dampMtxFld, precompStiffMtxFld, precomResVctFld, ...
        computeMtxSteadyStateFld, DOFNumberingFld, freeDOFsFld, ...
        homDOFsFld, inhomDOFsFld, valuesInhomDOFsFld, uMeshALEFld, ...
        solve_LinearSystemFld, propFldDynamics, t, propNLinearAnalysisFld, ...
        propGaussIntFld, strcat(tab,'\t'), '');
    
    %% 1ix. Update the time derivatives of the CFD solution fields
    [upDot, upDDot] = propFldDynamics.computeUpdatedVct ...
        (up, upSaved, upDotSaved, upDDotSaved, propFldDynamics);
    
    %% 1x. Compute the forces acting on the FSI interface
    propPostProcFld = computePostProc ...
        (resVctFld, propAnalysisFld, propParametersFld, propPostProcFld);
    FFSIStr = propPostProcFld.valuePostProc{1}(2,1);
    
    %% 1xi. Save the solution of the CSD problem at the current Gauss-Seidel iteration
    d_k = d;
    
    %% 1xii. Update the coupling iteration counter
    counterCI = counterCI + 1;
end

%% 2. Check if the maximum number of Gauss-Seidel iterations has been exceeded
if counterCI == propFSI.maxIter + 2    
     if strcmp(outMsg, 'outputEnabled')
        warning('Maximum number of coupling iterations has been reached');
        fprintf('\n');
     end
end

end