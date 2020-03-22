function [propALE,lift,drag] = nominalSolution(fldMsh,up,homDOFs,inhomDOFs,valuesInhomDBCModified,...
        parameters,computeBodyForces,analysis,computeInitialConditions,...
        VTKResultFile,propALE,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
        i,propVTK_true,gaussInt,caseName,p1,p2,p3,postProc)
%% Licensing
%
%  License:         BSD License
%                   cane Multiphysics default license: cane/license.txt
%
%  Main authors:    Matthew Keller
%
%% Function documentation
%
%  Returns the lift and drag of the nominal mesh
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
%          parameters : Flow parameters
%   computeBodyForces : Function handle to the computation of the body
%                       force vector
%            analysis : .type : The analysis type
%     computeInitCnds : Function handle to the initial boundary conditions 
%                       computation
%       VTKResultFile : The name of the result file in the output folder
%                        where to get the initial conditions for the
%                        transient simulation
%             propALE : Properties regarding the ALE boundary
%                         .nodes : The sequence of the nodal coordinates
%                                  on the ALE boundary
%                     .fcthandle : Function handle to the computation of
%                                  the ALE motion
%                       propUser : Extra user-defined parameters
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
%                   i : The current optimization iteration step
%            gaussInt : On the Gauss Point integration
%                             .type : 'default', 'user'
%                       .domainNoGP : Number of Gauss Points for the domain
%                                     integration
%                     .boundaryNoGP : Number of Gauss Points for the
%                                     boundary integration
%            caseName : String defining the case name
%            p1,p2,p3 : Internal variables from previous step
%            postProc : Post processing of force values across the body
%
%              Output :
%             propALE : Extra user-defined parameters
%                Lift : Lift on the body
%                Drag : Drag on the body

%% Function main body

% Update internal variables with updated values from previous iteration
propALE.propUser.p1 = p1;
propALE.propUser.p2 = p2;
propALE.propUser.p3 = p3;

% Solve the CFD problem in nominal state   
[~,FComplete,~,~] = solve_FEMVMSStabSteadyStateNSE2D...
    (fldMsh,up,homDOFs,inhomDOFs,valuesInhomDBCModified,'undefined',parameters,...
    computeBodyForces,analysis,computeInitialConditions,...
    VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
    i,propVTK_true,gaussInt,caseName,'outputEnabled');

% Calculate drag and lift force from the nodal forces
postProc_update = computePostProc(FComplete,analysis,parameters,postProc);

% Retrieve Fx and Fy from post processing
forcesOnDomain = postProc_update.valuePostProc{1};
Fx = forcesOnDomain(1,1);
Fy = forcesOnDomain(2,1);

% Set drag and lift
lift = Fy;
drag = Fx;
end    