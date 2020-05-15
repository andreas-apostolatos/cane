function testFluidStructureInteractionSpringMass(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the finite element formulation for the Fluid-Structure interaction
% analysis on the problem of a rigid cylinder represented by a mass which 
% is attached to a spring and whose motion is excited by a fluid flow in a
% rectangular channel.
%
% Function layout :
%
% 0. Read input
%
% 1. Parse the data from the GiD input file
%
% 2. Choose the equation system solver for the fluid problem
%
% 3. UI
%
% 4. Define the structural SDOF model
%
% 5. Solve the coupled transient Fluid-Structure Interaction problem
%
% 6. Define the expected solutions
%
% 7. Verify the results
%
%% Function main body

%% 0. Read input

% Define absolute tolerance
absTol = 1e-15;

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidStructureInteraction/';
caseNameFld = 'unitTest_flowAroundCylinderAdaptiveALE_FSI';

%% 1. Parse the data from the GiD input file
[fldMsh, homDOFsFld, inhomDOFsFld, valuesInhomDOFsFld, propALE, ~, ...
    propAnalysisFld, propParametersFld, propNLinearAnalysisFld, ...
    propFldDynamics, propGaussIntFld, propPostProcFld] = ...
    parse_FluidModelFromGid ...
    (pathToCase, caseNameFld, '');

%% 2. Choose the equation system solver for the fluid problem
if strcmp(propAnalysisFld.type, 'NAVIER_STOKES_2D')
    solve_LinearSystemFld = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysisFld.type, 'NAVIER_STOKES_3D')
    solve_LinearSystemFld = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen')
end

%% 3. UI

% Structural analysis
propAnalysisStr.type = 'SDOF';
propAnalysisStr.dir = 'y';

% Function handle to the update of the values of the DOFs on the
% inhomogeneous Dirichlet fluid boundary
updateInhomDOFsFld = 'undefined';

% Fluid-structure interaction properties
propFSI.relaxation = .3;
propFSI.tol = 1e-3;
propFSI.maxIter = 6e1;

% Choose solver for the fluid finite element equation system
solve_FEMSystemFld = @solve_FEMNLinearSystem;

% Choose function handles to the computation of the transient matrices and
% to the updates of the solution vector for the fluid problem
if strcmp(propFldDynamics.method, 'BOSSAK')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE;
else
    error('Time integration %s has not yet been implemented for the CFD system', ...
        propFldDynamics.method);
end

% Function handle to the computation of the fluid mass matrix
computeMassMtxFld = @computeMassMtx4FEMVMSStabNSE2D;

% Function handle to the computation of the steady-state matrix for the CFD
% problem
computeMtxSteadyStateFld = @computeFEMVMSStabMtxAndVct4NLinear4NSE;

% Define properties for the transient Dirichlet boundary conditions
propIDBCFld = 'undefined';

% On the computation of the body forces
computeBodyForcesFld = @computeConstantVerticalFluidBodyForceVct;

% Function handle to the computation of the initial conditions
% computeInitCndsFld = @computeNullInitialConditionsFEM4NSE;
computeInitCndsFld = @computeInitialConditionsVMS4NSEBurnedIn;

% On the writing the output function
propOutputFld.isOutput = false;
propOutputFld.writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;
propOutputFld.VTKResultFile = 'undefined';

%% 4. Define the structural SDOF model
strMsh.nodes = [1 0 0 0];
strMsh.elements = [];
homDOFsStr = [];
inhomDOFsStr = [];
valuesInhomDOFsStr = [];
updateInhomDOFsStr = 'undefined';
propParametersStr.k = 1e0;
propParametersStr.m = 3e1;
computeBodyForcesStr = 'undefined';
solve_LinearSystemStr = @solve_LinearSystemMatlabBackslashSolver;
solve_FEMSystemStr = @solve_FEMLinearSystem;
propStrDynamics.method = 'EXPLICIT-EULER';
if strcmp(propStrDynamics.method, 'EXPLICIT-EULER')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsEESDOF;
    propStrDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldEESDOF;
else
    error('Time integration %s has not yet been implemented for the SDOF system', ...
        propStrDynamics.method);
end
propStrDynamics.timeDependence =  propFldDynamics.timeDependence;
propStrDynamics.T0 = propFldDynamics.T0;
propStrDynamics.TEnd = propFldDynamics.TEnd;
propStrDynamics.noTimeSteps = propFldDynamics.noTimeSteps;
propStrDynamics.dt = propFldDynamics.dt;
computeMassMtxStr = @computeMassMtxSDOF;
computeMtxSteadyStateStr = @computeStiffMtxSDOF;
propNLinearAnalysisStr = 'undefined';
propIDBCStr = 'undefined';
propGaussIntStr = 'undefined';
computeInitCndsStr = @computeNullInitialConditionsSDOF;
propOutputStr.isOutput = false;
propOutputStr.VTKResultFile = 'undefined';
caseNameStr = 'sdof';
propPostProcStr = 'undefined';

%% 5. Solve the coupled transient Fluid-Structure Interaction problem
[upHistory, dHistory, resVctHistoryFld, resVctHistoryStr, ...
    minElSizeFld, ~] = ...
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
    '');

%% 6. Define the expected solutions

% Load the expected solutions
data_expSol = ...
    importdata('../../unitTest/@testFEMComputationalFluidStructureInteractionAnalysis/data_unitTest_expSolFSI.mat');

% Define the expected solution in terms of the velocity/pressure nodal
% solution field
expSolUpHistory = data_expSol.upHistory;

% Define the expected solution in terms of the displacement
expSolDHistory = data_expSol.dHistory;
           
% Define the expected solution in terms of the right-hand side residual
% field corresponding to the fluid problem
expSolResVctHistoryFld = data_expSol.resVctHistoryFld;

% Define the expected solution in terms of the right-hand side residual
% field corresponding to the structural problem
expSolResVctHistoryStr = data_expSol.resVctHistoryStr;

% Define the expected solution in terms of the minium element edge size in
% the fluid mesh
expSolMinElSizeFld = data_expSol.minElSizeFld;

%% 7. Verify the results
testCase.verifyEqual(upHistory, expSolUpHistory, 'AbsTol', absTol);
testCase.verifyEqual(dHistory, expSolDHistory, 'AbsTol', absTol);
testCase.verifyEqual(resVctHistoryFld, expSolResVctHistoryFld, 'AbsTol', absTol);
testCase.verifyEqual(resVctHistoryStr, expSolResVctHistoryStr, 'AbsTol', absTol);
testCase.verifyEqual(minElSizeFld, expSolMinElSizeFld, 'AbsTol', absTol);

end

%% Custom functions
function [tanMtx, resVct] = computeProblemMtrcsEESDOF ...
    (d, dSaved, dDot, dDotSaved, dDDot, dDDotSaved, massMtx, damMtx, tanMtx, ...
    resVct, propTransientAnalysis)
    
    % Update the tangent stiffness matrix for the SDOF system corresponding
    % to the Explicit-Euler time integration method
    tanMtx = tanMtx + massMtx*(1/propTransientAnalysis.dt^2);
    
    % Update the residual vector of the SDOF system system corresponding to
    % the Explicit-Euler time integration method
    resVct = resVct + massMtx*dSaved*(1/propTransientAnalysis.dt^2);
end
function  [dDot, dDDot] = computeBossakTIUpdatedVctAccelerationFieldEESDOF ...
    (d, dSaved, dDotSaved, dDDotSaved, propTransientAnalysis)

    % This function is dummy
    dDot = 'undefined';
    dDDot = 'undefined';

end
function massMtx = computeMassMtxSDOF ...
    (propAnalysis, strMsh, PropParameters, propGaussInt)

    % For the SDOF system the mass matrix is simply its mass
    massMtx = PropParameters.m;

end
function [stiffMtx, resVct, minElSize] = computeStiffMtxSDOF ...
    (propAnalysis, d, dSaved, dDot, dDotSaved, uMeshALE, ...
    precompStiffMtx, precomResVct, DOFNumbering, strMsh, F, ...
    loadFactor, computeBodyForces, propStrDynamics, t, ...
    propParameters, propGaussInt)

    % For the SDOF system the stiffness matrix is simply the spring 
    % stiffness
    stiffMtx = propParameters.k;
    
    % No body forces are assumed
    resVct = F;
    
    % A point has no size
    minElSize = 'undefined';

end
function [up, upDot, upDDot, numTimeStep] = ...
    computeInitialConditionsVMS4NSEBurnedIn ...
    (propAnalysis, fldMsh, DOF4Output, propParameters, propFlsDynamics, ... 
    VTKResultFile, caseName, pathToFile)
    
    path_prefix = '../../';
    data_burn_in_analysis = ...
        load([path_prefix 'preComputedData/FEMComputationalFluidDynamicsAnalysis/data_' caseName]);
    up = data_burn_in_analysis.u;
    upDot = data_burn_in_analysis.uDot;
    upDDot = data_burn_in_analysis.uDDot;
    numTimeStep = 0;
    
end
function [d, dDot, dDDot, numTimeStep] = ...
    computeNullInitialConditionsSDOF ...
    (propAnalysisStr, strMsh, DOF4Output, propParameters, propStrDynamics, ... 
    VTKResultFile, caseName, pathToFile)

    d = 0;
    dDot = 'undefined';
    dDDot = 'undefined';
    numTimeStep = 0;

end
