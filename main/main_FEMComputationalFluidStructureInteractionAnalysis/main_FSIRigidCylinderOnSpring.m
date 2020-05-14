%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Solves the transient incompressible Navier-Stokes equations
%
% Date : 02.05.2020
%
%% Preamble
clear;
clc;
close all;

%% Includes

% path prefix
path_prefix = '../../';

% Add transient analysis functions
addpath([path_prefix 'transientAnalysis/']);

% Add functions related to equation system solvers
addpath([path_prefix 'equationSystemSolvers/']);

% Add general math functions
addpath([path_prefix 'generalMath/']);

% Add the classical finite element basis functions
addpath([path_prefix 'basisFunctions/']);

% Add all functions related to plate in membrane action analysis
addpath([path_prefix 'FEMPlateInMembraneActionAnalysis/solvers/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/loads/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/graphics/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/output/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/postprocessing/']);

% Add all functions related to the Finite Element Methods for Computational
% Fluid Dynamics problems
addpath([path_prefix 'FEMComputationalFluidDynamicsAnalysis/solutionMatricesAndVectors/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/initialConditions'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/solvers/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/loads/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/ALEMotion/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/transientAnalysis/'],...
        [path_prefix 'FEMComputationalFluidDynamicsAnalysis/postProcessing/']);

% Add all functions related to Fluid-Structure interaction
addpath([path_prefix 'FEMComputationalFluidStructureInteractionAnalysis/ALEMotion'], ...
        [path_prefix 'FEMComputationalFluidStructureInteractionAnalysis/solvers']);
    
% Add all functions related to parsing
addpath([path_prefix 'parsers/']);

% Add all functions related to the efficient computation functions
addpath([path_prefix 'efficientComputation/']);

%% Parse the data from the GiD input file
pathToCase = '../../inputGiD/FEMComputationalFluidStructureInteraction/';
caseNameFld = 'flowAroundCylinderAdaptiveALE_FSI';
[fldMsh, homDOFsFld, inhomDOFsFld, valuesInhomDOFsFld, propALE, ~, ...
    propAnalysisFld, propParametersFld, propNLinearAnalysisFld, ...
    propFldDynamics, propGaussIntFld, propPostProcFld] = ...
    parse_FluidModelFromGid ...
    (pathToCase, caseNameFld, 'outputEnabled');

%% Choose the equation system solver for the fluid problem
if strcmp(propAnalysisFld.type, 'NAVIER_STOKES_2D')
    solve_LinearSystemFld = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysisFld.type, 'NAVIER_STOKES_3D')
    solve_LinearSystemFld = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen')
end

%% UI

% Structural analysis
propAnalysisStr.type = '2DOF'; %'SDOF'
if strcmp(propAnalysisStr.type , 'SDOF')
    propAnalysisStr.dir = 'y';
end

% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMFluidStructureInteraction/';

% Function handle to the update of the values of the DOFs on the
% inhomogeneous Dirichlet fluid boundary
updateInhomDOFsFld = 'undefined';

% Fluid-structure interaction properties
propFSI.relaxation = .3;
propFSI.tol = 1e-3;
propFSI.maxIter = 1e2;

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
propOutputFld.isOutput = true;
propOutputFld.writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;
propOutputFld.VTKResultFile = 'undefined'; % '_contourPlots_75'

%% Define the structural SDOF model
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
        @computeProblemMtrcsEESpringMassSystem;
    propStrDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldEESpringMassSystem;
else
    error('Time integration %s has not yet been implemented for the SDOF system', ...
        propStrDynamics.method);
end
propStrDynamics.timeDependence =  propFldDynamics.timeDependence;
propStrDynamics.T0 = propFldDynamics.T0;
propStrDynamics.TEnd = propFldDynamics.TEnd;
propStrDynamics.noTimeSteps = propFldDynamics.noTimeSteps;
propStrDynamics.dt = propFldDynamics.dt;
if strcmp(propAnalysisStr.type, 'SDOF')
    computeMassMtxStr = @computeMassMtxSDOF;
elseif strcmp(propAnalysisStr.type, '2DOF')
    computeMassMtxStr = @computeMassMtx2DOF;
else
    error('Analysis type %s is not yet supported', propAnalysisStr.type);
end
if strcmp(propAnalysisStr.type, 'SDOF')
    computeMtxSteadyStateStr = @computeStiffMtxSDOF;
elseif strcmp(propAnalysisStr.type, '2DOF')
    computeMtxSteadyStateStr = @computeStiffMtx2DOF;
else
    error('Analysis type %s is not yet supported', propAnalysisStr.type);
end
propNLinearAnalysisStr = 'undefined';
propIDBCStr = 'undefined';
propGaussIntStr = 'undefined';
if strcmp(propAnalysisStr.type, 'SDOF')
    computeInitCndsStr = @computeNullInitialConditionsSDOF;
elseif strcmp(propAnalysisStr.type, '2DOF')
    computeInitCndsStr = @computeNullInitialConditions2DOF;
else
    error('Analysis type %s is not yet supported', propAnalysisStr.type);
end
propOutputStr.isOutput = false;
propOutputStr.VTKResultFile = 'undefined';
caseNameStr = 'sdof';
propPostProcStr = 'undefined';

%% Solve the coupled transient Fluid-Structure Interaction problem
[upHistory, dHistory, resVctHistoryFld, resVctHistoryStr, ...
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
    'outputEnabled');

%% Custom functions
function [tanMtx, resVct] = computeProblemMtrcsEESpringMassSystem ...
    (d, dSaved, dDot, dDotSaved, dDDot, dDDotSaved, massMtx, damMtx, tanMtx, ...
    resVct, propTransientAnalysis)
    
    % Update the tangent stiffness matrix for the SDOF system corresponding
    % to the Explicit-Euler time integration method
    tanMtx = tanMtx + massMtx*(1/propTransientAnalysis.dt^2);
    
    % Update the residual vector of the SDOF system system corresponding to
    % the Explicit-Euler time integration method
    resVct = resVct + massMtx*dSaved*(1/propTransientAnalysis.dt^2);
end
function  [dDot, dDDot] = computeBossakTIUpdatedVctAccelerationFieldEESpringMassSystem ...
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
function massMtx = computeMassMtx2DOF ...
    (propAnalysis, strMsh, PropParameters, propGaussInt)

    % For the SDOF system the mass matrix is simply its mass
    massMtx = eye(2, 2)*PropParameters.m;

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
function [stiffMtx, resVct, minElSize] = computeStiffMtx2DOF ...
    (propAnalysis, d, dSaved, dDot, dDotSaved, uMeshALE, ...
    precompStiffMtx, precomResVct, DOFNumbering, strMsh, F, ...
    loadFactor, computeBodyForces, propStrDynamics, t, ...
    propParameters, propGaussInt)

    % For the SDOF system the stiffness matrix is simply the spring 
    % stiffness
    stiffMtx = [1e5 0
                0  1]*propParameters.k;
    
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
function [d, dDot, dDDot, numTimeStep] = ...
    computeNullInitialConditions2DOF ...
    (propAnalysisStr, strMsh, DOF4Output, propParameters, propStrDynamics, ... 
    VTKResultFile, caseName, pathToFile)

    d = zeros(2, 1);
    dDot = 'undefined';
    dDDot = 'undefined';
    numTimeStep = 0;

end

%% END OF SCRIPT
