%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Solves the coupled Fluid-Structure Interaction problem for the
%        Turek FSI benchmark.
%
% Date : 17.05.2020
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
        [path_prefix 'FEMPlateInMembraneActionAnalysis/postprocessing/'],...
        [path_prefix 'FEMPlateInMembraneActionAnalysis/initialConditions/']);

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

%% Case name
pathToCase = '../../inputGiD/FEMComputationalFluidStructureInteraction/';
% caseName = 'turek_fsi';
caseName = 'building_fsi';

%% Parse the fluid setup from the GiD input file
caseNameFld = caseName;
[fldMsh, homDOFsFld, inhomDOFsFld, valuesInhomDOFsFld, propALE, propNBCFld, ...
    propAnalysisFld, propParametersFld, propNLinearAnalysisFld, ...
    propFldDynamics, propGaussIntFld, propPostProcFld, propFSIFld] = ...
    parse_FluidModelFromGid ...
    (pathToCase, caseNameFld, 'outputEnabled');
caseNameFld = horzcat(caseNameFld, '_fluid');

%% Parse the structural setup from the GiD input file
caseNameStr = caseName;
[strMsh, homDOFsStr, inhomDOFsStr, valuesInhomDOFsStr, propNBCStr, ...
    propAnalysisStr, propParametersStr, propNLinearAnalysisStr, ...
    propStrDynamics, propGaussIntStr, ~, propFSIStr] = ...
    parse_StructuralModelFromGid...
    (pathToCase, caseNameStr, 'outputEnabled');
caseNameStr = horzcat(caseNameStr, '_structure');

%% Apply a non-constant inlet
if strcmp(caseNameFld, 'turek_fsi_fluid') || ...
        strcmp(caseNameFld, 'building_fsi_fluid')
        
    % Define number of DOFs per node
    noDOFsPerNode = 3;
    
    % Define the law for the inlet velocity
    if strcmp(caseNameFld, 'turek_fsi_fluid')
        % Parabolic law
        u_max = 2;
        computeInletVelocity = @(x,y,z) u_max*6*y*(0.41 - y)/0.1681;
    elseif strcmp(caseNameFld, 'building_fsi_fluid')
        % 1/7-th power law
        u_max = 5;
        computeInletVelocity = @(x,y,z) u_max*y^(1/7);
    end
    
    % Loop over all inlet DOFs
    for i = 1:length(inhomDOFsFld)
        % Find the inlet DOF
        idxDOF = inhomDOFsFld(1, i);

        % Find the corresponding node
        idxNode = ceil(idxDOF/noDOFsPerNode);

        % Get the nodal coordinates
        XYZ = fldMsh.nodes(idxNode, 2:end);

        % Compute the value according to the law
        presValue = computeInletVelocity(XYZ(1, 1), XYZ(1, 2), XYZ(1, 3));

        % Cartesian direction
        cartDir = idxDOF - (noDOFsPerNode*ceil(idxDOF/noDOFsPerNode) - noDOFsPerNode);

        if cartDir == 1
            valuesInhomDOFsFld(1, i) = presValue;
        end
    end
end

%% UI

% Choose the equation system solver for the fluid and the structural problems
if strcmp(propAnalysisFld.type, 'NAVIER_STOKES_3D')
    error('Fluid-Structure Interaction analysis is not implemented for 3D fluid flows');
end
solve_LinearSystemFld = @solve_LinearSystemMatlabBackslashSolver;
solve_LinearSystemStr = @solve_LinearSystemMatlabBackslashSolver;

% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMFluidStructureInteraction/';

% Fluid-structure interaction properties
propFSI.relaxationMethod = 'constant'; % 'constant', 'aitken'
propFSI.relaxation = 3e-1;
propFSI.tol = 1e-5;
propFSI.maxIter = 1e2;

%% Function handle to the solution procedure of the coupled system
solve_FEMCoupledSystem = ...
    @solve_FEMFluidStructureInteractionPartitionedGaussSeidel;

%% Fluid dynamics

% Choose solver for the fluid finite element equation system
solve_FEMSystemFld = @solve_FEMNLinearSystem;

% Function handle to the update of the values of the DOFs on the
% inhomogeneous Dirichlet fluid boundary
updateInhomDOFsFld = 'undefined';

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
computeInitCndsFld = @computeInitialConditionsFromVTKFileFEM4NSEWrapper;
% computeInitCndsFld = @computeInitialConditionsFromVTKFileFEM4NSE;
% computeInitCndsFld = @computeInitialConditionsVMS4NSEBurnedIn;

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

% On the writing the output function
propOutputFld.isOutput = true;
propOutputFld.writeOutputToFile = @writeOutputFEMIncompressibleFlowToVTK;
propOutputFld.VTKResultFile = '_contourPlots_100'; % '_contourPlots_21'

%% Structural dynamics

% Choose solver for the fluid finite element equation system
solve_FEMSystemStr = @solve_FEMNLinearSystem;

% Function handle to the update of the values of the DOFs on the
% inhomogeneous Dirichlet structural boundary
updateInhomDOFsStr = 'undefined';

% Function handle to the computation of the fluid mass matrix
computeMassMtxStr = @computeMassMtxFEMPlateInMembraneActionCST;

% Function handle to the computation of the steady-state matrix for the CFD
% problem
computeMtxSteadyStateStr = @computeTangentStiffMtxResVctFEMPlateInMembraneActionCST;

% Define properties for the transient Dirichlet boundary conditions
propIDBCStr = 'undefined';

% On the computation of the body forces
computeBodyForcesStr = @computeConstantVerticalStructureBodyForceVct;

% Function handle to the computation of the initial conditions
computeInitCndsStr = @computeInitCndsFEMPlateInMembraneAction;
% computeInitCndsStr = ...
%     @computeInitialConditionsFromVTKFileFEMPlateInMembraneAction;

% Choose function handles to the computation of the transient matrices and
% to the updates of the solution vector for the fluid problem
if strcmp(propStrDynamics.method, 'EXPLICIT_EULER')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsExplicitEuler;
    propStrDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propStrDynamics.method, 'BOSSAK')
    propStrDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossak;
    propStrDynamics.computeUpdatedVct = ...
        @computeBossakTransientUpdatedVctAccelerationField;
else
    error('Invalid time integration method selected in propStrDynamics.method as %s', ...
        propStrDynamics.method);
end

% On the writing the output function
propOutputStr.isOutput = true;
propOutputStr.writeOutputToFile = @writeOutputFEMPlateInMembraneActionToVTK;
propOutputStr.VTKResultFile = '_contourPlots_139'; % '_contourPlots_2'

% Postprocessing properties for the CSD problem
propPostProcStr = 'undefined';

%% Solve the coupled transient Fluid-Structure Interaction problem
[upHistory, dHistory, resVctHistoryFld, resVctHistoryStr, ...
    minElSizeFld, minElSizeStr] = ...
    solve_FEMTransientFluidStructureInteraction ...
    (propAnalysisFld, propAnalysisStr, fldMsh, strMsh, homDOFsFld, homDOFsStr, ...
    inhomDOFsFld, inhomDOFsStr, valuesInhomDOFsFld, valuesInhomDOFsStr, ...
    updateInhomDOFsFld, updateInhomDOFsStr, propALE, propParametersFld, ...
    propParametersStr, computeBodyForcesFld, computeBodyForcesStr, ...
    propNBCFld, propNBCStr, computeInitCndsFld, computeInitCndsStr, ...
    solve_FEMSystemFld, solve_FEMSystemStr, solve_FEMCoupledSystem, ...
    computeMassMtxFld, computeMassMtxStr, computeMtxSteadyStateFld, ...
    computeMtxSteadyStateStr, solve_LinearSystemFld, solve_LinearSystemStr, ...
    propFldDynamics, propStrDynamics, propNLinearAnalysisFld, ...
    propNLinearAnalysisStr, propFSIFld, propFSIStr, propFSI, propIDBCFld, ...
    propIDBCStr, propPostProcFld, propPostProcStr, propGaussIntFld, ...
    propGaussIntStr, propOutputFld, propOutputStr, caseNameFld, caseNameStr, ...
    'outputEnabled');

%% Custom function
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

function [up, upDot, upDDot, numTimeStep] = ...
    computeInitialConditionsFromVTKFileFEM4NSEWrapper ...
    (analysis, fldMsh, DOF4Output, parameters, fldDynamics, ...
    VTKResultFile, caseName, pathToFileDummy)

    caseName = erase(caseName, '_fluid');
    
    pathToFile = '../../outputVTK/FEMComputationalFluidDynamicsAnalysis/';

    [up, upDot, upDDot, numTimeStep] = ... 
        computeInitialConditionsFromVTKFileFEM4NSE ...
        (analysis, fldMsh, DOF4Output, parameters, fldDynamics, ...
        VTKResultFile, caseName, pathToFile);
    
    numTimeStep = 0;

end

%% END OF SCRIPT