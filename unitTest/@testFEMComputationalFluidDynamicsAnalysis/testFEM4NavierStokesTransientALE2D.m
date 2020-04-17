function testFEM4NavierStokesTransientALE2D(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% Tests the ALE formulation and the mesh motion algorithm when using the
% pseudostructural solver for the Variational Multi-Scale Finite Element 
% formulation of the 2D incmpressible Navier-Stokes equation
%
% Function layout :
%
% 0. Read input
%
% 1. Parse the data from the GiD input file
%
% 2. UI
%
% 3. Choose the equation system solver
%
% 4. Solve the CFD problem
%
% 5. Calculate drag and lift force for each time step
%
% 6. Define the expected solutions
%
% 7. Verify the results
%
%% Function main body

%% 0. Read input

% Define absolute tolerance
absTol = 1e-15; % tolerance w.r.t our value

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_ALE';

%% 1. Parse the data from the GiD input file
[fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propALE, ~, propAnalysis, ...
    parameters, propNLinearAnalysis, propFldDynamics, propGaussInt, ...
    propPostproc] = ...
    parse_FluidModelFromGid...
    (pathToCase, caseName, '');

%% 2. UI
% On the body forces
computeBodyForces = @computeConstantVerticalFluidBodyForceVct;

% Properties for the VTK visualization
propVTK.isOutput = false;
propVTK.VTKResultFile = 'undefined';

% On the initial conditions
% computeInitialConditions = @computeInitialConditionsFromVTKFileFEM4NSE2D;
computeInitialConditions = @computeNullInitialConditionsFEM4NSE;

% On the transient analysis properties
if strcmp(propFldDynamics.method, 'BOSSAK')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE;
end

% On transient inhomogeneous Dirichlet boundary conditions
updateInhomDOFs = 'undefined';
propIDBC = [];

%% 3. Choose the equation system solver
if strcmp(propAnalysis.type, 'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysis.type, 'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
end

%% 4. Solve the CFD problem
[~, FHistory, minElSize] = solve_FEMVMSStabTransientNSEBossakTI ...
    (fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, updateInhomDOFs, ...
    propALE, parameters, computeBodyForces, propAnalysis, ...
    computeInitialConditions, solve_LinearSystem, propFldDynamics, ...
    propNLinearAnalysis, propIDBC, propGaussInt, propVTK, caseName, ...
    '');

%% 5. Calculate drag and lift force for each time step
forcesOnCylinder = zeros(propFldDynamics.noTimeSteps + 1, 2);
for iTimeStep = 1:propFldDynamics.noTimeSteps + 1
    propPostproc = computePostProc ...
        (FHistory(:, iTimeStep), propAnalysis, parameters, propPostproc);
    forcesOnCylinder(iTimeStep, :) = propPostproc.valuePostProc{1}'; 
end
    
%% 6. Define the expected solutions

% Define the expected solution in terms of the forces acting on the
% cylinder
expSolForcesOnCylinder = [ 0                   0
                          -0.015679093695312  -0.000045980271542
                           0.027009475069121  -0.000023780449891
                           0.001235463820661   0.000007137558898
                          -0.019415271861986  -0.000125997913955];
                      
% Define the expected solution in terms of the minimum element area size
expSolMinElSize = 0.002;

%% 7. Verify the results
testCase.verifyEqual(forcesOnCylinder, expSolForcesOnCylinder, 'AbsTol', absTol);
testCase.verifyEqual(minElSize, expSolMinElSize);

end