function testFEM4NavierStokesSteadyStateFlowAroundCylinder2D(testCase)
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
% Tests the stabilized finite element analysis for the 2D steady-state
% Navier-Stokes equations, over the flow around a cylinder problem with 
% respect to reference paper by:
% http://www.mathematik.tu-dortmund.de/lsiii/cms/papers/SchaeferTurek1996.pdf
%
% Function layout :
%
% 0. Read input
%
% 1. Parse the data from the GiD input file
%
% 2. GUI
%
% 3. Choose the equation system solver
%
% 4. Define the name of the vtk file from where to resume the simulation
%
% 5. Change the input velocity to match the reference paper - parabolic input
%
% 6. Solve the CFD problem
%
% 7. Calculate drag and lift force from the nodal forces
%
% 8. Calculate drag and lift coefficient based on drag and lift force
%
% 9. Define the expected solutions
%
% 10. Verify the results
%
%% Function main body

%% 0. Read input
% Define absolute tolerance
absTol = 0.04; % tolerance w.r.t the reference paper
absTol2 = 1e-10; % tolerance w.r.t our value

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_testFEM4NavierStokesSteadyStateFlowAroundCylinder2D';

%% 1. Parse the data from the GiD input file
[fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,~,analysis,parameters,...
    propNLinearAnalysis,propFldDynamics,gaussInt,postProc] = ...
    parse_FluidModelFromGid...
    (pathToCase,caseName,'');

%% 2. GUI
% On the body forces
computeBodyForces = @computeConstantVerticalBodyForceVct;

% On the initial conditions
% computeInitialConditions = @computeInitialConditionsFromVTKFileFEM4NSE2D;
computeInitialConditions = @computeNullInitialConditionsFEM4NSE2D;

% On the transient analysis properties
if strcmp(propFldDynamics.method,'bossak')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE2D;
end

%% 3. Choose the equation system solver
if strcmp(analysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
end

%% 4. Define the name of the vtk file from where to resume the simulation
VTKResultFile = 'undefined';

%% 5. Change the input velocity to match the reference paper - parabolic input

% max input velocity defined in the reference paper
Umax = 0.3;

% change the input velocity to have the parabolic distribution
valuesInhomDBCModified = computeInletVelocityParabolic_unitTest(fldMsh, inhomDBC, valuesInhomDBC, Umax);

%% 6. Solve the CFD problem
[~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
    (fldMsh,homDBC,inhomDBC,valuesInhomDBCModified,nodesALE,parameters,...
    computeBodyForces,analysis,computeInitialConditions,...
    VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
    gaussInt,caseName,'');

%% 7. Calculate drag and lift force from the nodal forces
postProc = computePostProc(FComplete, analysis, parameters, postProc);

%% 8. Calculate drag and lift coefficient based on drag and lift force

% define parameters used in reference paper and simualiton
Ubar = 0.2; % mid velocity 
D = 0.1;    % diameter of the body
rho = parameters.rho; % density

% get Fx and Fy from post processing
forcesOnDomain = postProc.valuePostProc{1};
Fx = forcesOnDomain(1,1);
Fy = forcesOnDomain(2,1);

% calculate drag and lift coefficiet
dragCoefficient = (2 * Fx)/(rho * Ubar * Ubar * D);
liftCoefficient = (2 * Fy)/(rho * Ubar * Ubar * D);

% find absolute value so we don't get negative coefficients
liftCoefficient = abs(liftCoefficient);

%% 9. Define the expected solutions

% 5.5700 - 5.5900 (upper and lower bound in the paper, 2D case)
expSolDragCoefficientFromLiterature = 5.58;
expSolDragCoefficient = 5.61079238629658;

% 0.0104 - 0.0110 (upper and lower bound in the paper, 2D case)
% in reality it should be 0, since it's a symmetric case
expSolLiftCoefficientFromLiterature = 0.0107;
expSolLiftCoefficient = 0.00224889649350511;

% Define the expected solution in terms of the convergence flag
expSolHasConverged = true;

%% 10. Verify the results
testCase.verifyEqual(dragCoefficient,expSolDragCoefficientFromLiterature,'AbsTol',absTol);
testCase.verifyEqual(liftCoefficient,expSolLiftCoefficientFromLiterature,'AbsTol',absTol);
testCase.verifyEqual(dragCoefficient,expSolDragCoefficient,'AbsTol',absTol2);
testCase.verifyEqual(liftCoefficient,expSolLiftCoefficient,'AbsTol',absTol2);
testCase.verifyEqual(hasConverged,expSolHasConverged,'AbsTol',absTol);

end
