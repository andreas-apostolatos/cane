function testFEM4TransientTaylorGreenVortices2D(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the stabilized finite element analysis for the 2D transient Navier-
% Stokes equations for the benchmark problem of the Taylor-Green vortices. 
% For this problem there is an analytical solution as a function of time 
% and space.
%
% Function layout :
%
% 0. Read input
%
% 1. Parse the data from the GiD input file
%
% 2. GUI
%
% 3. Define the transient inhomogeneous Dirichlet boundary conditions
%
% 4. Choose the equation system solver
%
% 5. Solve the CFD problem
%
% 6. Compute the 2-norm of the velocity at given Cartesian location over time
%
% 7. Compute the pressure at given Cartesian location over time
%
% 8. Define the expected solutions
%
% 9. Verify the results
%
%% Function main body

%% 0. Read input

% Define absolute tolerance
absTol = 1e-15;

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_taylorGreenVortices_2pi_domain';

%% 1. Parse the data from the GiD input file
% Parse the data
[fldMsh, homDOFs, inhomDOFs, ~, nodesALE, ~, propAnalysis, parameters, ...
    propNLinearAnalysis, propFldDynamics, propGaussInt] = ....
    parse_FluidModelFromGid...
    (pathToCase, caseName, '');

%% 2. GUI

% On the computation of the body forces
computeBodyForces = @computeConstantVerticalFluidBodyForceVct;

% On the writing the output function
propVTK.isOutput = false;
propVTK.writeOutputToFile = 'undefined';
propVTK.VTKResultFile = 'undefined'; % '_contourPlots_75'

% On the initial conditions
computeInitialConditions = ...
    @computeInitialConditionsForTaylorGreenVorticesFEM4NSE2D;

% On the transient analysis properties
if strcmp(propFldDynamics.method, 'BOSSAK')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE;
end

% On the postprocessing properties
%'xVelocity','yVelocity','pressure','2normVelocity','velocityVectorPlot'
propPostprocVelocity.postProcComponent = '2normVelocity';
propPostprocVelocity.computeAnalytical = 'undefined';
propPostprocPressure.postProcComponent = 'pressure';
propPostprocPressure.computeAnalytical = 'undefined';

%% 3. Define the transient inhomogeneous Dirichlet boundary conditions
computeTaylorGreenBCs = @(fldMsh,propIDBC,t) reshape([-cos(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)).*sin(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2))*exp(-2*t*propIDBC.nue),...
                                                       sin(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)).*cos(fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2))*exp(-2*t*propIDBC.nue),...
                                                      -0.25*( cos(2*fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),1)) + cos(2*fldMsh.nodes(unique(ceil(inhomDOFs./propAnalysis.noFields)),2)) )*exp(-4*t*propIDBC.nue)]',1,[]);
% Assign anonymous function as a function handle
updateInhomDOFs = computeTaylorGreenBCs;
propIDBC = [];

% Assign the taylor-Green boundary conditions function
valuesInhomDOFs = computeTaylorGreenBCs(fldMsh,parameters,propFldDynamics.T0);

%% 4. Choose the equation system solver
if strcmp(propAnalysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen')
end

%% 5. Solve the CFD problem
[upHistory, minElSize] = solve_FEMVMSStabTransientNSEBossakTI ...
    (fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, updateInhomDOFs, ...
    nodesALE, parameters, computeBodyForces, propAnalysis, ...
    computeInitialConditions, solve_LinearSystem, propFldDynamics, ...
    propNLinearAnalysis, propIDBC, propGaussInt, propVTK, caseName, ...
    '');

%% 6. Compute the 2-norm of the velocity at given Cartesian location over time
[timeSpaceDiscrete, velocity2NormNumerical, ~] = ...
    computeResultantAtPointOverTime...
    (pi/3, pi/3, fldMsh, parameters, upHistory, ...
    propFldDynamics, propPostprocVelocity, '');

%% 7. Compute the pressure at given Cartesian location over time
[~, pressureNumerical, ~] = ...
    computeResultantAtPointOverTime...
    (pi/2, -pi/2, fldMsh, parameters, upHistory, ...
    propFldDynamics, propPostprocPressure, '');


%% 8. Define the expected solutions

% Define the expected solution in terms of the time discretization
expTimeSpaceDiscrete = [   0
                           0.500000000000000
                           1.000000000000000];

% Define the expected solution in terms of the 2-norm of the velocity at
% (pi/3, pi/3) over time
expVelocity2NormNumerical = [0.602187820036510
                             0.601105320349726
                             0.599657792178951];
                         
% Define the expected solution in terms of the pressure at (pi/2, -pi/2) 
% over time
expPressureNumerical = [0.499999288000128
                        0.495393690707311
                        0.495066261748330];

% Define the expected solution in terms of the minimum element size
expSolMinElEdgeSize = 0.196340000000000;

%% 9. Verify the results
testCase.verifyEqual(timeSpaceDiscrete, expTimeSpaceDiscrete, 'AbsTol', absTol);
testCase.verifyEqual(velocity2NormNumerical, expVelocity2NormNumerical, 'AbsTol', absTol);
testCase.verifyEqual(pressureNumerical, expPressureNumerical, 'AbsTol', absTol);
testCase.verifyEqual(minElSize, expSolMinElEdgeSize, 'AbsTol', absTol);

end