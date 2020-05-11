function testTransientWallHeating(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the transient thermal conduction analysis over the benchmark case
% of the transient wall heating
%
% Function layout :
%
% 0. Read input
%
% 1. Parse data from GiD input file
%
% 2. UI
%
% 3. Define intial temperature and applied heat flux
%
% 4. Solve the transient heat transfer problem
%
% 5. Compute the evolution of the temperature at a given Cartesian location
%
% 6. Define the expected solution
%
% 7. Verify the results
%
%% Function main body

%% 0. Read input

% Absolute tolerances
absTol = 1e-12;
absTol3 = 1e-12*1e3;

%% 1. Parse data from GiD input file
pathToCase = '../../inputGiD/FEMThermalConductionAnalysis/';
caseName = 'unitTest_transientWallHeating';
[thermalMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    propParameters, propNLinearAnalysis, propThermalDynamics, propGaussInt] = ...
    parse_ThermalModelFromGid(pathToCase, caseName, '');

%% 2. UI

% On the computation of the body forces
computeBodyForces = 'undefined';

% Equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% On the writing the output function
propVTK.isOutput = false;
propVTK.VTKResultFile = 'undefined';

% On transient inhomogeneous Dirichlet boundary conditions
updateInhomDOFs = 'undefined';
propIDBC = [];

% On the postprocessing
propPostproc.computeAnalytical = 'undefined';

% Choose the appropriate matrix update computation corresponding to the
% chosen time integration scheme
if strcmp(propThermalDynamics.method,'IMPLICIT_EULER')
    propThermalDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsImplicitEulerThermalConduction;
    propThermalDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propThermalDynamics.method,'GALERKIN')
    propThermalDynamics.computeProblemMtrcsTransient =  ...
        @computeProblemMtrcsGalerkinThermalConduction;
    propThermalDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
elseif strcmp(propThermalDynamics.method,'CRANK_NICOLSON')
    propThermalDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsCrankNicolsonThermalConduction;
    propThermalDynamics.computeUpdatedVct = ...
        @computeBETITransientUpdatedVctAccelerationField;
else
    error('Invalid time integration method selected in propStrDynamics.method as %s',propThermalDynamics.method);
end

% Define the initial condition function 
computeInitialConditions = @computeInitCndsFEMThermalConductionAnalysis;

%% 3. Define intial temperature and applied heat flux
propThermalDynamics.temperatureInit = 200;
propNBC.flux = 0;

%% 4. Solve the transient heat transfer problem
[THistory, ~, minElSize] = solve_FEMThermalConductionTransient ...
    (thermalMsh, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    updateInhomDOFs, propNBC, @computeLoadVctFEMThermalConductionAnalysis, ...
    propParameters, computeBodyForces, propAnalysis, computeInitialConditions, ...
    @computeStiffMtxAndLoadVctFEMThermalConductionAnalysisCST, ...
    propNLinearAnalysis, propIDBC, propThermalDynamics, solve_LinearSystem, ...
    @solve_FEMLinearSystem, propGaussInt, propVTK, caseName,'');

%% 5. Compute the evolution of the temperature at a given Cartesian location
x = 0.6;
y = 0.6;
[~, resultantNumerical, ~] = ...
    computeTemperatureAtPointOverTime ...
    (x, y, thermalMsh, THistory, ...
    propThermalDynamics, propPostproc, '');

%% 6. Define the expected solution

% Define the expected solution in terms of the temperature evolution at a
% given Cartesian location
expSolResultantNumerical = 1.0e+02*[   2.000000000000000
                                       2.384579833841253
                                       3.183390813663897
                                       3.622456346352308
                                       3.807643851271036
                                       4.044361657517795
                                       4.189430065024848
                                       4.337857164014472
                                       4.446967514532731
                                       4.542487135477225
                                       4.621857627119399
                                       4.684453857367032
                                       4.741059788293426
                                       4.782640465585473
                                       4.822496503940799
                                       4.850415356443687
                                       4.878217290974058
                                       4.897138793656419
                                       4.916381766137766
                                       4.929321589953344
                                       4.942540815468968
                                       4.951474498676134
                                       4.960481606321835
                                       4.966714980153099
                                       4.972792712138925
                                       4.977194310928822
                                       4.981245491452873
                                       4.984395631691392
                                       4.987052964916552
                                       4.989340847846956
                                       4.991046190196660];

% Define the expected solution in terms of the minimum element area size
expSolMinElSize = 0.062500000000000;

%% 7. Verify the results
testCase.verifyEqual(resultantNumerical, expSolResultantNumerical,'AbsTol', absTol3);
testCase.verifyEqual(minElSize,expSolMinElSize,'AbsTol',absTol);

end