function testTransientSquareCavity(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the transient thermal conduction analysis over the cavity benchmark
%
% Function layout :
%
% 0. Read input
%
% 1. Parse data from GiD input file
%
% 2. UI
%
% 3. Define intial temperature and/or applied heat flux
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
absTol = 1e-15;
absTol2 = 1e-15*1e2;

%% 1. Parse data from GiD input file
pathToCase = '../../inputGiD/FEMThermalConductionAnalysis/';
caseName = 'unitTest_transientSquareCavity';
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

%% 3. Define intial temperature and/or applied heat flux
propThermalDynamics.temperatureInit = 300;
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
expSolResultantNumerical = 1.0e+02*[   3.000000000000000
                                       3.019944914972513
                                       3.058518686695047
                                       3.103137320020808
                                       3.144943734234173
                                       3.180600894781213
                                       3.209781822733505
                                       3.233265094170095
                                       3.252062550864475
                                       3.267107487460751
                                       3.279174400515569
                                       3.288880547147621
                                       3.296710718612609
                                       3.303044515601801
                                       3.308179868359673
                                       3.312351696404856
                                       3.315746232657875
                                       3.318511896471344
                                       3.320767542454492
                                       3.322608749176758
                                       3.324112649677426
                                       3.325341673456490
                                       3.326346469812379
                                       3.327168209645261
                                       3.327840410594292
                                       3.328390392976242
                                       3.328840447127539
                                       3.329208773301406
                                       3.329510241039584
                                       3.329757004397274
                                       3.329959001493815];

% Define the expected solution in terms of the minimum element area size
expSolMinElSize = 0.062500000000000;

%% 7. Verify the results
testCase.verifyEqual(resultantNumerical, expSolResultantNumerical,'AbsTol', absTol2);
testCase.verifyEqual(minElSize,expSolMinElSize,'AbsTol',absTol);

end