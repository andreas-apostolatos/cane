%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Performs Multilevel Monte-Carlo for the 2D Navier-Stokes equations
%        in 2D (to be implemented)
%
% Date : 16.12.2019
%
%% Function main body

%% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_testFEM4NavierStokesSteadyStateFlowAroundCylinder2D';

%% Parse the data from the GiD input file
[fldMsh,homDBC,inhomDBC,valuesInhomDBC,nodesALE,~,analysis,parameters,...
    propNLinearAnalysis,propFldDynamics,gaussInt,postProc] = ...
    parse_FluidModelFromGid...
    (pathToCase,caseName,'');

%% GUI
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

%% Choose the equation system solver
if strcmp(analysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(analysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
end

%% Define the name of the vtk file from where to resume the simulation
VTKResultFile = 'undefined';

%% Change the input velocity to match the reference paper - parabolic input

% max input velocity defined in the reference paper
Umax = 0.3;

% change the input velocity to have the parabolic distribution
valuesInhomDBCModified = computeInletVelocityParabolic_unitTest(fldMsh, inhomDBC, valuesInhomDBC, Umax);

%% Solve the CFD problem
[~,FComplete,hasConverged,~] = solve_FEMVMSStabSteadyStateNSE2D...
    (fldMsh,homDBC,inhomDBC,valuesInhomDBCModified,nodesALE,parameters,...
    computeBodyForces,analysis,computeInitialConditions,...
    VTKResultFile,solve_LinearSystem,propFldDynamics,propNLinearAnalysis,...
    gaussInt,caseName,'');

%% Calculate drag and lift force from the nodal forces
postProc = computePostProc(FComplete, analysis, parameters, postProc);

%% Calculate drag and lift coefficient based on drag and lift force

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

%% END OF THE SCRIPT