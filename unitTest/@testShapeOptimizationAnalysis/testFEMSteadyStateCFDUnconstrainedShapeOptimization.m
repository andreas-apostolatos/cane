function testFEMSteadyStateCFDUnconstrainedShapeOptimization(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the unconstrained shape optimization algorithm using the steepest
% descent method on the case of a 2D cylinder in flow with object function
% the drag force on the cylinder and one design variable which is the
% radius of the cylinder. The expected solution is the radius to converge
% to zero.
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
% 4. Get the center of the cylinder in flow
%
% 5. Initializations
%
% 6. Loop over all the optimization iterations
% ->
%    6i. Update the mesh for the nominal state according to the updated design
%
%   6ii. Solve the nominal steady-state CFD problem
%
%  6iii. Compute the nominal drag force
%
%   6iv. Check if the objective function is below a threshold
%
%    6v. Move the mesh according to the prescribed perturbation
%
%   6vi. Solve the perturbed steady-state CFD problem
%
%  6vii. Compute the perturbed drag force
%
% 6viii. Compute sensitivity via finite differencing
%
%   6ix. Compute Hessian of the system
%
%    6x. Compute the design update5
%
%   6xi. Increment iteration counter
% <-
%
% 7. Define the expected solutions
%
% 8. Verify the results
%
%% Function main body

%% 0. Read input

% Define absolute tolerance
absTol = 1e-15;

% Define the path to the case
pathToCase = '../../inputGiD/FEMComputationalFluidDynamicsAnalysis/';
caseName = 'unitTest_cylinder2D_CFD_unconstrainedOptimization';

%% 1. Parse the data from the GiD input file
[fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propALE, ~, propAnalysis, ...
    propParameters, propNLinearAnalysis, propFldDynamics, propGaussInt, ...
    propPostproc] = parse_FluidModelFromGid(pathToCase, caseName, '');

%% 2. UI

% On the body forces
computeBodyForces = @computeConstantVerticalBodyForceVct;

% On the transient analysis properties
if strcmp(propFldDynamics.method,'bossak')
    propFldDynamics.computeProblemMtrcsTransient = ...
        @computeProblemMtrcsBossakFEM4NSE;
    propFldDynamics.computeUpdatedVct = ...
        @computeBossakTIUpdatedVctAccelerationFieldFEM4NSE2D;
end

% Free DOFs of the system (actual DOFs over which the solution is computed)
numNodes = length(fldMsh.nodes(:, 1));
numDOFs = propAnalysis.noFields*numNodes;
freeDOFs = 1:numDOFs;
prescribedDOFs = mergesorted(homDOFs, inhomDOFs);
prescribedDOFs = unique(prescribedDOFs);
freeDOFs(ismember(freeDOFs, prescribedDOFs)) = [];

% Define the name of the vtk file from where to resume the simulation
propOutput.isOutput = false;
propOutput.writeOutputToFile = 'undefined';

% Maximum optimization iterations and drag tolerance
maxIter = 10;
tolDrag = 1e-4;

% Perturbation size and design update scaling
epsilonTilde = 1e-2;
alphaTilde = 1e-2;

% Function handle for the computation of the Hessian
djdp1_ = @(djdp1_, p1) djdp1_;

% Dummy parameters
nodesSaved = 'undefined';
uMeshALE = 'undefined';

%% 3. Choose the equation system solver
if strcmp(propAnalysis.type,'NAVIER_STOKES_2D')
    solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
elseif strcmp(propAnalysis.type,'NAVIER_STOKES_3D')
    solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;
else
    error('Neither NAVIER_STOKES_2D or NAVIER_STOKES_3D has been chosen');
end

%% 4. Get the center of the cylinder in flow
propALE.propUser.x0 = mean(fldMsh.nodes(propPostproc.nodesDomain{1}(:, 1), 1));
propALE.propUser.y0 = mean(fldMsh.nodes(propPostproc.nodesDomain{1}(:, 1), 2));
radiusInit = max(fldMsh.nodes(propPostproc.nodesDomain{1}(:, 1), 1)) - propALE.propUser.x0;

%% 5. Initializations

% Initialize the solution of the state equations
up = computeNullInitialConditionsFEM4NSE ...
    (propAnalysis, fldMsh, 'undefined', 'undefined', 'undefined', ... 
    'undefined', 'undefined', 'undefined');

% Initialize objective and design history
optHistory = zeros(maxIter, 2);

% Initialize the minimum element edge size history
minElSizeHistory = zeros(maxIter, 1);

% Initialize the design update
dr = 0;

% Initialize the design
radius = radiusInit + dr;

% Initialize the optimization counter
counterOpt = 1;

%% 6. Loop over all the optimization iterations
while counterOpt <= maxIter
    %% 6i. Update the mesh for the nominal state according to the updated design
    propALE.propUser.dr = dr;
    radius = radius + dr;
    [fldMsh, ~, ~, ~, ~, ~] = computeUpdatedMeshAndVelocitiesPseudoStrALE2D ...
        (fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, freeDOFs, nodesSaved, ...
        propALE, solve_LinearSystem, propFldDynamics, counterOpt);
    fldMsh.initialNodes = fldMsh.nodes;
    
    %% 6ii. Solve the nominal steady-state CFD problem
    [up, FComplete, ~, minElSizeHistory(counterOpt, 1)] = ...
        solve_FEMVMSStabSteadyStateNSE ...
        (fldMsh, up, homDOFs, inhomDOFs, valuesInhomDOFs, uMeshALE, ...
        propParameters, computeBodyForces, propAnalysis, ...
        solve_LinearSystem, propFldDynamics, propNLinearAnalysis, ...
        counterOpt, propGaussInt, propOutput, caseName, '');
    
    %% 6iii. Compute the nominal drag force
    postProc_update = computePostProc ...
        (FComplete, propAnalysis, propParameters, propPostproc);
    forcesOnDomain = postProc_update.valuePostProc{1};
    drag_nom = forcesOnDomain(1, 1);
    
    %% 6iv. Check if the objective function is below a threshold
    optHistory(counterOpt, 1) = drag_nom;
    optHistory(counterOpt, 2) = radius;
    if abs(drag_nom) < tolDrag
        break;
    end
 
    %% 6v. Move the mesh according to the prescribed perturbation
    propALE.propUser.dr = epsilonTilde;
    [fldMsh_p1, ~, ~, ~, ~, ~] = ...
        computeUpdatedMeshAndVelocitiesPseudoStrALE2D ...
        (fldMsh, homDOFs, inhomDOFs, valuesInhomDOFs, freeDOFs, nodesSaved, ....
        propALE, solve_LinearSystem, propFldDynamics, counterOpt);
    fldMsh_p1.initialNodes = fldMsh_p1.nodes;
    
    %% 6vi. Solve the perturbed steady-state CFD problem
    [~, FComplete, ~, ~] = solve_FEMVMSStabSteadyStateNSE ...
        (fldMsh_p1, up, homDOFs, inhomDOFs, valuesInhomDOFs, uMeshALE, ...
        propParameters, computeBodyForces, propAnalysis, ...
        solve_LinearSystem, propFldDynamics, propNLinearAnalysis, ...
        counterOpt, propGaussInt, propOutput, caseName, '');
    
    %% 6vii. Compute the perturbed drag force
    postProc_update = computePostProc ...
        (FComplete, propAnalysis, propParameters, propPostproc);
    forcesOnDomain = postProc_update.valuePostProc{1};
    drag = forcesOnDomain(1, 1);

    %% 6viii. Compute sensitivity via finite differencing
    drag_dp1 = (drag - drag_nom)/epsilonTilde;
    
    %% 6ix. Compute Hessian of the system
    djd1 = djdp1_(drag_dp1);
    
    %% 6x. Compute the design update
    dr = -alphaTilde*djd1;

    %% 6xi. Increment iteration counter
    counterOpt = counterOpt + 1;
end

%% 7. Define the expected solutions

% Define the expected solution in terms of the drag and radius history
expSolOptHistory = [   0.009744743326393   0.050000000000000
                       0.009066984758399   0.047193504458965
                       0.008476041739562   0.044566754172162
                       0.007954501804678   0.042095712459956
                       0.007493423224948   0.039775482244403
                       0.007082681925176   0.037588870546007
                       0.006713976333818   0.035520761373707
                       0.006380851272066   0.033558144886823
                       0.006077626927841   0.031690435787415
                       0.005800091172037   0.029907505557754];
                   
% Define the expected solution in terms of the minimum element edge size
expSolMinElEdgeSizeHistory = [ 0.015641201999846
                               0.014763258266501
                               0.013941543455398
                               0.013168538277109
                               0.012442710800144
                               0.011758682619875
                               0.011111725066889
                               0.010497768329589
                               0.009913501055488
                               0.009355754784654];

%% 8. Verify the results
testCase.verifyEqual(optHistory, expSolOptHistory, 'AbsTol', absTol);
testCase.verifyEqual(minElSizeHistory, expSolMinElEdgeSizeHistory, 'AbsTol', absTol);

end