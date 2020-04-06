function testIGA4TransientTaylorGreenVortices2D(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the stabilized isogeometric analysis for the 2D transient Stokes
% and Navier-Stokes equations for the benchmark problem of the Taylor-Green
% vortices. For this problem there is an analytical solution as a function
% of time and space.
%
% Function layout :
%
% 0. Read input
%
% 1. Define NURBS parameters
%
% 2. Define the flow parameters
%
% 3. UI
%
% 4. Perform h- and p-refinement
%
% 5. Define the Dirichlet boundary conditions
%
% 6. Define the Neumann boundary conditions
%
% 7. Create computational information for the patch
%
% 8. Define the transient analysis parameters
%
% 9. Define the nonlinear analysis parameters
%
% 10. Solve the the transient Stokes Problem using the Boosak scheme
%
% 11. Solve the transient Navier-Stokes Problem using the Bossak scheme
%
% 12. Define the expected solutions
%
% 13. Verify the results
%
%% Function main body

%% 0. Read input

% Define absolute tolerances
absTol = 1e-15;

%% 1. Define NURBS parameters

% Geometrical parameters
channelLength = 2*pi;
channelHeight = 2*pi;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

factor = 1/2;

% Control Point coordinates

% x-coordinate
CP(:,:,1) = [-factor*channelLength -factor*channelLength
             factor*channelLength factor*channelLength];
        
% y-coordinate
CP(:,:,2) = [-factor*channelHeight factor*channelHeight
             -factor*channelHeight factor*channelHeight];
        
% z-coordinate
CP(:,:,3) = [0 0
             0 0];
         
% Control Point weights
CP(:,:,4) = [1 1
             1 1];
         
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = false;
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));
for i = 1:nxi
    for j = 1:neta
        if CP(i, j, 4) ~= 1
            isNURBS = true;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% 2. Define the flow parameters

% Kinematic viscosity
parameters.nue = 1e-3;

%% 3. UI

% Analysis type
analysis.type = 'isogeometricIncompressibleFlowAnalysis';

% Function handle to the linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Define the body-force vector
amplification = 0;
computeBodyForces = amplification*[1 0]';

% Define the function handle for the initial conditions
computeInitCnds = @computeInitialConditionsForTaylorGreenVorticesIGA4NSE2D;

% Integration parameters

% type: 'default', 'user'
propInt.type = 'default';
if strcmp(propInt.type,'user')
    propInt.xiNGP = 6;
    propInt.etaNGP = 3;
    propInt.xiNGPForLoad = 6;
    propInt.etaNGPForLoad = 3;
    propInt.nGPForLoad = 6;
end

%% 4. Perform h- and p-refinement

% Degree elevation
tp = 0;  
tq = 0;
[Xi, Eta, CP, p, q] = degreeElevateBSplineSurface ...
    (p, q, Xi, Eta, CP, tp, tq, '');

% Knot insertion
xiRef = 4;
etaRef = 4;
[Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, xiRef, etaRef, '');

%% 5. Define the Dirichlet boundary conditions

% Homogeneous Dirichlet Boundary Conditions
homDOFs = [];

% Inhomogeneous Dirichlet Boundary Conditions
inhomDOFs = [];
propIDBC.numCnd = 12;
propIDBC.xiExtension = zeros(propIDBC.numCnd,2);
propIDBC.etaExtension = zeros(propIDBC.numCnd,2);
propIDBC.prescribedDirection = zeros(propIDBC.numCnd,1);
propIDBC.isUniqueOnBoundary = zeros(propIDBC.numCnd,1);

% Iterate over all the boundary conditions and assign their values
for i = 1:4
    % Find the correct span for the application of the inhomogeneous DBC
    if i == 1
        xiExtension = [0 1];
        etaExtension = [0 0];
    elseif i == 2
        xiExtension = [0 1];
        etaExtension = [1 1];
    elseif i == 3
        xiExtension = [0 0];
        etaExtension = [0 1];
    elseif i == 4
        xiExtension = [1 1];
        etaExtension = [0 1]; 
    end
        
    % x-component of the velocity field
    propIDBC.xiExtension(3*i-2,:) = xiExtension;
    propIDBC.etaExtension(3*i-2,:) = etaExtension;
    propIDBC.prescribedDirection(3*i-2) = 1;
    propIDBC.isUniqueOnBoundary(3*i-2) = true;
    propIDBC.irb(3*i-2,:) = ...
        findDofs3D(inhomDOFs,propIDBC.xiExtension(3*i-2,:), ...
        propIDBC.etaExtension(3*i-2,:), propIDBC.prescribedDirection(3*i-2),CP);
    
    % y-component of the velocity field
    propIDBC.xiExtension(3*i-1,:) = xiExtension;
    propIDBC.etaExtension(3*i-1,:) = etaExtension;
    propIDBC.prescribedDirection(3*i-1) = 2;
    propIDBC.isUniqueOnBoundary(3*i-1) = true;
    propIDBC.irb(3*i-1,:) = ...
        findDofs3D(inhomDOFs,propIDBC.xiExtension(3*i-1,:), ...
        propIDBC.etaExtension(3*i-1,:), propIDBC.prescribedDirection(3*i-1),CP);
    
    % pressure field
    propIDBC.xiExtension(3*i,:) = xiExtension;
    propIDBC.etaExtension(3*i,:) = etaExtension;
    propIDBC.prescribedDirection(3*i) = 3;
    propIDBC.isUniqueOnBoundary(3*i) = true;
    propIDBC.irb(3*i,:) = ...
        findDofs3D(inhomDOFs,propIDBC.xiExtension(3*i,:), ...
        propIDBC.etaExtension(3*i,:), propIDBC.prescribedDirection(3*i),CP);
end

% The prescribed values by function pointers
propIDBC.prescribedValue = ...
    {@computeXVelocityComponentForTaylorGreenVortices2D,...
    @computeYVelocityComponentForTaylorGreenVortices2D,...
    @computePressureFieldForTaylorGreenVortices2D,...
    @computeXVelocityComponentForTaylorGreenVortices2D,...
    @computeYVelocityComponentForTaylorGreenVortices2D,...
    @computePressureFieldForTaylorGreenVortices2D,...
    @computeXVelocityComponentForTaylorGreenVortices2D,...
    @computeYVelocityComponentForTaylorGreenVortices2D,...
    @computePressureFieldForTaylorGreenVortices2D,...
    @computeXVelocityComponentForTaylorGreenVortices2D,...
    @computeYVelocityComponentForTaylorGreenVortices2D,...
    @computePressureFieldForTaylorGreenVortices2D};

% Flag on the dominance of the inhomogeneous bc's to the homogeneous
propIDBC.isDominant = false;

% Find the DOFs where inhomogeneous Dirichlet boundary conditions are
% applied
for i = 1:propIDBC.numCnd
    inhomDOFs = mergesorted(inhomDOFs,propIDBC.irb(i,:));
end
inhomDOFs = unique(inhomDOFs);

%% 6. Define the Neumann boundary conditions
propNBC.noCnd = 2;
propNBC.xiLoadExtension = {[1 1] [1 1]};
propNBC.etaLoadExtension = {[0 1] [0 1]};
propNBC.loadAmplitude = {0 0};
propNBC.loadDirection = {1 2};
propNBC.isFollower = [false
                      false];
propNBC.computeLoadVct = ...
    {'computeLoadVctLineIGAIncompressibleNavierStokesFlow' ...
     'computeLoadVctLineIGAIncompressibleNavierStokesFlow'};

%% 7. Create computational information for the patch
BSplinePatch = fillUpPatch ...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, ...
    inhomDOFs, [], [], [], propNBC, [], [], [], [], [], propInt);

%% 8. Define the transient analysis parameters

% Select time integration scheme
% method: 'explicitEuler', 'Bossak'
propFldDynamics.method = 'Bossak';

% Computation of transient problem matrices
propFldDynamics.timeDependence = 'transient';
propFldDynamics.computeProblemMtrcsTransient = @computeProblemMtrcsBossakIGA4NSE;
propFldDynamics.computeUpdatedVct = @computeBossakTIUpdatedVctAccelerationFieldIGA4NSE;

% Parameters of the selected scheme
% alphaBeta <= .5 (for unconditional stability)
propFldDynamics.alphaBeta = -.1;

% alphaBeta + gamma >= .25 (for unconditional stability)
propFldDynamics.gamma = .5 - propFldDynamics.alphaBeta;

% The start and the end time of the simulation
propFldDynamics.TStart = 0;
propFldDynamics.TEnd = 1;

% The number of time steps
propFldDynamics.noTimeSteps = 2; % 1e1, 1e3

% The time step
propFldDynamics.dt = ...
    (propFldDynamics.TEnd - propFldDynamics.TStart)/propFldDynamics.noTimeSteps;
propFldDynamics.isAdaptive = true;

%% 9. Define the nonlinear analysis parameters
propNLinearAnalysis.method = 'Newton';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 1e-9;
propNLinearAnalysis.maxIter = 50;

%% 10. Solve the the transient Stokes Problem using the Boosak scheme
[upHistoryStokes, ~, ~] = solve_IGATransientFlow ...
    (analysis, BSplinePatch, computeBodyForces, computeInitCnds, ...
    @computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4StokesE2D, ...
    solve_LinearSystem, @solve_IGALinearSystem, propIDBC, ...
    propFldDynamics, propNLinearAnalysis,'');

%% 11. Solve the transient Navier-Stokes Problem using the Bossak scheme
[upHistoryNavierStokes, resHistNavierStokes, minElSize] = solve_IGATransientFlow ...
    (analysis, BSplinePatch, computeBodyForces, computeInitCnds, ...
    @computeIGAVMSStabMtxAndVct4BossakTINewtonNLinear4NSE2D, ...
    solve_LinearSystem, @solve_IGANLinearSystem, propIDBC, ...
    propFldDynamics, propNLinearAnalysis,'');

%% 12. Define the expected solutions

% Expected solution in terms of the solution history using Stokes
expUpHistoryStokes = [-0.000000000000000  -0.000000000000000  -0.000000000000000
                       0.000000000000000   0.000000000000000   0.000000000000000
                      -0.500000000000000  -0.499000999333667  -0.498003994671996
                       0.000000000000000   0.000000000000000   0.000000000000000
                       1.000000000000000   0.999000499833375   0.998001998667333
                                       0                   0                   0
                       0.000000000000000   0.000000000000000   0.000000000000000
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                      -0.500000000000000  -0.499000999333667  -0.498003994671996
                       0.000000000000000   0.000000000000000   0.000000000000000
                      -1.000000000000000  -0.999000499833375  -0.998001998667333
                                       0                   0                   0
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                      -0.500000000000000  -0.499000999333667  -0.498003994671996
                      -1.000000000000000  -0.999000499833375  -0.998001998667333
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                                       0                   0                   0
                       0.000000000000000   0.284154373797529   0.070496233125829
                      -0.000000000000000  -0.284264193674137  -0.070554077267265
                       0.500000000000000   0.500541735678236   0.501094382637720
                       1.000000000000000  -0.140422416214688   0.719809817641996
                       0.000000000000000   0.000065875886402   0.000034704626120
                                       0   0.000422363997633   0.000963598000606
                       0.000000000000000   0.284264193674137   0.070554077267265
                       0.000000000000000   0.284154373797529   0.070496233125829
                       0.500000000000000   0.500541735678236   0.501094382637721
                      -1.000000000000000  -0.999000499833375  -0.998001998667333
                       0.000000000000000   0.000000000000000   0.000000000000000
                                       0                   0                   0
                       0.000000000000000   0.000000000000000   0.000000000000000
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                      -0.500000000000000  -0.499000999333667  -0.498003994671996
                      -0.000000000000000   0.000065875886402   0.000034704626120
                      -1.000000000000000   0.140422416214688  -0.719809817641996
                                       0   0.000422363997634   0.000963598000607
                      -0.000000000000000   0.000000000000000  -0.000000000000000
                       0.000000000000000   0.000000000000000  -0.000000000000000
                      -0.500000000000000  -0.499507844568381  -0.499160150883463
                      -0.000000000000000  -0.000065875886402  -0.000034704626120
                       1.000000000000000  -0.140422416214688   0.719809817641995
                                       0   0.000422363997634   0.000963598000607
                       0.000000000000000   0.000000000000000   0.000000000000000
                       0.000000000000000   0.000000000000000   0.000000000000000
                      -0.500000000000000  -0.499000999333667  -0.498003994671996
                       1.000000000000000   0.999000499833375   0.998001998667333
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                                       0                   0                   0
                      -0.000000000000000  -0.284264193674137  -0.070554077267265
                      -0.000000000000000  -0.284154373797529  -0.070496233125829
                       0.500000000000000   0.500541735678235   0.501094382637720
                      -1.000000000000000   0.140422416214687  -0.719809817641995
                       0.000000000000000  -0.000065875886402  -0.000034704626120
                                       0   0.000422363997633   0.000963598000607
                      -0.000000000000000  -0.284154373797529  -0.070496233125829
                       0.000000000000000   0.284264193674137   0.070554077267265
                       0.500000000000000   0.500541735678236   0.501094382637722
                       1.000000000000000   0.999000499833375   0.998001998667333
                       0.000000000000000   0.000000000000000   0.000000000000000
                                       0                   0                   0
                       0.000000000000000   0.000000000000000   0.000000000000000
                       0.000000000000000   0.000000000000000   0.000000000000000
                      -0.500000000000000  -0.499000999333667  -0.498003994671996
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                       1.000000000000000   0.999000499833375   0.998001998667333
                                       0                   0                   0
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                      -0.500000000000000  -0.499000999333667  -0.498003994671996
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                      -1.000000000000000  -0.999000499833375  -0.998001998667333
                                       0                   0                   0
                       0.000000000000000   0.000000000000000   0.000000000000000
                      -0.000000000000000  -0.000000000000000  -0.000000000000000
                      -0.500000000000000  -0.499000999333667  -0.498003994671996];
                  
% Expected solution in terms of the solution history using Navier-Stokes
expUpHistoryNavierStokes = [  -0.000000000000000  -0.000000000000000  -0.000000000000000
                               0.000000000000000   0.000000000000000   0.000000000000000
                              -0.500000000000000  -0.499000999333667  -0.498003994671996
                               0.000000000000000   0.000000000000000   0.000000000000000
                               1.000000000000000   0.999000499833375   0.998001998667333
                                               0                   0                   0
                               0.000000000000000   0.000000000000000   0.000000000000000
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                              -0.500000000000000  -0.499000999333667  -0.498003994671996
                               0.000000000000000   0.000000000000000   0.000000000000000
                              -1.000000000000000  -0.999000499833375  -0.998001998667333
                                               0                   0                   0
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                              -0.500000000000000  -0.499000999333667  -0.498003994671996
                              -1.000000000000000  -0.999000499833375  -0.998001998667333
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                                               0                   0                   0
                               0.000000000000000   0.017378251030926   0.034965912112617
                              -0.000000000000000   0.025886589370529   0.067685343779850
                               0.500000000000000   0.162336994613165   0.173989332457320
                               1.000000000000000   0.999973803398018   1.003180364124872
                               0.000000000000000  -0.024062765096446  -0.058135386777936
                                               0  -0.168701666055728  -0.104229984750095
                               0.000000000000000  -0.025886589370529  -0.067685343779850
                               0.000000000000000   0.017378251030927   0.034965912112617
                               0.500000000000000   0.162336994613164   0.173989332457321
                              -1.000000000000000  -0.999000499833375  -0.998001998667333
                               0.000000000000000   0.000000000000000   0.000000000000000
                                               0                   0                   0
                               0.000000000000000   0.000000000000000   0.000000000000000
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                              -0.500000000000000  -0.499000999333667  -0.498003994671996
                              -0.000000000000000  -0.024062765096446  -0.058135386777935
                              -1.000000000000000  -0.999973803398018  -1.003180364124873
                                               0  -0.168701666055728  -0.104229984750096
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                               0.000000000000000   0.000000000000000   0.000000000000001
                              -0.500000000000000  -0.304228916667252  -0.345468715723775
                              -0.000000000000000   0.024062765096446   0.058135386777935
                               1.000000000000000   0.999973803398017   1.003180364124871
                                               0  -0.168701666055728  -0.104229984750095
                               0.000000000000000   0.000000000000000   0.000000000000000
                               0.000000000000000   0.000000000000000   0.000000000000000
                              -0.500000000000000  -0.499000999333667  -0.498003994671996
                               1.000000000000000   0.999000499833375   0.998001998667333
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                                               0                   0                   0
                              -0.000000000000000   0.025886589370529   0.067685343779850
                              -0.000000000000000  -0.017378251030926  -0.034965912112617
                               0.500000000000000   0.162336994613165   0.173989332457320
                              -1.000000000000000  -0.999973803398017  -1.003180364124872
                               0.000000000000000   0.024062765096446   0.058135386777935
                                               0  -0.168701666055728  -0.104229984750095
                              -0.000000000000000  -0.017378251030926  -0.034965912112617
                               0.000000000000000  -0.025886589370529  -0.067685343779850
                               0.500000000000000   0.162336994613164   0.173989332457320
                               1.000000000000000   0.999000499833375   0.998001998667333
                               0.000000000000000   0.000000000000000   0.000000000000000
                                               0                   0                   0
                               0.000000000000000   0.000000000000000   0.000000000000000
                               0.000000000000000   0.000000000000000   0.000000000000000
                              -0.500000000000000  -0.499000999333667  -0.498003994671996
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                               1.000000000000000   0.999000499833375   0.998001998667333
                                               0                   0                   0
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                              -0.500000000000000  -0.499000999333667  -0.498003994671996
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                              -1.000000000000000  -0.999000499833375  -0.998001998667333
                                               0                   0                   0
                               0.000000000000000   0.000000000000000   0.000000000000000
                              -0.000000000000000  -0.000000000000000  -0.000000000000000
                              -0.500000000000000  -0.499000999333667  -0.498003994671996];

% Expected solution in terms of the residual history using Navier-Stokes
expResHistNavierStokes = [     0.445656451269541   0.378999245425081
                               2.816551117756964   2.803527200414355
                               0.056709445089952   0.057061569432210
                               0.001083698325197   0.001104427805278
                               0.000016370640691   0.000019240846167
                               0.000000248869656   0.000000350049663
                               0.000000003853803   0.000000006387773
                               0.000000000059977   0.000000000116614
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0
                                               0                   0];

% Expected solution in terms of the minimum element area size
expMinElSize = 2.467401100272339;

%% 13. Verify the results
testCase.verifyEqual(upHistoryStokes, expUpHistoryStokes, 'AbsTol', absTol);
testCase.verifyEqual(upHistoryNavierStokes, expUpHistoryNavierStokes, 'AbsTol', absTol);
testCase.verifyEqual(resHistNavierStokes, expResHistNavierStokes, 'AbsTol', absTol);
testCase.verifyEqual(minElSize, expMinElSize, 'AbsTol', absTol);

end