function testIGA4StokesSteadyState2D(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the stabilized isogeometric analysis for the 2D steady-state
% Stokes equations for the flow in unit square domain.
%
% Function layout :
%
% 0. Read input
%
% 1. Define NURBS parameters
%
% 2. Define the material constants
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
% 8. Solve the steady-state Stokes problem
%
% 9. Compute the relative errors in the L2-norm
%
% 10. Define the expected solutions
%
% 11. Verify the results
%
%% Function main body

%% 0. Read input

% Define absolute tolerances
absTol = 1e-15;
absTol2 = absTol*1e2;

%% 1. Define NURBS parameters

% Geometrical parameters
channelLength = 1;
channelHeight = 1;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinate
CP(:,:,1) = [0             0
             channelLength channelLength];
        
% y-coordinate
CP(:,:,2) = [0 channelHeight
             0 channelHeight];
        
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
        if CP(i, j, 4)~=1
            isNURBS = true;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% 2. Define the material constants

% Kinematic viscosity
parameters.nue = 1e2;

%% 3. UI

% Analysis type
analysis.type = 'isogeometricIncompressibleFlowAnalysis';

% Function handle to the linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = @solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Function handle to the computation of the analytical body force vector
computeBodyForces = @bodyForcesForAnalyticalSolutionToStokesProblemInUnitSquare;

% Integration parameters

% type: 'default', 'user'
propInt.type = 'user';
if strcmp(propInt.type,'user')
    propInt.xiNGP = 6;
    propInt.etaNGP = 3;
    propInt.xiNGPForLoad = 6;
    propInt.etaNGPForLoad = 3;
    propInt.nGPForLoad = 6;
end
propIntError.type = 'user';
propIntError.xiNGP = 10;
propIntError.etaNGP = 10;

%% 4. Perform h- and p-refinement

% Degree elevation
tp = 0;  
tq = 0;
[Xi, Eta, CP, p, q] = degreeElevateBSplineSurface ...
    (p, q, Xi, Eta, CP, tp, tq, '');

% Knot insertion
xiRef = 5;
etaRef = 5;
[Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, xiRef, etaRef, '');

%% 5. Define the Dirichlet boundary conditions

% Homogeneous Dirichlet Boundary Conditions
homDOFs = [];

% No Slip condition at the lower wall of the Channel
xiSupp = [0 1];
etaSupp = [0 0];   
dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);
xiSupp = [0 1];   
etaSupp = [0 0];   
dirSupp = 2;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% We have to constraint the pressure at one location so that we ensure
% uniqueness to the pressure space due to the ex
    % finite element analysis for the 2D and 3D Navier-Stokes equations.istence of the grad.
xiSupp = [0 0];   etaSupp = [0 0];   dirSupp = 3;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% No Slip condition at the upper wall of the Channel
xiSupp = [0 1];   
etaSupp = [1 1];   
dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSupp,etaSupp, dirSupp, CP);
xiSupp = [0 1];   
etaSupp = [1 1];   
dirSupp = 2;
homDOFs = findDofs3D(homDOFs, xiSupp, etaSupp, dirSupp, CP);

% No Slip condition at the left wall of the Channel
xiSupp = [0 0];   etaSupp = [0 1];   dirSupp = 1;
homDOFs = findDofs3D(homDOFs,xiSupp,etaSupp,dirSupp,CP);
xiSupp = [0 0];   etaSupp = [0 1];   dirSupp = 2;
homDOFs = findDofs3D(homDOFs,xiSupp,etaSupp,dirSupp,CP);

% No Slip condition at the right wall of the Channel
xiSupp = [1 1];   etaSupp = [0 1];   dirSupp = 1;
homDOFs = findDofs3D(homDOFs,xiSupp,etaSupp,dirSupp,CP);
xiSupp = [1 1];   etaSupp = [0 1];   dirSupp = 2;
homDOFs = findDofs3D(homDOFs,xiSupp,etaSupp,dirSupp,CP);

% Inhomogeneous Dirichlet Boundary Conditions
inhomDOFs = [];

propIDBC.numCnd = 0;
propIDBC.xiSpan = zeros(propIDBC.numCnd, 2);
propIDBC.etaSpan = zeros(propIDBC.numCnd, 2);
propIDBC.prescribedDirection = zeros(propIDBC.numCnd, 1);
propIDBC.isUniqueOnBoundary = zeros(propIDBC.numCnd, 1);

% The prescribed values by function pointers
% IDBC.prescribedValue = {@quadraticInletDistributionForVectorTransportProblems2D};
propIDBC.prescribedValue = {};

% Flag on the dominance of the inhomogeneous bc's to the homogeneous
propIDBC.isDominant = 1;

% Find the DOFs where inhomogeneous Dirichlet boundary conditions are
% applied
for i = 1:propIDBC.numCnd
    inhomDOFs = mergesorted(inhomDOFs,propIDBC.irb(i,:));
end
inhomDOFs = unique(inhomDOFs);

%% 6. Define the Neumann boundary conditions

% Transient Neumann boundary conditions

% Initialize the boundary conditions
propNBC.noCnd = 0; 
propNBC.xiSpan = zeros(propNBC.noCnd, 2);
propNBC.etaSpan = zeros(propNBC.noCnd, 2);
propNBC.loadAmplitude = zeros(propNBC.noCnd, 1);
propNBC.loadDirection = zeros(propNBC.noCnd, 1);

% Iterate over all the boundary conditions and assign their values
propNBC.xiSpan(1, :) = [1 1];
propNBC.etaSpan(1, :) = [0 1];
propNBC.loadAmplitude(1) = 1000;
propNBC.loadDirection(1) = 1;

% Assign the pointers to the load vector function computations
propNBC.loadVctComputation = {@computeLoadVctLineIGAIncompressibleNavierStokesFlow};

%% 7. Create computational information for the patch
BSplinePatch = fillUpPatch ...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, ...
    inhomDOFs, [], [], [], propNBC, [], [], [], [], [], propInt);

%% 8. Solve the steady-state Stokes problem
[up, ~, ~] = solve_IGAVMSStabSteadyStateStokesE2D...
    (analysis, BSplinePatch, computeBodyForces, solve_LinearSystem, ...
    propIDBC, propNBC, propInt, '');

%% 9. Compute the relative errors in the L2-norm
[errL2Velocity, errL2Pressure, minElArea] = ...
    computeIGAErrUnitDomainStokesE2D ...
    (p, Xi, q, Eta, CP, isNURBS, up, propIntError);

%% 10. Define the expected solutions

% Expected velocity error
expErrL2Velocity = 1.129823711239979e+02;

% Expected pressure error
expErrL2Pressure = 0.092524723211858;

% Expected minimum element area size
expMinElArea = 0.040000000000000;

%% 11. Verify the results
testCase.verifyEqual(errL2Velocity, expErrL2Velocity, 'AbsTol', absTol2);
testCase.verifyEqual(errL2Pressure, expErrL2Pressure, 'AbsTol', absTol2);
testCase.verifyEqual(minElArea, expMinElArea, 'AbsTol', absTol);

end