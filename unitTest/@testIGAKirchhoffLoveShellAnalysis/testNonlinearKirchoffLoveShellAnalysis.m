function testNonlinearKirchoffLoveShellAnalysis(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function definition
%
% Test the geometrically nonlinear Kirchoff-Love shell isogeometric 
% analysis over the single patch Scordelis-Lo roof benchmark example
%
% Function layout :
%
% 0. Read input
%
% 1. Define NURBS parameters
%
% 2. Define material constants
%
% 3. Define Gauss quadrature
%
% 4. Define h- and p-refinement
%
% 5. Define the expected solutions
%
% 6. Define Dirichlet and Neumann boundary conditions
%
% 7. Create the B-Spline patch array
%
% 8. Define nonlinear analysis parameters
%
% 9. Solve the nonlinear system
%
% 10. Verify the results
%
%% Function main body

%% 0. Read input

% Define a tolerance for the verification of the results
relaxationFactor = 1e2;
absTol = relaxationFactor*1e-15;
absTolRelaxed = relaxationFactor*absTol*1e5;
   
%% 1. Define NURBS parameters

% Global variables
Length = 50;
Radius = 25;

% Polynomial degrees
p = 1;
q = 2;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP(:, :, 1) = [-Length/2 -Length/2 -Length/2
               Length/2  Length/2  Length/2];
         
% y-coordinates
CP(:, :, 2) = [-Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)
               -Radius*sin(2*pi/9) 0 Radius*sin(2*pi/9)];
         
% z-coordinates
CP(:, :, 3) = [Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)
               Radius*cos(2*pi/9) Radius/cos(2*pi/9) Radius*cos(2*pi/9)];
       
% Weights
weight = cos(2*pi/9);
CP(:, :, 4) = [1 weight 1
               1 weight 1];

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

%% 2. Define material constants

% Young's modulus
parameters.E = 4.32e8;

% Poisson ratio
parameters.nue = .0;

% Thickness of the shell
parameters.t = .25;

% Density of the shell (used only for dynamics)
parameters.rho = 7850;

%% 3. Define Gauss quadrature

% Analysis type
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'user';
if strcmp(int.type,'user')
    int.xiNGP = p + 2;
    int.etaNGP = q + 1;
    int.xiNGPForLoad = ceil((p + 2)/2);
    int.etaNGPForLoad = ceil((q + 1)/2);
end

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

%% 4. Define h- and p-refinement

% Degree by which to elevate
tp = 1;
tq = 0;
[Xi, Eta, CP, p, q] = degreeElevateBSplineSurface ...
    (p, q, Xi, Eta, CP, tp, tq, '');

% Number of knots to exist in both directions
scaling = 1;
edgeRatio = ceil(Length/Radius/(sin(4*pi/9)));
refXi = edgeRatio*scaling;
refEta = ceil(4/3)*scaling;
[Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, refXi, refEta, '');

%% 5. Define the expected solutions
% Define reference solution on the displacements
expSolDisp = [     0
                   0
                   0
                   0.000307416611017
                   0.001828491247823
                  -0.017844084787584
                   0.005424875291000
                   0.001144679441999
                  -0.030425014925422
                   0.010542333970979
                   0.001828491247827
                  -0.017844084787585
                   0.010849750581998
                   0
                   0
                   0.008231504025116
                   0
                   0
                   0.008051548791940
                   0.000062932576136
                  -0.015351206806173
                   0.005424875290998
                  -0.000016108057216
                  -0.028613784217258
                   0.002798201790055
                   0.000062932576133
                  -0.015351206806163
                   0.002618246556879
                   0
                   0
                   0.008231504025113
                   0
                   0
                   0.008051548791942
                  -0.000062932576136
                  -0.015351206806166
                   0.005424875291000
                   0.000016108057215
                  -0.028613784217276
                   0.002798201790055
                  -0.000062932576133
                  -0.015351206806172
                   0.002618246556886
                   0
                   0
                   0.000000000000001
                   0
                   0
                   0.000307416611016
                  -0.001828491247826
                  -0.017844084787587
                   0.005424875290998
                  -0.001144679442027
                  -0.030425014925467
                   0.010542333970981
                  -0.001828491247822
                  -0.017844084787594
                   0.010849750582001
                   0
                   0];
               
% Define the expected solution in terms of the Control Point deformation
expSolCPHist = struct([]);
for iComp = 1:4
    expSolCPHist{1}(:, :, iComp, 1) = CP(:, :, iComp);
end
expSolCPHist{1}(:, :, 1, 2) = [-25.000000000000000 -24.991768495974885 -24.991768495974888 -25.000000000000000
                               -16.666359250055653 -16.658615117874728 -16.658615117874724 -16.666359250055653
                               0.005424875290998   0.005424875290995   0.005424875290997   0.005424875290995
                               16.677209000637642  16.669464868456725  16.669464868456725  16.677209000637646
                               25.010849750581997  25.002618246556878  25.002618246556885  25.010849750582000];
expSolCPHist{1}(:, :, 2, 2) = [ -16.069690242163482  -9.099255856655059   9.099255856655059  16.069690242163482
                                -16.067861750915657  -9.099192924078924   9.099192924078924  16.067861750915657
                                -16.068545562721482  -9.099271964712274   9.099271964712274  16.068545562721454
                                -16.067861750915654  -9.099192924078926   9.099192924078926  16.067861750915661
                                -16.069690242163482  -9.099255856655059   9.099255856655059  16.069690242163482];
expSolCPHist{1}(:, :, 3, 2) = [19.151111077974452  25.000000000000000  25.000000000000000  19.151111077974452
                               19.133266993186869  24.984648793193827  24.984648793193834  19.133266993186865
                               19.120686063049032  24.971386215782744  24.971386215782726  19.120686063048986
                               19.133266993186869  24.984648793193838  24.984648793193827  19.133266993186858
                               19.151111077974452  25.000000000000000  25.000000000000000  19.151111077974452];
expSolCPHist{1}(:, :, 4, 2) = [1.000000000000000   0.883022221559489   0.883022221559489   1.000000000000000
                               1.000000000000000   0.883022221559489   0.883022221559489   1.000000000000000
                               1.000000000000000   0.883022221559489   0.883022221559489   1.000000000000000
                               1.000000000000000   0.883022221559489   0.883022221559489   1.000000000000000
                               1.000000000000000   0.883022221559489   0.883022221559489   1.000000000000000];

% Residual history
expSolResHistory = 1.0e+04 *[
                               3.778629700973396
                               0.406254542836505
                               0.00000744496092483343
                               0.0000000000831335726027163
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0
                               0];
                     
% Define the reference solution on the minimum element area size
expSolMinElArea = 50;

%% 6. Define Dirichlet and Neumann boundary conditions

% supports (Dirichlet boundary conditions)

% back and front curved edges are a rigid diaphragm
homDOFs = [];
xiSup = [0 0];
etaSup = [0 1];
for dirSupp = [2 3]
    homDOFs = findDofs3D(homDOFs, xiSup, etaSup, dirSupp, CP);
end
xiSup = [1 1];
etaSup = [0 1];
for dirSupp = [2 3]
    homDOFs = findDofs3D(homDOFs, xiSup, etaSup, dirSupp, CP);
end

% Fix the back left corner of the shell to avoid rigid body motions
xiSup = [0 0];
etaSup = [0 0];
dirSupp = 1;
homDOFs = findDofs3D(homDOFs, xiSup, etaSup, dirSupp, CP);

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables = [];

% load (Neuman boundary conditions)
FAmp = - 9e1;
NBC.noCnd = 1;
xib = [0 1];   etab = [0 1];   dirForce = 'z';
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC.isFollower(1, 1) = false;
NBC.isTimeDependent(1, 1) = false;

%% 7. Create the B-Spline patch array
BSplinePatch = fillUpPatch ...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, weakDBC, cables, NBC, [], [], [], [], [], int);

%% 8. Define nonlinear analysis parameters
propNLinearAnalysis.scheme = 'newtonRaphson';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 10e-7;
propNLinearAnalysis.maxIter = 120;
noPatches = 1;
propNLinearAnalysis.conservativeLoad = true(noPatches, 1);

%% 9. Solve the nonlinear system
graph = [];
plot_IGANonlinear = 'undefined';
[dHatNLinear, CPHistory, resHistory, isConverged, ~, minElASize] = ...
    solve_IGAKirchhoffLoveShellNLinear ...
    (BSplinePatch, propNLinearAnalysis, solve_LinearSystem, ...
    plot_IGANonlinear, graph, '');

%% 10. Verify the results
testCase.verifyEqual(dHatNLinear, expSolDisp, 'AbsTol', absTol);
testCase.verifyEqual(CPHistory, expSolCPHist, 'AbsTol', absTol);
testCase.verifyEqual(resHistory, expSolResHistory, 'AbsTol', absTolRelaxed);
testCase.verifyEqual(isConverged, true, 'AbsTol', absTol);
testCase.verifyEqual(minElASize, expSolMinElArea, 'AbsTol', absTol);

end
