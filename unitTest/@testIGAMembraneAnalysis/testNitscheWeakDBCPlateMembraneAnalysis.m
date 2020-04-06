function testNitscheWeakDBCPlateMembraneAnalysis(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the Nitsche method for a plate discretized with one bilinear 
% element subject to tip horizontal load and for which the analytical
% solution of the linear and the nonlinear problem at the tip displacement 
% is equal to 1.0 and 0.8803391469128941469, respectively.
% 
% Function layout :
%
% 0. Read input
%
% 1. Define the multipatch geometry in terms of NURBS
%
% 2. Define the material parameters
%
% 3. GUI
%
% 4. Define the refinement
%
% 5. Define the boundary conditions
%
% 6. Create the patch cell array
%
% 7. Define the expected solutions for the linear problem
%
% 8. Define the expected solutions for the nonlinear problem
%
% 9. Define the linear analysis parameters
%
% 10. Define the linear analysis parameters
%
% 11. Solve the steady-state linear problem
%
% 12. Solve the steady-state nonlinear problem
%
% 13. Verify the solution
%
%% Function main body

%% 0. Read input

% Define absolute tolerances
absTolDisp = 1e-14;
absTolResHistory = 1e-15*1e3;

%% 1. Define the multipatch geometry in terms of NURBS

% Global variables:
Length = 10;
Height = 1;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [-Length/2 -Length/2
             Length/2  Length/2];
         
% y-coordinates
CP(:,:,2) = [-Height/2 Height/2
             -Height/2 Height/2];
         
% z-coordinates
CP(:,:,3) = [0 0
             0 0];
         
% Weights
CP(:,:,4) = [1 1
             1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = false;
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));
for i= 1:nxi
    for j=1:neta
        if CP(i, j, 4) ~= true
            isNURBS = 1;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% 2. Define the material parameters

% Young's modulus
parameters.E = 1e2;

% Poisson ratio
parameters.nue = 0.0;

% Thickness of the membrane
parameters.t = 1;

% Density of the membrane (used only for dynamics)
parameters.rho = 7810;

% Prestress for the membrane
parameters.prestress.voigtVector = [0
                                    0
                                    0];

%% 3. GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'default';
if strcmp(int.type, 'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.xetaNGPForLoad = 6;
end

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','curvature','force','moment','shearForce'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x', 'y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

%% 4. Define the refinement

% Degree by which to elevate ("p-refinement")
a = 1; % default value = 2
tp = a;
tq = a;
[Xi, Eta, CP, p, q] = degreeElevateBSplineSurface ...
    (p, q, Xi, Eta, CP, tp, tq, '');

% Number of knots to exist in both directions ("h-refinement")
scaling = 2;
refXi = scaling;
refEta = scaling;
[Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, refXi, refEta, '');

%% 5. Define the boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Homogeneous Dirichlet boundary conditions
homDOFs = [];
xiSup = [0 1];
etaSup = [0 1];
for dirSupp = [2 3]
    homDOFs = findDofs3D(homDOFs, xiSup, etaSup, dirSupp, CP);
end

% Weak homogeneous Dirichlet boundary conditions
weakDBC.noCnd = 1;
weakDBC.method = 'nitsche';
weakDBC.estimationStabilPrm = false;
% weakDBC.alpha = 0;
weakDBC.xiExtension = {[0 0]};
weakDBC.etaExtension = {[0 1]};
weakDBC.int.type = 'default';

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Embedded cables
cables.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xib = 1;
etab = [0 1];
dirForce = 'x';
NBC.noCnd = 1;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude{1} = 10;
NBC.loadDirection = {dirForce};
NBC.computeLoadVct{1} = 'computeLoadVctLineIGAThinStructure';
NBC.isFollower(1, 1) = false;
NBC.isTimeDependent(1, 1) = false;

%% 6. Create the patch cell array
BSplinePatch = fillUpPatch ...
    (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, inhomDOFs, ...
    valuesInhomDOFs, weakDBC, cables, NBC, [], [], [], [], [], int);

%% 7. Define the expected solutions for the linear problem

% Define the expected solution for the Penalty method in terms of the
% Control Point displacements
expSolDispLinear =      [ -0.000000000000000
                                           0
                                           0
                           0.250000000000005
                                           0
                                           0
                           0.750000000000011
                                           0
                                           0
                           1.000000000000012
                                           0
                                           0
                          -0.000000000000000
                                           0
                                           0
                           0.250000000000004
                                           0
                                           0
                           0.750000000000011
                                           0
                                           0
                           1.000000000000012
                                           0
                                           0
                          -0.000000000000000
                                           0
                                           0
                           0.250000000000005
                                           0
                                           0
                           0.750000000000011
                                           0
                                           0
                           1.000000000000013
                                           0
                                           0
                           0.000000000000000
                                           0
                                           0
                           0.250000000000005
                                           0
                                           0
                           0.750000000000011
                                           0
                                           0
                           1.000000000000013
                                           0
                                           0];
                               
% Define the expected solution in terms of the Control Point displacements
expSolCPHistoryLinear{1}(:, :, :, 1) = CP;
expSolCPHistoryLinear{1}(:, :, 1, 2) = [ -5.000000000000000  -5.000000000000000  -5.000000000000000  -5.000000000000000
                                         -2.249999999999996  -2.249999999999996  -2.249999999999995  -2.249999999999996
                                          3.250000000000012   3.250000000000012   3.250000000000011   3.250000000000012
                                          6.000000000000012   6.000000000000012   6.000000000000012   6.000000000000012];
expSolCPHistoryLinear{1}(:, :, 2, 2) = [ -0.500000000000000  -0.250000000000000   0.250000000000000   0.500000000000000
                                         -0.500000000000000  -0.250000000000000   0.250000000000000   0.500000000000000
                                         -0.500000000000000  -0.250000000000000   0.250000000000000   0.500000000000000
                                         -0.500000000000000  -0.250000000000000   0.250000000000000   0.500000000000000];
expSolCPHistoryLinear{1}(:, :, 3, 2) = [0     0     0     0
                                        0     0     0     0
                                        0     0     0     0
                                        0     0     0     0];
expSolCPHistoryLinear{1}(:, :, 4, 2) = [1     1     1     1
                                        1     1     1     1
                                        1     1     1     1
                                        1     1     1     1];
                           
% Define the expected solution in terms of the residual history
expSolResHistoryLinear = 5.270462766947300;
                           
% Define the expected solution in terms of the minimum element area size
expSolMinElArea = 2.500000000000000;

%% 8. Define the expected solutions for the nonlinear problem

% Define the expected solution for the Penalty method in terms of the
% Control Point displacements
expSolDispNonlinear = [   -0.000000000000000
                                           0
                                           0
                           0.220084786728229
                                           0
                                           0
                           0.660254360184687
                                           0
                                           0
                           0.880339146912915
                                           0
                                           0
                          -0.000000000000000
                                           0
                                           0
                           0.220084786728230
                                           0
                                           0
                           0.660254360184687
                                           0
                                           0
                           0.880339146912916
                                           0
                                           0
                           0.000000000000000
                                           0
                                           0
                           0.220084786728229
                                           0
                                           0
                           0.660254360184687
                                           0
                                           0
                           0.880339146912916
                                           0
                                           0
                           0.000000000000001
                                           0
                                           0
                           0.220084786728229
                                           0
                                           0
                           0.660254360184687
                                           0
                                           0
                           0.880339146912915
                                           0
                                           0];
                                   
% Define the expected solution in terms of the Control Point displacements
expSolCPHistoryNonlinear{1}(:, :, :, 1) = CP;
expSolCPHistoryNonlinear{1}(:, :, 1, 2) = [  -5.000000000000000  -5.000000000000000  -5.000000000000000  -4.999999999999999
                                             -2.279915213271771  -2.279915213271770  -2.279915213271771  -2.279915213271771
                                              3.160254360184687   3.160254360184687   3.160254360184687   3.160254360184687
                                              5.880339146912915   5.880339146912915   5.880339146912917   5.880339146912915];
expSolCPHistoryNonlinear{1}(:, :, 2, 2) = [  -0.500000000000000  -0.250000000000000   0.250000000000000   0.500000000000000
                                             -0.500000000000000  -0.250000000000000   0.250000000000000   0.500000000000000
                                             -0.500000000000000  -0.250000000000000   0.250000000000000   0.500000000000000
                                             -0.500000000000000  -0.250000000000000   0.250000000000000   0.500000000000000];
expSolCPHistoryNonlinear{1}(:, :, 3, 2) = [ 0     0     0     0
                                            0     0     0     0
                                            0     0     0     0
                                            0     0     0     0];
expSolCPHistoryNonlinear{1}(:, :, 4, 2) = [ 1     1     1     1
                                            1     1     1     1
                                            1     1     1     1
                                            1     1     1     1];
                                    
% Define the expected solution in terms of the residual history
expSolResHistoryNonlinear = [  5.270462766947300
                               0.816921728876840
                               0.012039004508595
                               0.000002755501612
                               0.000000000000410
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

%% 9. Define the linear analysis parameters
propLinearAnalysis.scheme = 'newtonRapshon';
propLinearAnalysis.noLoadSteps = 1;
propLinearAnalysis.eps = 10e-9;
propLinearAnalysis.maxIter = 1;

%% 10. Define the linear analysis parameters
propNonlinearAnalysis.scheme = 'newtonRapshon';
propNonlinearAnalysis.noLoadSteps = 1;
propNonlinearAnalysis.eps = 10e-9;
propNonlinearAnalysis.maxIter = 120;

%% 11. Solve the steady-state linear problem
plot_IGANLinear = 'undefined';
[dHatLinear, CPHistoryLinear, resHistoryLinear, ~, ~] = ...
    solve_IGAMembraneNLinear ...
    (BSplinePatch, propLinearAnalysis, solve_LinearSystem, ...
    plot_IGANLinear, graph,'');

%% 12. Solve the steady-state nonlinear problem
plot_IGANLinear = 'undefined';
[dHatNonlinear, CPHistoryNonlinear, resHistoryNonlinear, hasConverged, ~] = ...
    solve_IGAMembraneNLinear ...
    (BSplinePatch, propNonlinearAnalysis, solve_LinearSystem, ...
    plot_IGANLinear, graph, '');

%% 13. Verify the solution
testCase.verifyEqual(dHatLinear,expSolDispLinear, 'AbsTol', absTolDisp);
testCase.verifyEqual(CPHistoryLinear,expSolCPHistoryLinear, 'AbsTol', absTolDisp);
testCase.verifyEqual(resHistoryLinear,expSolResHistoryLinear, 'AbsTol', absTolResHistory);
testCase.verifyEqual(BSplinePatch.minElArea,expSolMinElArea, 'AbsTol', 0);
testCase.verifyEqual(dHatNonlinear,expSolDispNonlinear, 'AbsTol', absTolDisp);
testCase.verifyEqual(CPHistoryNonlinear,expSolCPHistoryNonlinear, 'AbsTol', absTolDisp);
testCase.verifyEqual(resHistoryNonlinear,expSolResHistoryNonlinear, 'AbsTol', absTolResHistory);
testCase.verifyEqual(hasConverged, true, 'AbsTol', 0);

end
