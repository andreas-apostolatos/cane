function testDDMNitschePlateMembraneAnalysis(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests the Nitsche method for both the multipatch coupling and the
% application of weak boundary conditions for a two-patch geometry each
% patch of which is discretized with one bilinear element. The analytical
% solution of the tip displacement is known to be 0.8803391469128941469 for
% the geometrically nonlinear problem.
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
% 6. Define the interface parametrization
%
% 7. Fill up the patch arrays
%
% 8. Define the expected solutions for the linear problem
%
% 9. Define the linear analysis parameters
%
% 10. Solve the steady-state nonlinear problem with the Penalty method
%
% 11. Solve the steady-state nonlinear problem using the Nitsche method
%
% 12. Verify the solution
%
%% Function main body

%% 0. Read input
absTol = 1e-12;

%% 1. Define the multipatch geometry in terms of NURBS

% Global variables:
Length = 10;
Width = 1;

% Patch 1 :
% _________

% Polynomial degrees
p1 = 1;
q1 = 1;

% Knot vectors
Xi1 = [0 0 1 1];
Eta1 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP1(:, :, 1) = [-Length/2 -Length/2
                0         0];

% y-coordinates
CP1(:, :, 2) = [-Width/2 Width/2
                -Width/2 Width/2];

% z-coordinates
CP1(:, :, 3) = [0 0
                0 0];
          
% Weights
CP1(:, :, 4) = [1 1
                1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = 0;
nxi1 = length(CP1(:, 1, 1));
neta1 = length(CP1(1, :, 1));
for i= 1:nxi1
    for j = 1:neta1
        if CP1(i, j, 4) ~= 1
            isNURBS1 = 1;
            break;
        end
    end
    if isNURBS1
        break;
    end
end

% Patch 2 :
% _________

% Polynomial degrees
p2 = 1;
q2 = 1;

% Knot vectors
Xi2 = [0 0 1 1];
Eta2 = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP2(:, :, 1) = [0        0
                Length/2 Length/2];

% y-coordinates
CP2(:, :, 2) = [-Width/2 Width/2
                -Width/2 Width/2];

% z-coordinates
CP2(:, :, 3) = [0 0
                0 0];
          
% Weights
CP2(:, :, 4) = [1 1
                1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS2 = 0;
nxi2 = length(CP2(:, 1, 1));
neta2 = length(CP2(1, :, 1));
for i = 1:nxi2
    for j = 1:neta2
        if CP2(i, j, 4) ~= 1
            isNURBS2 = 1;
            break;
        end
    end
    if isNURBS2
        break;
    end
end

%% 2. Define the material parameters

% Parameters
EYoung = 1e2;
poissonRatio = 0.0;
thickness = 1;
density = 7810;
prestress.voigtVector = [0
                         0
                         0];

% Patch 1 :
% _________

% Young's modulus
parameters1.E = EYoung;

% Poisson ratio
parameters1.nue = poissonRatio;

% Thickness of the plate
parameters1.t = thickness;

% Density of the membrane (used only for dynamics)
parameters1.rho = density;

% Prestress for the membrane
parameters1.prestress = prestress;

% Patch 2 :
% _________

% Young's modulus
parameters2.E = EYoung;

% Poisson ratio
parameters2.nue = poissonRatio;

% Thickness of the plate
parameters2.t = thickness;

% Density of the membrane (used only for dynamics)
parameters2.rho = density;

% Prestress for the membrane
parameters2.prestress = prestress;

%% 3. GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points

% Patch 1 :
% _________

int1.type = 'default';
if strcmp(int1.type, 'user')
    int1.xiNGP = 6;
    int1.etaNGP = 6;
    int1.xiNGPForLoad = 6;
    int1.etaNGPForLoad = 6;
    int1.nGPForLoad = 6;
    int1.nGPError = 12;
end

% Patch 2 :
% _________

int2.type = 'default';
if strcmp(int2.type, 'user')
    int2.xiNGP = 6;
    int2.etaNGP = 6;
    int2.xiNGPForLoad = 6;
    int2.etaNGPForLoad = 6;
    int2.nGPForLoad = 6;
    int2.nGPError = 12;
end

% Interface integration :
% _______________________

intC.type = 'user';
intC.method = 'Nitsche';
if strcmp(intC.type, 'user')
    if strcmp(intC.method, 'lagrangeMultipliers')
        intC.noGPsLambda = 12;
        intC.noGPsMu = 12;
    else
        intC.noGPs = 12;
    end
    intC.nGPsError = 12;
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

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

tp1 = 0;
tq1 = 0;
[Xi1, Eta1, CP1, p1, q1] = degreeElevateBSplineSurface ...
    (p1, q1, Xi1, Eta1, CP1, tp1, tq1, '');

% Patch 2 :
% _________

tp2 = 0;
tq2 = 0;
[Xi2, Eta2, CP2, p2, q2] = degreeElevateBSplineSurface ...
    (p2, q2, Xi2, Eta2, CP2, tp2, tq2, '');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = 1;
nKnotsXi1 = a;
nKnotsEta1 = ceil(a*Width/Length);
[Xi1, Eta1, CP1] = knotRefineUniformlyBSplineSurface ...
    (p1, Xi1, q1, Eta1, CP1, nKnotsXi1, nKnotsEta1, '');

% Patch 2 :
% _________

b = 1;
nKnotsXi2 = b;
nKnotsEta2 = ceil(b*Width/Length);
[Xi2, Eta2, CP2] = knotRefineUniformlyBSplineSurface ...
    (p2, Xi2, q2, Eta2, CP2, nKnotsXi2, nKnotsEta2, '');

%% 5. Define the boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous DOFs
homDOFs1 = [];
% xisup1 = [0 0];   
% etasup1 = [0 1];
% for dirSupp = 1:3
%     homDOFs1 = findDofs3D(homDOFs1,xisup1,etasup1,dirSupp,CP1);
% end
xisup1 = [0 1];   
etasup1 = [0 1];
for dirSupp = [2 3]
    homDOFs1 = findDofs3D(homDOFs1, xisup1, etasup1, dirSupp, CP1);
end

% Inhomogeneous DOFs
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak Dirichlet boundary conditions
weakDBC1 = [];
weakDBC1.noCnd = 1;
weakDBC1.method = 'nitsche';
weakDBC1.estimationStabilPrm = false;
weakDBC1.alpha = 0.0;
weakDBC1.xiExtension = {[0 0]};
weakDBC1.etaExtension = {[0 1]};
weakDBC1.int.type = 'default';
weakDBC1.int.noGPs = 16;

% Embedded cables
cables1.No = 0;

% Patch 2 :
% _________

homDOFs2 = [];
xisup2 = [0 1];
etasup2 = [0 1];
for dirSupp = [2 3]
    homDOFs2 = findDofs3D(homDOFs2, xisup2, etasup2, dirSupp, CP2);
end

% Inhomogeneous DOFs
inhomDOFs2 = [];
valuesInhomDOFs2 = [];

% Weak Dirichlet boundary conditions
weakDBC2 = [];

% Embedded cables
cables2.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load amplitude
FAmp = + 1e1;

% Patch 1 :
% _________

FAmp1 = FAmp*0;
NBC1.noCnd = 1;
xib1 = 1;   etab1 = [0 1];   dirForce1 = 'x';
NBC1.xiLoadExtension = {xib1};
NBC1.etaLoadExtension = {etab1};
NBC1.loadAmplitude{1} = FAmp1;
NBC1.loadDirection = {dirForce1};
NBC1.computeLoadVct{1} = 'computeLoadVctLineIGAThinStructure';
NBC1.isFollower(1, 1) = false;
NBC1.isTimeDependent(1, 1) = false;

% Patch 2 :
% _________

FAmp2 = FAmp;
NBC2.noCnd = 1;
xib2 = 1;   etab2 = [0 1];   dirForce2 = 'x';
NBC2.xiLoadExtension = {xib2};
NBC2.etaLoadExtension = {etab2};
NBC2.loadAmplitude{1} = FAmp2;
NBC2.loadDirection = {dirForce2};
NBC2.computeLoadVct{1} = 'computeLoadVctLineIGAThinStructure';
NBC2.isFollower(1, 1) = false;
NBC2.isTimeDependent(1, 1) = false;

%% 6. Define the interface parametrization

% Patch 1 :
% _________

% Connection with patch 2:
xicoup12 = [1 1];
etacoup12 = [0 1];

% Collect all interfaces into arrays:
xicoup1 = xicoup12;
etacoup1 = etacoup12;

% Patch 2 :
% _________

% Connection with patch 1:
xicoup21 = [0 0];
etacoup21 = [0 1];

% Collect all interfaces into arrays:
xicoup2 = xicoup21;
etacoup2 = etacoup21;

% Define connections by patch numbers
connections.xiEtaCoup(:,:) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21];
connections.No = length(connections.xiEtaCoup(:, 1));

%% 7. Fill up the patch arrays

% 1st patch :
% ___________

patch1 = fillUpPatch...
    (analysis, p1, Xi1, q1, Eta1, CP1, isNURBS1, parameters1, homDOFs1, ...
    inhomDOFs1, valuesInhomDOFs1, weakDBC1, cables1, NBC1, [], [], [], ...
    xicoup1, etacoup1, int1);

% 2nd patch :
% ___________

patch2 = fillUpPatch...
    (analysis, p2, Xi2, q2, Eta2, CP2, isNURBS2, parameters2, homDOFs2, ...
    inhomDOFs2, valuesInhomDOFs2, weakDBC2, cables2, NBC2, [], [], [], ...
    xicoup2, etacoup2, int2);

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch1 patch2};
noPatches = length(BSplinePatches);

%% 8. Define the expected solutions for the linear problem

% Define the expected solution for the Penalty method in terms of the
% Control Point displacements
expSolDispNonlinearPenalty = [ 0.000000000000000
                                               0
                                               0
                               0.440169573456458
                                               0
                                               0
                              -0.000000000000000
                                               0
                                               0
                               0.440169573456458
                                               0
                                               0
                               0.940169573456458
                                               0
                                               0
                               1.380339146912916
                                               0
                                               0
                               0.940169573456458
                                               0
                                               0
                               1.380339146912916
                                               0
                                               0];
                                           
% Define the expected solution for the Nitsche method in terms of the 
% Control Point displacements
expSolCPHistoryLinearPenalty{1}(:, :, :, 1) = CP1;
expSolCPHistoryLinearPenalty{1}(:, :, 1, 2) = [-5.000000000000000  -5.000000000000000
                                               0.440169573456458   0.440169573456458];                                
expSolCPHistoryLinearPenalty{1}(:, :, 2, 2) = [-0.500000000000000   0.500000000000000
                                               -0.500000000000000   0.500000000000000];
expSolCPHistoryLinearPenalty{1}(:, :, 3, 2) = [0 0
                                               0 0];
expSolCPHistoryLinearPenalty{1}(:, :, 4, 2) = [1 1
                                               1 1];
expSolCPHistoryLinearPenalty{2}(:, :, :, 1) = CP2;
expSolCPHistoryLinearPenalty{2}(:, :, 1, 2) = [0.940169573456458   0.940169573456458
                                               6.380339146912916   6.380339146912916];                                
expSolCPHistoryLinearPenalty{2}(:, :, 2, 2) = [-0.500000000000000   0.500000000000000
                                               -0.500000000000000   0.500000000000000];
expSolCPHistoryLinearPenalty{2}(:, :, 3, 2) = [0 0
                                               0 0];
expSolCPHistoryLinearPenalty{2}(:, :, 4, 2) = [1 1
                                               1 1];
                                
% Define the expected solution for the Nitsche method in terms of the 
% residual history
expSolResHistoryNonlinearPenalty = [   7.071067811865476
                                       1.898354550657023
                                       0.027976118380038
                                       0.000006403207080
                                       0.000000000000344];

% Define the expected solution for the Nitsche method in terms of the
% Control Point displacements
expSolDispNonlinearNitsche = [ 0.000000000000000
                                               0
                                               0
                               0.440169573456458
                                               0
                                               0
                              -0.000000000000000
                                               0
                                               0
                               0.440169573456458
                                               0
                                               0
                               0.440169573456458
                                               0
                                               0
                               0.880339146912915
                                               0
                                               0
                               0.440169573456458
                                               0
                                               0
                               0.880339146912916
                                               0
                                               0];
                                   
% Define the expected solution for the Nitsche method in terms of the 
% Control Point displacements
expSolCPHistoryLinearNitsche{1}(:, :, :, 1) = CP1;
expSolCPHistoryLinearNitsche{1}(:, :, 1, 2) = [-5.000000000000000  -4.999999999999999
                                               0.440169573456458   0.440169573456458];                                 
expSolCPHistoryLinearNitsche{1}(:, :, 2, 2) = [-0.500000000000000   0.500000000000000
                                               -0.500000000000000   0.500000000000000]; 
expSolCPHistoryLinearNitsche{1}(:, :, 3, 2) = [0 0
                                               0 0];
expSolCPHistoryLinearNitsche{1}(:, :, 4, 2) = [1 1
                                               1 1];
expSolCPHistoryLinearNitsche{2}(:, :, :, 1) = CP2;
expSolCPHistoryLinearNitsche{2}(:, :, 1, 2) = [0.440169573456457   0.440169573456457
                                               5.880339146912915   5.880339146912915];                                 
expSolCPHistoryLinearNitsche{2}(:, :, 2, 2) = [-0.500000000000000   0.500000000000000
                                               -0.500000000000000   0.500000000000000]; 
expSolCPHistoryLinearNitsche{2}(:, :, 3, 2) = [0 0
                                               0 0];
expSolCPHistoryLinearNitsche{2}(:, :, 4, 2) = [1 1
                                               1 1];

% Define the expected solution for the Nitsche method in terms of the 
% residual history
expSolResHistoryNonlinearNitsche = [   7.071067811865476
                                       1.096015510839187
                                       0.016152019477600
                                       0.000003696893340
                                       0.000000000000189];
                                  
%% 9. Define the linear analysis parameters
propNLinearAnalysis.scheme = 'newtonRapshon';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 1e-9;
propNLinearAnalysis.maxIter = 5;
propNLinearAnalysis.conservativeLoad = true(noPatches, 1);

%% 10. Solve the steady-state nonlinear problem with the Penalty method
Dm1 = parameters1.E*parameters1.t/(1-parameters1.nue^2)*...
      [1              parameters1.nue 0
       parameters1.nue 1              0
       0               0              (1-parameters1.nue)/2];
Dm2 = parameters2.E*parameters2.t/(1-parameters2.nue^2)*...
      [1              parameters2.nue 0
       parameters2.nue 1              0
       0               0              (1-parameters2.nue)/2];
propCouplingPenalty.method = 'penalty';
propCouplingPenalty.alphaD = max([norm(Dm1) norm(Dm2)])*...
    (1/min(BSplinePatches{1}.minElArea,BSplinePatches{2}.minElArea));
propCouplingPenalty.alphaR = 0;
propCouplingPenalty.intC = intC;
plot_IGANLinear = 'undefined';
[dHatPenalty, CPHistoryPenalty, resHistoryPenalty, ~, ~, ~, ~, ~, ~, ~, ~] = ...
    solve_DDMIGAMembraneMultipatchesNLinear...
    (BSplinePatches, connections, propCouplingPenalty, ...
    propNLinearAnalysis, solve_LinearSystem, plot_IGANLinear, graph, '');

%% 11. Solve the steady-state nonlinear problem using the Nitsche method
propCouplingNitsche.method = 'nitsche';
propCouplingNitsche.estimationStabilPrm = false;
propCouplingNitsche.alphaD = 0;
propCouplingNitsche.alphaR = 0;
propCouplingNitsche.gammaTilde = 0.5;
propCouplingNitsche.intC = intC;
plot_IGANLinear = 'undefined';
[dHatNitsche, CPHistoryNitsche, resHistoryNitsche, ~, ~, ~, ~, ~, ~, ~, ~] = ...
    solve_DDMIGAMembraneMultipatchesNLinear ...
    (BSplinePatches, connections, propCouplingNitsche, ...
    propNLinearAnalysis, solve_LinearSystem, plot_IGANLinear, graph, '');

%% 12. Verify the solution
testCase.verifyEqual(dHatPenalty, expSolDispNonlinearPenalty, 'AbsTol', absTol);
testCase.verifyEqual(CPHistoryPenalty, expSolCPHistoryLinearPenalty, 'AbsTol', absTol);
testCase.verifyEqual(resHistoryPenalty, expSolResHistoryNonlinearPenalty, 'AbsTol', absTol);
testCase.verifyEqual(dHatNitsche, expSolDispNonlinearNitsche, 'AbsTol', absTol);
testCase.verifyEqual(CPHistoryNitsche, expSolCPHistoryLinearNitsche, 'AbsTol', absTol);
testCase.verifyEqual(resHistoryNitsche, expSolResHistoryNonlinearNitsche, 'AbsTol', absTol);

end
