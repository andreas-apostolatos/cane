function testNonlinearMembraneAnalysis(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Test the isogeometric membrane geometically nonlinear analysis over a
% circular membrane subject to pressure load.
%
% 0. Read input
%
% 1. Define the geometry in terms of NURBS
%
% 2. Define the material constants
%
% 3. GUI
%
% 4. Define the refinement
%
% 5. Define the Dirichlet and the Neumann boundary conditions
%
% 6. Create the B-Spline patch array
%
% 7. Define the expected solutions
%
% 8. Define the nonlinear analysis properties
%
% 9. Solve the nonlinear system
%
% 10. Verify the results
%
%% Function main body

%% 0. Read input

% Define absolute tolerance
absTol = 1e-12;
absTolRelaxed8 = absTol*1e8;

%% 1. Define the geometry in terms of NURBS

% Global variables:
Radius = .5;
edgeSquare = Radius*sqrt(2);

% Polynomial degrees
p = 2;
q = 2;

% Knot vectors
Xi = [0 0 0 1 1 1];
Eta = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [-edgeSquare/2 0 edgeSquare/2
             -edgeSquare   0 edgeSquare
             -edgeSquare/2 0 edgeSquare/2];
         
% y-coordinates
CP(:,:,2) = [-edgeSquare/2 -edgeSquare   -edgeSquare/2
             0             0             0
             edgeSquare/2  edgeSquare    edgeSquare/2];
         
% z-coordinates
CP(:,:,3) = [0 0 0
             0 0 0
             0 0 0];
       
% Weights
weight = sqrt(2)/2;
CP(:,:,4) = [1      weight 1
             weight 1      weight
             1      weight 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i= 1:nxi
    for j=1:neta
        if CP(i,j,4)~=1
            isNURBS = 1;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% 2. Define the material constants

% Young's modulus
parameters.E = 1.8e6;

% Poisson ratio
parameters.nue = 0.3;

% Thickness of the shell
parameters.t = .01;

% Density of the shell (used only for dynamics)
parameters.rho = 7810;

% Prestress for the membrane
parameters.prestress.voigtVector = [10/parameters.t
                                    10/parameters.t
                                    0];
                    
%% 3. GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'user';
if strcmp(int.type,'user')
    int.xiNGP = p + 1;
    int.etaNGP = q + 1;
    int.xiNGPForLoad = ceil((p + 1)/2);
    int.etaNGPForLoad = ceil((q + 1)/2);
end

% Choose linear equation solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

%% 4. Define the refinement

% Degree by which to elevate
a = 0;
tp = a;
tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'');

% Number of knots to exist in both directions
scaling = 3;
refXi = scaling;
refEta = scaling;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'');

%% 5. Define the Dirichlet and the Neumann boundary conditions

% supports (Dirichlet boundary conditions)

% Homogeneous Dirichlet boundary conditions
homDOFs = [];
xiSup = [0 0];   etaSup = [0 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [0 1];   etaSup = [0 0];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [0 1];   etaSup = [1 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC = [];

% Embedded cables
cables.No = 0;

% load (Neuman boundary conditions)
FAmp = + 1e3;
xib = [0 1];   etab = [0 1];   dirForce = 'z';
NBC.noCnd = 1;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.computeLoadVct{1} = 'computeLoadVctAreaIGAThinStructure';
NBC.isFollower(1,1) = false;
NBC.isTimeDependent(1,1) = false;

%% 6. Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,...
    valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);

%% 7. Define the expected solutions

% Define the expected solution in terms of Control Point displacements
expSolDisp = [                 0
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
              -0.005083184917132
              -0.005083184917132
               0.043960556784785
              -0.009047081812731
              -0.000000000000000
               0.076306274644207
              -0.005083184917132
               0.005083184917132
               0.043960556784785
                               0
                               0
                               0
                               0
                               0
                               0
              -0.000000000000000
              -0.009047081812731
               0.076306274644207
               0.000000000000000
              -0.000000000000000
               0.119878327387778
              -0.000000000000000
               0.009047081812731
               0.076306274644207
                               0
                               0
                               0
                               0
                               0
                               0
               0.005083184917132
              -0.005083184917132
               0.043960556784785
               0.009047081812731
               0.000000000000000
               0.076306274644207
               0.005083184917132
               0.005083184917132
               0.043960556784785
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

% Define the expected solution in terms of the Control Point deformation
expSolCPHist = struc([]);
for counterPatches = 1
    for counterComp = 1:4
        expSolCPHist{counterPatches}(:,:,counterComp,1) = CP(:,:,counterComp);
    end
end
expSolCPHist{1}(:,:,1,2) = [  -0.353553390593274  -0.261203874963741  -0.000000000000000   0.261203874963741   0.353553390593274
                              -0.445902906222806  -0.313473813571207  -0.000000000000000   0.313473813571207   0.445902906222806
                              -0.519434138474439  -0.350907334469663  -0.000000000000000   0.350907334469663   0.519434138474439
                              -0.445902906222806  -0.313473813571207  -0.000000000000000   0.313473813571207   0.445902906222806
                              -0.353553390593274  -0.261203874963741  -0.000000000000000   0.261203874963741   0.353553390593274];
expSolCPHist{1}(:,:,2,2) = [  -0.353553390593274  -0.445902906222806  -0.519434138474439  -0.445902906222806  -0.353553390593274
                              -0.261203874963741  -0.313473813571207  -0.350907334469663  -0.313473813571207  -0.261203874963741
                              -0.000000000000000  -0.000000000000000  -0.000000000000000  -0.000000000000000  -0.000000000000000
                               0.261203874963741   0.313473813571207   0.350907334469663   0.313473813571207   0.261203874963741
                               0.353553390593274   0.445902906222806   0.519434138474439   0.445902906222806   0.353553390593274];
expSolCPHist{1}(:,:,3,2) = [                   0                   0                   0                   0                   0
                                               0   0.043960556784785   0.076306274644207   0.043960556784785                   0
                                               0   0.076306274644207   0.119878327387778   0.076306274644207                   0
                                               0   0.043960556784785   0.076306274644207   0.043960556784785                   0
                                               0                   0                   0                   0                   0];
expSolCPHist{1}(:,:,4,2) = [   1.000000000000000   0.902368927062183   0.837281545103638   0.902368927062183   1.000000000000000
                               0.902368927062183   0.869825236082910   0.848129442096728   0.869825236082910   0.902368927062183
                               0.837281545103638   0.848129442096728   0.855361373425455   0.848129442096728   0.837281545103638
                               0.902368927062183   0.869825236082910   0.848129442096728   0.869825236082910   0.902368927062183
                               1.000000000000000   0.902368927062183   0.837281545103638   0.902368927062183   1.000000000000000];

% Define the reference solution in terms of the residual history
expSolResHistory = 1e7*[   0.000018575626951
                           8.097896618204144
                           2.398409238396341
                           0.709990385646885
                           0.209931787911792
                           0.061912724337103
                           0.018161526319530
                           0.005281204403340
                           0.001526917578831
                           0.000450294792041
                           0.000141217192007
                           0.000046882200547
                           0.000014115973783
                           0.000002756428872
                           0.000000184102697
                           0.000000001105621
                           0.000000000000047
                           0.000000000000000
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
                                       
% Define the expected solution in terms of the minimum element area size
expSolMinElAreaSize = 0.069836579539244;
                           
%% 8. Define the nonlinear analysis properties
propNLinearAnalysis.method = 'newtonRapshon';
propNLinearAnalysis.noLoadSteps = 1;
propNLinearAnalysis.eps = 10e-9;
propNLinearAnalysis.maxIter = 120;

%% 9. Solve the nonlinear system
plot_IGANLinear = '';
[dHatNLinear,CPHistory,resHistory,hasConverged,~] = ...
    solve_IGAMembraneNLinear...
    (BSplinePatch,propNLinearAnalysis,solve_LinearSystem,...
    plot_IGANLinear,[],'');

%% 10. Verify the results
testCase.verifyEqual(dHatNLinear,expSolDisp,'AbsTol',absTol);
testCase.verifyEqual(CPHistory,expSolCPHist,'AbsTol',absTol);
testCase.verifyEqual(resHistory,expSolResHistory,'AbsTol',absTolRelaxed8);
testCase.verifyEqual(hasConverged,true);
testCase.verifyEqual(BSplinePatch.minElArea,expSolMinElAreaSize,'AbsTol',absTol);

end
