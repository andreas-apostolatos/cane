function testDDMFourPointSailAnalysis(testCase)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Tests all methods for the form-finding analysis over the 3-patch four
% point sail.
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
% 6. Create the patch cell arrays
%
% 7. Create the Lagrange Multipliers fields for all interfaces
%
% 8. Define the expected solutions
%
% 9. Define the form-finding analysis properties
%
% 10. Solve the multipatch problem with Penalty
%
% 11. Solve the multipatch problem with Lagrange Multipliers
%
% 12. Solve the multipatch problem with Mortar
%
% 13. Solve the multipatch problem with Nitsche
%
% 14. Verify the solution
%
%% Function main body

%% 0. Read input

% Define absolute tolerances
absTolDisp = 1e-11;
absTolResHistory = absTolDisp*1e1;
absTolStabil = absTolDisp*1e7;

%% 1. Define the multipatch geometry in terms of NURBS

% Global variables:
Length = 20;
Width = Length;
Height = Length/2;

% Patch 1 :
% _________

% Polynomial degrees
p1 = 1;
q1 = 2;

% Knot vectors
Xi1 = [0 0 1 1];
Eta1 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP1(:, :, 1) = [-Length/2 -Length/2 -Length/2
                -Length/2 -Length/3 Length/2];
         
% y-coordinates
CP1(:, :, 2) = [-Width/2 -Width/2 -Width/2
                Width/2  -Width/3 -Width/2];
         
% z-coordinates
CP1(:, :, 3) = [0      0                 0
                Height abs(CP1(2,2,1))/2 Height];
       
% Weights
CP1(:, :, 4) = [1 1 1
                1 1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = false;
nxi1 = length(CP1(:, 1, 1));
neta1 = length(CP1(1, :, 1));
for i = 1:nxi1
    for j = 1:neta1
        if CP1(i, j, 4) ~= 1
            isNURBS1 = true;
            break;
        end
    end
    if isNURBS1
        break;
    end
end

% Patch 3 :
% _________

% Polynomial degrees
p3 = 1;
q3 = 2;

% Knot vectors
Xi3 = [0 0 1 1];
Eta3 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP3(:, :, 1) = [-Length/2 Length/3 Length/2
                Length/2  Length/2 Length/2];
         
% y-coordinates
CP3(:, :, 2) = [Width/2 Width/3 -Width/2
                Width/2 Width/2 Width/2];
         
% z-coordinates
CP3(:, :, 3) = [Height abs(CP3(1,2,1))/2 Height
                0      0                 0];
       
% Weights
CP3(:, :, 4) = [1 1 1
                1 1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS3 = false;
nxi3 = length(CP3(:, 1, 1));
neta3 = length(CP3(1, :, 1));
for i = 1:nxi3
    for j = 1:neta3
        if CP3(i, j, 4) ~= 1
            isNURBS3 = true;
            break;
        end
    end
    if isNURBS3
        break;
    end
end

% Patch 2 :
% _________

% Polynomial degrees
p2 = 1;
q2 = 2;

% Knot vectors
Xi2 = [0 0 1 1];
Eta2 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP2(:, :, 1) = [CP1(2, 1, 1) CP1(2, 2, 1) CP1(2, 3, 1)
                CP3(1, 1, 1) CP3(1, 2, 1) CP3(1, 3, 1)];
         
% y-coordinates
CP2(:, :, 2) = [CP1(2, 1, 2) CP1(2, 2, 2) CP1(2, 3, 2)
                CP3(1, 1, 2) CP3(1, 2, 2) CP3(1, 3, 2)];
         
% z-coordinates
CP2(:, :, 3) = [CP1(2, 1, 3) CP1(2, 2, 3) CP1(2, 3, 3)
                CP3(1, 1, 3) CP3(1, 2, 3) CP3(1, 3, 3)];
       
% Weights
CP2(: ,:, 4) = [CP1(2, 1, 4) CP1(2, 2, 4) CP1(2, 3, 4)
                CP3(1, 1, 4) CP3(1, 2, 4) CP3(1, 3, 4)];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS2 = false;
nxi2 = length(CP2(:, 1, 1));
neta2 = length(CP2(1, :, 1));
for i = 1:nxi2
    for j = 1:neta2
        if CP2(i, j, 4)~=1
            isNURBS2 = true;
            break;
        end
    end
    if isNURBS2
        break;
    end
end

%% 2. Define the material parameters

% general parameters
EYoung = 8e+8;
nue = .4;
thickness = 1e-3;
sigma0 = 3e+3/thickness;
prestress.voigtVector = [sigma0
                         sigma0
                         0];
density = 8050;

% Patch 1 :
% _________

% Young's modulus
parameters1.E = EYoung;

% Poisson ratio
parameters1.nue = nue;

% Thickness of the membrane
parameters1.t = thickness;

% Density of the membrane
parameters1.rho = density;

% Prestress of the membrane
parameters1.prestress = prestress;

% Patch 2 :
% _________

% Young's modulus
parameters2.E = EYoung;

% Poisson ratio
parameters2.nue = nue;

% Thickness of the plate
parameters2.t = thickness;

% Density of the membrane
parameters2.rho = density;

% Prestress of the membrane
parameters2.prestress = prestress;

% Patch 3 :
% _________

% Young's modulus
parameters3.E = EYoung;

% Poisson ratio
parameters3.nue = nue;

% Thickness of the plate
parameters3.t = thickness;

% Density of the membrane
parameters3.rho = density;

% Prestress of the membrane
parameters3.prestress = prestress;

% Cable :
% _______

parametersCable.E = 1.6e+11;
parametersCable.radiusCS = 12e-3/2;
parametersCable.areaCS = pi*parametersCable.radiusCS^2;
parametersCable.rho = 8050;
parametersCable.prestress = 6e+4/parametersCable.areaCS;

%% 3. GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points

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

% Patch 3 :
% _________

int3.type = 'default';
if strcmp(int3.type, 'user')
    int3.xiNGP = 6;
    int3.etaNGP = 6;
    int3.xiNGPForLoad = 6;
    int3.etaNGPForLoad = 6;
    int3.nGPForLoad = 6;
    int3.nGPError = 12;
end

% Interface integration :
% _______________________

intC.type = 'user';
intC.method = 'Nitsche';
if strcmp(intC.type, 'user')
    if strcmp(intC.method, 'lagrangeMultipliers')
        intC.nGP1 = 16;
        intC.nGP2 = 16;
    else
        intC.noGPs = 16;
    end
    intC.nGPError = 16;
end

% Define the coupling properties

% Patch 1 :
% _________

Dm1 = parameters1.E*parameters1.t/(1 - parameters1.nue^2)*...
      [1              parameters1.nue 0
       parameters1.nue 1              0
       0               0              (1 - parameters1.nue)/2];
   
% Patch 2 :
% _________

Dm2 = parameters2.E*parameters2.t/(1 - parameters2.nue^2)*...
      [1              parameters2.nue 0
       parameters2.nue 1              0
       0               0              (1 - parameters2.nue)/2];

% Patch 3 :
% _________

Dm3 = parameters3.E*parameters3.t/(1 - parameters3.nue^2)*...
      [1              parameters3.nue 0
       parameters3.nue 1              0
       0               0              (1 - parameters3.nue)/2];

%% 4. Define the refinement

% Select p-refinement level
iPRef = 0;

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = iPRef;
tp1 = a + 1;
tq1 = a;
[Xi1, Eta1, CP1, p1, q1] = degreeElevateBSplineSurface ...
    (p1, q1, Xi1, Eta1, CP1, tp1, tq1, '');

% Patch 2 :
% _________

b = iPRef + 1;
tp2 = b + 1;
tq2 = b;
[Xi2, Eta2, CP2, p2, q2] = ...
    degreeElevateBSplineSurface(p2, q2, Xi2, Eta2, CP2, tp2, tq2, '');

% Patch 3 :
% _________

c = iPRef;
tp3 = c + 1;
tq3 = c;
[Xi3, Eta3, CP3, p3, q3] = degreeElevateBSplineSurface ...
    (p3, q3, Xi3, Eta3, CP3, tp3, tq3, '');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% General parameter used for preserving the produced elements 
% shaped as close as possible to a rectangle
edgeRatio = ceil((2/3)*Length/Width/(sin(4*pi/9))/2);

% Select h-refinement level
iHRef = 1;

% Patch 1 :
% _________

noKnotsXi1 = ceil((14/9)*(iHRef + 1));
noKnotsEta1 = edgeRatio*noKnotsXi1;
[Xi1, Eta1, CP1] = knotRefineUniformlyBSplineSurface ...
    (p1, Xi1, q1, Eta1, CP1, noKnotsXi1, noKnotsEta1, '');

% Patch 2 :
% _________

noKnotsXi2 = ceil((iHRef + 1)*(10/15));
noKnotsEta2 = ceil(4*edgeRatio*noKnotsXi2*(2/3));
[Xi2, Eta2, CP2] = knotRefineUniformlyBSplineSurface ...
    (p2, Xi2, q2, Eta2, CP2, noKnotsXi2, noKnotsEta2,'');

% Patch 3 :
% _________

noKnotsXi3 = ceil((9/6)*(iHRef + 1));
noKnotsEta3 = edgeRatio*noKnotsXi3;
[Xi3, Eta3, CP3] = knotRefineUniformlyBSplineSurface ...
    (p3, Xi3, q3, Eta3, CP3, noKnotsXi3, noKnotsEta3, '');

%% 5. Define the boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs1 = [];
xisup1 = [0 0];
etasup1 = [0 1];
for dir = 1:3
    homDOFs1 = findDofs3D ...
        (homDOFs1, xisup1, etasup1, dir, CP1);
end
xisup1 = [1 1];
etasup1 = [0 0];
for dir = 1:3
    homDOFs1 = findDofs3D ...
        (homDOFs1, xisup1, etasup1, dir, CP1);
end
xisup1 = [1 1];
etasup1 = [1 1];
for dir = 1:3
    homDOFs1 = findDofs3D ...
        (homDOFs1, xisup1, etasup1, dir, CP1);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC1.noCnd = 0;

% Embedded cables
cables1.No = 2;
cables1.xiExtension = {[0 1] [0 1]};
cables1.etaExtension = {[0 0] [1 1]};
cables1.parameters = {parametersCable parametersCable parametersCable parametersCable};
cables1.int.type = 'default';
% cables1.int.type = 'user';
cables1.int.noGPs = 16;

% Patch 2 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs2 = [];
xisup2 = [0 1];
etasup2 = [0 0];
for dir = 1:3
    homDOFs2 = findDofs3D ...
        (homDOFs2, xisup2, etasup2, dir, CP2);
end
xisup2 = [0 1];
etasup2 = [1 1];
for dir = 1:3
    homDOFs2 = findDofs3D ...
        (homDOFs2, xisup2, etasup2, dir, CP2);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs2 = [];
valuesInhomDOFs2 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC2.noCnd = 0;

% Embedded cables
cables2.No = 0;

% Patch 3 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs3 = [];
xisup3 = [0 0];
etasup3 = [0 0];
for dir = 1:3
    homDOFs3 = findDofs3D ...
        (homDOFs3, xisup3, etasup3, dir, CP3);
end
xisup3 = [0 0];
etasup3 = [1 1];
for dir = 1:3
    homDOFs3 = findDofs3D ...
        (homDOFs3, xisup3, etasup3, dir, CP3);
end
xisup3 = [1 1];
etasup3 = [0 1];
for dir = 1:3
    homDOFs3 = findDofs3D ...
        (homDOFs3, xisup3, etasup3, dir, CP3);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs3 = [];
valuesInhomDOFs3 = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC3.noCnd = 0;

% Embedded cables
cables3.No = 2;
cables3.xiExtension = {[0 1] [0 1]};
cables3.etaExtension = {[0 0] [1 1]};
cables3.parameters = {parametersCable parametersCable parametersCable parametersCable};
cables3.int.type = 'default';
% cables3.int.type = 'user';
cables3.int.noGPs = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parameter
loadAmplitude = 0;

% Patch 1 :
% _________

FAmp1 = loadAmplitude;
NBC1.noCnd = 1;
xib1 = [0 1];   etab1 = [0 1];   dirForce1 = 'z';
NBC1.xiLoadExtension = {xib1};
NBC1.etaLoadExtension = {etab1};
NBC1.loadAmplitude = {FAmp1};
NBC1.loadDirection = {dirForce1};
NBC1.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC1.isFollower = false;
NBC1.isTimeDependent = true;

% Patch 2 :
% _________

FAmp2 = - loadAmplitude;
NBC2.noCnd = 1;
xib2 = [0 1];   etab2 = [0 1];   dirForce2 = 'z';
NBC2.xiLoadExtension = {xib2};
NBC2.etaLoadExtension = {etab2};
NBC2.loadAmplitude = {FAmp2};
NBC2.loadDirection = {dirForce2};
NBC2.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC2.isFollower = false;
NBC2.isTimeDependent = true;

% Patch 3 :
% _________

FAmp3 = - loadAmplitude;
NBC3.noCnd = 1;
xib3 = [0 1];   etab3 = [0 1];   dirForce3 = 'z';
NBC3.xiLoadExtension = {xib3};
NBC3.etaLoadExtension = {etab3};
NBC3.loadAmplitude = {FAmp3};
NBC3.loadDirection = {dirForce3};
NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC3.isFollower = false;
NBC3.isTimeDependent = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface parametrizations     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Connection with patch 2:
xicoup12 = [1 1];   etacoup12 = [0 1];

% Collect all interfaces into arrays:
xicoup1 = xicoup12;
etacoup1 = etacoup12;

% Patch 2 :
% _________

% Connection with patch 1:
xicoup21 = [0 0];   etacoup21 = [0 1];

% Connection with patch 3:
xicoup23 = [1 1];   etacoup23 = [0 1];

% Collect all interfaces into arrays:
xicoup2 = [xicoup21 xicoup23];
etacoup2 = [etacoup21 etacoup23];

% Patch 3 :
% _________

% Connection with patch 2:
xicoup32 = [0 0];   etacoup32 = [0 1];

% Collect all interfaces into arrays:
xicoup3 = xicoup32;
etacoup3 = etacoup32;

% Define connections :
% ____________________

% Number of connections
noConnections = 2;

% Define connections by patch numbers
connections.No = noConnections;
connections.xiEtaCoup = zeros(noConnections,10);
connections.xiEtaCoup(:,:) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21
                              2 3 xicoup23 etacoup23 xicoup32 etacoup32];
connections.masterSlave = [ false
                            true];
connectionsLM = connections;

%% 6. Create the patch cell arrays

% 1st patch :
% ___________

patch1 = fillUpPatch ...
    (analysis, p1, Xi1, q1, Eta1, CP1, isNURBS1, parameters1, homDOFs1, ...
    inhomDOFs1, valuesInhomDOFs1, weakDBC1, cables1, NBC1, [], [], [], ...
    xicoup1, etacoup1, int1);

% 2nd patch :
% ___________

patch2 = fillUpPatch ...
    (analysis, p2, Xi2, q2, Eta2, CP2, isNURBS2, parameters2, homDOFs2, ...
    inhomDOFs2, valuesInhomDOFs2, weakDBC2, cables2, NBC2, [], [], [], ...
    xicoup2, etacoup2, int2);

% 2nd patch :
% ___________

patch3 = fillUpPatch ...
    (analysis, p3, Xi3, q3, Eta3, CP3, isNURBS3, parameters3, homDOFs3, ...
    inhomDOFs3, valuesInhomDOFs3, weakDBC3, cables3, NBC3, [], [], [], ...
    xicoup3, etacoup3, int3);

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch1 patch2 patch3};
noPatches = length(BSplinePatches);

%% 7. Create the Lagrange Multipliers fields for all interfaces
for iConnections = 1:connections.No
    %% Get the IDs of the patches involved
    idI = connections.xiEtaCoup(iConnections,1);
    idJ = connections.xiEtaCoup(iConnections,2);
    
    %% Create a basic Lagrange Multipliers field
    pLambda = 0;
    XiLambda = [0 1];
    CPLambda(:,4) = [1];
    isNURBSLambda = 0;
    nxiLambda = length(CPLambda(:, 1, 1));
    for i = 1:nxiLambda
        if CPLambda(i,4) ~= 1
            isNURBSLambda = 1;
            break;
        end
    end
    
    %% Find the interface parametrization for the involved patches
    xiCoupI = connections.xiEtaCoup(iConnections, 3:4);
    etaCoupI = connections.xiEtaCoup(iConnections, 5:6);
    if xiCoupI(1, 1) ~= xiCoupI(1, 2) && etaCoupI(1, 1) == etaCoupI(1, 2)
        isOnXiI = true;
    elseif xiCoupI(1, 1) == xiCoupI(1, 2) && etaCoupI(1, 1) ~= etaCoupI(1, 2)
        isOnXiI = false;
    else
        error('Either the xi or the eta direction has to be fixed along the interface for patch %d', idI);
    end
    xiCoupJ = connections.xiEtaCoup(iConnections, 7:8);
    etaCoupJ = connections.xiEtaCoup(iConnections, 9:10);
    if xiCoupJ(1, 1) ~= xiCoupJ(1, 2) && etaCoupJ(1, 1) == etaCoupJ(1, 2)
        isOnXiJ = true;
    elseif xiCoupJ(1, 1) == xiCoupJ(1, 2) && etaCoupJ(1, 1) ~= etaCoupJ(1, 2)
        isOnXiJ = false;
    else
        error('Either the xi or the eta direction has to be fixed along the interface for patch %d', idJ);
    end

    %% Degree elevate the Lagrange Multipliers field
    if isOnXiI
        polOrderI = BSplinePatches{idI}.p;
    else
        polOrderI = BSplinePatches{idI}.q;
    end
    if isOnXiJ
        polOrderJ = BSplinePatches{idJ}.p;
    else
        polOrderJ = BSplinePatches{idJ}.q;
    end
    pLM = min(polOrderI,polOrderJ);
    if pLM > 0
        clear pLambda XiLambda CPLambda;
        pLambda = 1;
        XiLambda = [0 0 1 1];
        CPLambda(:, 4) = [1 1];
        isNURBSLambda = 0;
        nxiLambda = length(CPLambda(:, 1, 1));
        for i = 1:nxiLambda
            if CPLambda(i, 4) ~= 1
                isNURBSLambda = 1;
                break;
            end
        end
        
        % Perform accordingly a p-refinement
        tpLambda = pLM - pLambda;
        [XiLambda, CPLambda, pLambda] = degreeElevateBSplineCurve ...
            (pLambda, XiLambda, CPLambda, tpLambda, '');
    end
    
    %% Perform a knot insertion to the Lagrange Multipliers field
    if isOnXiI
        noKnotsI = length(unique(BSplinePatches{idI}.Xi)) - 2;
    else
        noKnotsI = length(unique(BSplinePatches{idI}.Eta)) - 2;
    end
    if isOnXiJ
        noKnotsJ = length(unique(BSplinePatches{idJ}.Xi)) - 2;
    else
        noKnotsJ = length(unique(BSplinePatches{idJ}.Eta)) - 2;
    end
    scaleLM = 2.0;
    noLambda = ceil(min([noKnotsI noKnotsJ])*scaleLM);
    [XiLambdaLM, CPLambdaLM] = knotRefineUniformlyBSplineCurve ...
        (noLambda, pLambda, XiLambda, CPLambda, '');
    
    %% Fill up the Lagrange Multipliers patch and add it to the array
    lambdaLM = fillUpLagrangeMultipliers ...
        (pLambda, XiLambdaLM, CPLambdaLM, isNURBSLambda);
    connectionsLM.lambda{iConnections} = lambdaLM;
end

%% 8. Define the expected solutions

% Define the expected solution in terms of the Penalty solution
expSolCPPenalty{1}(:, :, 1) = [ -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000
                                -8.228109707252313  -8.497892788990699  -8.697655980935862  -8.527729909501268  -7.997412987346291  -7.553914791799728
                                -6.723192887296725  -6.853343429091908  -6.856827804007471  -6.043815231350998  -4.095470376590916  -2.610256178381336
                                -6.684780663204148  -6.202775838161272  -5.530718649364709  -3.764070171918775  -0.191551495944527   2.363993975388313
                                -8.165055357378781  -6.231548606532571  -4.632523519751421  -1.469180789243062   3.657854215018534   7.439562716266470
                               -10.000000000000000  -6.622348648746447  -4.300183680106834  -0.303667745237445   5.439611031065374  10.000000000000000];
expSolCPPenalty{1}(:, :, 2) = [ -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000
                                -7.553914791799749  -7.997412987346309  -8.527729909501277  -8.697655980935858  -8.497892788990674  -8.228109707252299
                                -2.610256178381427  -4.095470376590973  -6.043815231351016  -6.856827804007462  -6.853343429091907  -6.723192887296705
                                2.363993975388058  -0.191551495944683  -3.764070171918764  -5.530718649364649  -6.202775838161246  -6.684780663204165
                                7.439562716266193   3.657854215018084  -1.469180789243007  -4.632523519751155  -6.231548606532903  -8.165055357378947
                                10.000000000000000   5.439611031064631  -0.303667745237386  -4.300183680106343  -6.622348648747567 -10.000000000000000];
expSolCPPenalty{1}(: ,: ,3) = [ 0                   0                   0                   0                   0                   0
                                1.552191531747978   1.293193472531260   1.020954317046266   1.020954317046271   1.293193472531275   1.552191531747992
                                3.910151036510069   3.303698680974408   2.606254681420186   2.606254681420198   3.303698680974430   3.910151036510110
                                5.987243124769289   4.876022596602849   3.749521904365741   3.749521904365732   4.876022596602915   5.987243124769402
                                8.370706990172344   6.266249873066756   4.477020284029108   4.477020284029009   6.266249873067017   8.370706990172474
                                10.000000000000000   7.006746318054729   4.687545937150913   4.687545937150716   7.006746318055358  10.000000000000000];
expSolCPPenalty{1}(:, :, 4) = [    1     1     1     1     1     1
                                   1     1     1     1     1     1
                                   1     1     1     1     1     1
                                   1     1     1     1     1     1
                                   1     1     1     1     1     1
                                   1     1     1     1     1     1];
expSolCPPenalty{2}(:, :, 1) = [-10.000000000000000  -8.485792391222450  -6.097805117887147  -4.743506998721749  -2.406903665584097   0.628736753958042   4.417021629154970   7.991439381480546  10.000000000000000
                               -10.000000000000000  -8.413524951789617  -5.989013583865464  -4.012273034120131  -1.603017925308404   1.323862940669414   4.850814300117984   8.094945654437128  10.000000000000000
                               -10.000000000000000  -8.280788099184663  -5.597215515050078  -2.748810950851591  -0.046973154443036   2.628746478697315   5.580455770857670   8.308546129122391  10.000000000000000
                               -10.000000000000000  -8.163534338901464  -4.927096812304078  -1.367267345989732   1.527426288656270   3.741019449075829   6.062638447109529   8.534948647172902  10.000000000000000
                               -10.000000000000000  -8.109245441054746  -4.512144839444513  -0.624446722716012   2.412722558496163   4.284891558562032   6.239667615039809   8.649971086156391  10.000000000000000];
expSolCPPenalty{2}(:, :, 2) = [ 10.000000000000000   7.991439381480221   4.417021629154303   0.628736753957972  -2.406903665583760  -4.743506998721412  -6.097805117888148  -8.485792391222926 -10.000000000000000
                                10.000000000000000   8.094945654436856   4.850814300117475   1.323862940669332  -1.603017925308194  -4.012273034120030  -5.989013583866205  -8.413524951790023 -10.000000000000000
                                10.000000000000000   8.308546129122229   5.580455770857428   2.628746478697277  -0.046973154442951  -2.748810950851608  -5.597215515050441  -8.280788099184933 -10.000000000000000
                                10.000000000000000   8.534948647172849   6.062638447109457   3.741019449075892   1.527426288656326  -1.367267345989716  -4.927096812304239  -8.163534338901574 -10.000000000000000
                                10.000000000000000   8.649971086156395   6.239667615039832   4.284891558562094   2.412722558496249  -0.624446722716014  -4.512144839444593  -8.109245441054782 -10.000000000000000];
expSolCPPenalty{2}(:, :, 3) = [ 10.000000000000000   8.669847235397885   6.538863359206387   5.044326941002879   4.510641974631483   5.044326941002738   6.538863359206973   8.669847235398144  10.000000000000000
                                10.000000000000000   8.681322725483943   6.632636605746702   5.187440994794531   4.754154527014120   5.187440994794462   6.632636605747138   8.681322725484163  10.000000000000000
                                10.000000000000000   8.704987414619296   6.736654654414263   5.324601625710137   4.965358874845513   5.324601625710116   6.736654654414476   8.704987414619429  10.000000000000000
                                10.000000000000000   8.725602667187237   6.699049746564246   5.161380419716784   4.771647973630958   5.161380419716774   6.699049746564318   8.725602667187282  10.000000000000000
                                10.000000000000000   8.734287250593610   6.646141125180145   4.982292188433463   4.536032661898719   4.982292188433473   6.646141125180155   8.734287250593612  10.000000000000000];
expSolCPPenalty{2}(:, :, 4) = [1     1     1     1     1     1     1     1     1
                               1     1     1     1     1     1     1     1     1
                               1     1     1     1     1     1     1     1     1
                               1     1     1     1     1     1     1     1     1
                               1     1     1     1     1     1     1     1     1];
expSolCPPenalty{3}(:, :, 1) = [ -10.000000000000000  -4.137529201659706   2.661235785498061   5.987768191144670  10.000000000000000
                                -6.651215150241988  -1.917123516963026   3.726498646495352   5.665788415273301   7.645724181068700
                                 0.099061524915535   2.815921522446965   5.807462639930961   6.312337087256048   6.397180384272028
                                 6.742157023858908   7.539554902788752   8.289379575090978   8.151539347950855   7.685718697708602
                                10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000];
expSolCPPenalty{3}(:, :, 2) = [  10.000000000000000   5.987768191144692   2.661235785498165  -4.137529201659809 -10.000000000000000
                                 7.645724181068673   5.665788415273289   3.726498646495411  -1.917123516963053  -6.651215150242018
                                 6.397180384272025   6.312337087256041   5.807462639930969   2.815921522446957   0.099061524915525
                                 7.685718697708601   8.151539347950850   8.289379575090974   7.539554902788752   6.742157023858906
                                 10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000];
expSolCPPenalty{3}(:, :, 3) = [  10.000000000000000   6.202281519000858   4.163607923688498   6.202281519000867  10.000000000000000
                                 7.878618705164055   5.459437789109303   4.075275891588737   5.459437789109312   7.878618705164069
                                 4.968956970908473   3.901740849025484   3.074438008650677   3.901740849025488   4.968956970908478
                                 2.045951408350618   1.601777847127133   1.261940803961638   1.601777847127130   2.045951408350620
                                                 0                   0                   0                   0                   0];
expSolCPPenalty{3}(:, :, 4) = [1     1     1     1     1
                               1     1     1     1     1
                               1     1     1     1     1
                               1     1     1     1     1
                               1     1     1     1     1];
expSolResHistoryPenalty = [   14.568073317508654
                              10.917362975892994
                               2.371581354888711
                               0.851633950048518
                               0.394919790321398
                               0.204089016843095
                               0.112264696861901
                               0.064257944172731
                               0.037729415349969
                               0.022520018481145];
                           
% Define the expected solution in terms of the Lagrange Multipliers 
% solution
expSolCPLM{1}(:, :, 1) = [ -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000
                           -8.229816514101145  -8.499956207168049  -8.699784531996727  -8.530221078282365  -8.001319974292969  -7.558718220978656
                           -6.722617062454779  -6.852797159960507  -6.859449421866599  -6.049986115977395  -4.104618993358851  -2.621679347386112
                           -6.679071757740918  -6.206270872378422  -5.535752905414360  -3.767384902691676  -0.203705819970992   2.335803183330495
                           -8.135806831658718  -6.179936754346971  -4.674600085648520  -1.433988326185997   3.603391492135026   7.397518516927594
                           -10.000000000000000  -6.454107478425250  -4.383720004854715  -0.227320032551111   5.312411068541386  10.000000000000000];
expSolCPLM{1}(:, :, 2) = [  -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000
                            -7.558718220978673  -8.001319974292985  -8.530221078282377  -8.699784531996714  -8.499956207168040  -8.229816514101126
                            -2.621679347386267  -4.104618993358955  -6.049986115977424  -6.859449421866622  -6.852797159960474  -6.722617062454705
                             2.335803183330334  -0.203705819971328  -3.767384902691869  -5.535752905414363  -6.206270872378259  -6.679071757740869
                             7.397518516927772   3.603391492134923  -1.433988326186929  -4.674600085648060  -6.179936754346877  -8.135806831658716
                          10.000000000000000   5.312411068541957  -0.227320032552679  -4.383720004853866  -6.454107478425363 -10.000000000000000];
expSolCPLM{1}(:, :, 3) = [   0                   0                   0                   0                   0                   0
                             1.550041912442032   1.290946944332028   1.019217700203102   1.019217700203114   1.290946944332043   1.550041912442051
                             3.905936756665364   3.301548199638172   2.603679703488517   2.603679703488520   3.301548199638250   3.905936756665495
                             5.973392702279196   4.869638020869913   3.753955727874678   3.753955727874794   4.869638020870203   5.973392702279352
                             8.353558787630563   6.230351530242373   4.493898273934236   4.493898273934940   6.230351530242475   8.353558787630405
                            10.000000000000000   6.920778053036258   4.712476493823951   4.712476493825128   6.920778053035813  10.000000000000000];
expSolCPLM{1}(:, :, 4) = [ 1     1     1     1     1     1
                           1     1     1     1     1     1
                           1     1     1     1     1     1
                           1     1     1     1     1     1
                           1     1     1     1     1     1
                           1     1     1     1     1     1];
expSolCPLM{2}(:, :, 1) = [ -10.000000000000000  -8.434030104885492  -5.963779557990398  -4.818671887630370  -2.423307238822679   0.688629828092513   4.299852605196446   7.920113632889332  10.000000000000000
                           -10.000000000000000  -8.343929516842378  -5.904589827998503  -3.994742399694924  -1.605815974081718   1.367024528913931   4.722137562547659   8.033856069692034  10.000000000000000
                           -10.000000000000000  -8.175594912215278  -5.579928801748433  -2.612361340130794  -0.025003288415208   2.677961431754960   5.397867377522750   8.273851360457735  10.000000000000000
                           -10.000000000000000  -7.998286646572189  -5.027632358821168  -0.987773213977366   1.491822656640681   3.939955578896182   5.693183997928230   8.556788280012110  10.000000000000000
                           -10.000000000000000  -7.897454480918208  -4.740948762148673   0.081980704470249   2.053804093994895   4.763614066230965   5.691424049217800   8.716252151243582  10.000000000000000];
expSolCPLM{2}(:, :, 2) = [  10.000000000000000   7.920113632889513   4.299852605196698   0.688629828091192  -2.423307238823381  -4.818671887629037  -5.963779557990669  -8.434030104885490 -10.000000000000000
                            10.000000000000000   8.033856069692188   4.722137562547822   1.367024528913023  -1.605815974082139  -3.994742399694147  -5.904589827998598  -8.343929516842413 -10.000000000000000
                            10.000000000000000   8.273851360457845   5.397867377522781   2.677961431754506  -0.025003288415381  -2.612361340130549  -5.579928801748280  -8.175594912215377 -10.000000000000000
                            10.000000000000000   8.556788280012096   5.693183997928355   3.939955578895754   1.491822656640797  -0.987773213977518  -5.027632358820777  -7.998286646572345 -10.000000000000000
                            10.000000000000000   8.716252151243484   5.691424049218068   4.763614066230304   2.053804093995574   0.081980704469646  -4.740948762148068  -7.897454480918415 -10.000000000000000];
expSolCPLM{2}(:, :, 3) = [  10.000000000000000   8.633887420028547   6.466843726316975   5.051971761072219   4.540886274336909   5.051971761073376   6.466843726316627   8.633887420028389  10.000000000000000
                            10.000000000000000   8.643950951328154   6.564450105806188   5.193424395912890   4.776521217466827   5.193424395913562   6.564450105805975   8.643950951328005  10.000000000000000
                            10.000000000000000   8.665244150653779   6.663207159590286   5.323046440487897   4.990783255692106   5.323046440488121   6.663207159590228   8.665244150653674  10.000000000000000
                            10.000000000000000   8.679311527882824   6.612286231010143   5.126005464083979   4.779604709332470   5.126005464084057   6.612286231010069   8.679311527882815  10.000000000000000
                            10.000000000000000   8.680406473308109   6.564342339074408   4.875114374260678   4.600484692869370   4.875114374260766   6.564342339074265   8.680406473308162  10.000000000000000];
expSolCPLM{2}(:, :, 4) = [ 1     1     1     1     1     1     1     1     1
                           1     1     1     1     1     1     1     1     1
                           1     1     1     1     1     1     1     1     1
                           1     1     1     1     1     1     1     1     1
                           1     1     1     1     1     1     1     1     1];
expSolCPLM{3}(:, :, 1) = [  -10.000000000000000  -4.011756841188636   2.904099230514266   5.672304562708411  10.000000000000000
                            -6.584683454008554  -1.776295276631318   3.855974828846238   5.591391127051978   7.592819875226418
                             0.190172407765815   2.894687702786316   5.847592491317050   6.334919867363173   6.395998586205607
                             6.770850402196975   7.558165206973351   8.301458510583931   8.163267126642348   7.702131158209299
                            10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000];
expSolCPLM{3}(:, :, 2) = [  10.000000000000000   5.672304562708327   2.904099230514235  -4.011756841188550 -10.000000000000000
                             7.592819875226414   5.591391127051921   3.855974828846221  -1.776295276631262  -6.584683454008533
                             6.395998586205587   6.334919867363165   5.847592491317048   2.894687702786312   0.190172407765828
                             7.702131158209292   8.163267126642351   8.301458510583933   7.558165206973361   6.770850402196973
                            10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000];
expSolCPLM{3}(:, :, 3) = [  10.000000000000000   6.052093312451190   4.202668343511523   6.052093312451201  10.000000000000000
                             7.846169072452724   5.400058400407203   4.065099642145016   5.400058400407198   7.846169072452720
                             4.930987065929906   3.870025870503435   3.049583114836333   3.870025870503441   4.930987065929910
                             2.030867417561272   1.590316426257705   1.252511770875857   1.590316426257703   2.030867417561272
                                           0                   0                   0                   0                   0];
expSolCPLM{3}(:, :, 4) = [ 1     1     1     1     1
                           1     1     1     1     1
                           1     1     1     1     1
                           1     1     1     1     1
                           1     1     1     1     1];
expSolResHistoryLM = [14.942294121403583
                      11.176225829326581
                       2.436235703459563
                       0.866971306477884
                       0.395877021938276
                       0.201395036930236
                       0.108835247678236
                       0.060985205670314
                       0.034890256353909
                       0.020169652695245];

% % Define the expected solution in terms of the Mortar solution
% expSolCPMortar{1}(:, :, 1) = [ -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000
%                                 -8.235130991262475  -8.503498417506318  -8.702805934398189  -8.535225420999550  -8.007634255328425  -7.566550612025906
%                                 -6.727855232279692  -6.862037256120201  -6.870270315340060  -6.062035055450008  -4.127642811861033  -2.650581774134907
%                                 -6.672083129035222  -6.177974811745685  -5.531925808909480  -3.808871751837193  -0.262322066987660   2.314256332967476
%                                 -8.150878255874144  -6.164670726847181  -4.535823511665825  -1.626582666347840   3.591157644049682   7.432229644526650
%                                -10.000000000000000  -6.565580334951189  -4.092623998014281  -0.585068314348617   5.472662571517183  10.000000000000000];
% expSolCPMortar{1}(:, :, 2) = [ -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000
%                                 -7.566550612025909  -8.007634255328421  -8.535225420999550  -8.702805934398189  -8.503498417506316  -8.235130991262473
%                                 -2.650581774134902  -4.127642811861030  -6.062035055450004  -6.870270315340056  -6.862037256120197  -6.727855232279692
%                                  2.314256332967475  -0.262322066987651  -3.808871751837195  -5.531925808909470  -6.177974811745693  -6.672083129035213
%                                  7.432229644526635   3.591157644049676  -1.626582666347833  -4.535823511665827  -6.164670726847198  -8.150878255874138
%                                 10.000000000000000   5.472662571517186  -0.585068314348617  -4.092623998014281  -6.565580334951219 -10.000000000000000];
% expSolCPMortar{1}(:, :, 3) = [   0                   0                   0                   0                   0                   0
%                                  1.544914903855130   1.287360729102798   1.016269702116835   1.016269702116835   1.287360729102797   1.544914903855131
%                                  3.893810261933955   3.288837441606556   2.594871554106681   2.594871554106681   3.288837441606553   3.893810261933953
%                                  5.971842927130913   4.854994713680650   3.733838299618284   3.733838299618283   4.854994713680648   5.971842927130910
%                                  8.358403184371870   6.233100548063092   4.450099874571729   4.450099874571723   6.233100548063097   8.358403184371872
%                                 10.000000000000000   6.973987172591881   4.644764509976969   4.644764509976962   6.973987172591893  10.000000000000000];
% expSolCPMortar{1}(:, :, 4) = [ 1     1     1     1     1     1
%                                1     1     1     1     1     1
%                                1     1     1     1     1     1
%                                1     1     1     1     1     1
%                                1     1     1     1     1     1
%                                1     1     1     1     1     1];
% expSolCPMortar{2}(:, :, 1) = [ -10.000000000000000  -8.476961253707083  -6.046623350620751  -4.535032998043901  -2.455591211184007   0.343450745658256   4.429808182056598   7.976494288977194  10.000000000000000
%                                -10.000000000000000  -8.394176242373097  -5.918344472664759  -3.858707505270488  -1.628458703894312   1.155274805419928   4.816545899455973   8.082365945938943  10.000000000000000
%                                -10.000000000000000  -8.237141871148578  -5.499425733647901  -2.631039895950564  -0.014876960519225   2.559753223100252   5.489130419566362   8.294207747985670  10.000000000000000
%                                -10.000000000000000  -8.091755136971832  -4.807932340218783  -1.211524402629854   1.627691475340489   3.700302197363143   5.912369918092145   8.518267963984339  10.000000000000000
%                                -10.000000000000000  -8.023465371513298  -4.380955176246090  -0.441606495107560   2.575417144755306   4.262731682331458   6.047811705586636   8.633842270398128  10.000000000000000];
% expSolCPMortar{2}(:, :, 2) = [  10.000000000000000   7.976494288977197   4.429808182056598   0.343450745658261  -2.455591211184012  -4.535032998043901  -6.046623350620782  -8.476961253707097 -10.000000000000000
%                                 10.000000000000000   8.082365945938946   4.816545899455961   1.155274805419948  -1.628458703894330  -3.858707505270476  -5.918344472664793  -8.394176242373106 -10.000000000000000
%                                 10.000000000000000   8.294207747985666   5.489130419566346   2.559753223100257  -0.014876960519237  -2.631039895950547  -5.499425733647925  -8.237141871148587 -10.000000000000000
%                                 10.000000000000000   8.518267963984343   5.912369918092109   3.700302197363155   1.627691475340473  -1.211524402629833  -4.807932340218793  -8.091755136971829 -10.000000000000000
%                                 10.000000000000000   8.633842270398123   6.047811705586602   4.262731682331458   2.575417144755310  -0.441606495107552  -4.380955176246087  -8.023465371513296 -10.000000000000000];
% expSolCPMortar{2}(:, :, 3) = [  10.000000000000000   8.659223121960736   6.503515933639705   5.002408698488646   4.463267806887126   5.002408698488643   6.503515933639715   8.659223121960741  10.000000000000000
%                                 10.000000000000000   8.667930784783794   6.583460612191028   5.164567910311030   4.749546560057021   5.164567910311026   6.583460612191040   8.667930784783795  10.000000000000000
%                                 10.000000000000000   8.683668613902428   6.670215616940223   5.309884606644733   4.984842364458477   5.309884606644732   6.670215616940236   8.683668613902430  10.000000000000000
%                                 10.000000000000000   8.695967007242254   6.611769905973011   5.126030622468861   4.786343324014972   5.126030622468858   6.611769905973027   8.695967007242251  10.000000000000000
%                                 10.000000000000000   8.700806027586173   6.545848203551560   4.924728537260633   4.546145605626920   4.924728537260631   6.545848203551578   8.700806027586170  10.000000000000000];
% expSolCPMortar{2}(:, :, 4) = [ 1     1     1     1     1     1     1     1     1
%                                1     1     1     1     1     1     1     1     1
%                                1     1     1     1     1     1     1     1     1
%                                1     1     1     1     1     1     1     1     1
%                                1     1     1     1     1     1     1     1     1];
% expSolCPMortar{3}(:, :, 1) = [ -10.000000000000000  -4.000131013864399   2.833252539602594   5.712707432274881  10.000000000000000
%                                 -6.586935168247373  -1.789585837639175   3.818842173890346   5.595990389497777   7.600185479411128
%                                  0.178297202858219   2.881659563737727   5.840607269836043   6.330007852263861   6.396502556544695
%                                  6.766341565724941   7.555742547751741   8.300015317455607   8.161687914834827   7.699542258437274
%                                 10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000];
% expSolCPMortar{3}(:, :, 2) = [  10.000000000000000   5.712707432274851   2.833252539602608  -4.000131013864393 -10.000000000000000
%                                  7.600185479411143   5.595990389497771   3.818842173890359  -1.789585837639190  -6.586935168247416
%                                  6.396502556544696   6.330007852263862   5.840607269836037   2.881659563737737   0.178297202858220
%                                  7.699542258437273   8.161687914834822   8.300015317455609   7.555742547751731   6.766341565724943
%                                 10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000];
% expSolCPMortar{3}(:, :, 3) = [  10.000000000000000   6.071057196528724   4.195669612337601   6.071057196528728  10.000000000000000
%                                  7.852932194503535   5.407460112193630   4.068546117312629   5.407460112193632   7.852932194503539
%                                  4.936623139303975   3.874806258649678   3.053833761354053   3.874806258649678   4.936623139303976
%                                  2.033098872331578   1.592017696195557   1.253719205251356   1.592017696195558   2.033098872331577
%                                                0                   0                   0                   0                   0];
% expSolCPMortar{3}(:, :, 4) = [ 1     1     1     1     1
%                                1     1     1     1     1
%                                1     1     1     1     1
%                                1     1     1     1     1
%                                1     1     1     1     1];
% expSolResHistoryMortar = [14.986606440334855
%                           11.186416933552758
%                            2.455039838307761
%                            0.870927352076550
%                            0.395368939735741
%                            0.200458869000483
%                            0.108268189826596
%                            0.060799889742575
%                            0.034969337268159
%                            0.020402225036517];
                       
% Define the expected solution in terms of the Nitsche solution
expSolCPNitsche{1}(:, :, 1) = [ -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000
                                -8.232754767529672  -8.501617429761888  -8.701133115602126  -8.532894548429255  -8.004530205714405  -7.562800672915575
                                -6.725493058234185  -6.858810381961397  -6.866123859568382  -6.056660593293731  -4.118040820161285  -2.638080598760276
                                -6.674853522883785  -6.185836169870142  -5.533718346627689  -3.797111099112687  -0.239872313089177   2.331440879137696
                                -8.154454714390008  -6.180600532652592  -4.575406128648621  -1.580444918972558   3.616891493257725   7.437038351051294
                               -10.000000000000000  -6.562040666604913  -4.175711411277876  -0.496528249193025   5.468422852932481  10.000000000000000];
expSolCPNitsche{1}(:, :, 2) = [ -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000 -10.000000000000000
                                -7.562800672915584  -8.004530205714413  -8.532894548429256  -8.701133115602127  -8.501617429761886  -8.232754767529666
                                -2.638080598760302  -4.118040820161307  -6.056660593293741  -6.866123859568391  -6.858810381961402  -6.725493058234181
                                 2.331440879137646  -0.239872313089221  -3.797111099112723  -5.533718346627694  -6.185836169870154  -6.674853522883790
                                 7.437038351051262   3.616891493257650  -1.580444918972623  -4.575406128648656  -6.180600532652644  -8.154454714390040
                                10.000000000000000   5.468422852932368  -0.496528249193096  -4.175711411277954  -6.562040666604974 -10.000000000000000];
expSolCPNitsche{1}(:, :, 3) = [  0                   0                   0                   0                   0                   0
                                 1.547278130214763   1.289240351845990   1.017788062865129   1.017788062865129   1.289240351845993   1.547278130214768
                                 3.899351169031767   3.293789602160240   2.598487481085360   2.598487481085362   3.293789602160249   3.899351169031778
                                 5.977346323984849   4.861934272259566   3.738495358625831   3.738495358625840   4.861934272259580   5.977346323984867
                                 8.364141444763034   6.244427117420158   4.458794407597715   4.458794407597723   6.244427117420203   8.364141444763058
                                10.000000000000000   6.981428942172371   4.660148878334259   4.660148878334277   6.981428942172435  10.000000000000000];
expSolCPNitsche{1}(:, :, 4) = [1     1     1     1     1     1
                               1     1     1     1     1     1
                               1     1     1     1     1     1
                               1     1     1     1     1     1
                               1     1     1     1     1     1
                               1     1     1     1     1     1];
expSolCPNitsche{2}(:, :, 1) = [-10.000000000000000  -8.478547467310403  -6.047442386523781  -4.619139547239155  -2.452896720666154   0.432886642536577   4.429391763965686   7.978730647247638  10.000000000000000
                               -10.000000000000000  -8.394336828190113  -5.935303123245781  -3.918546594896622  -1.630576481993583   1.206397532934877   4.835879243398063   8.085064988182236  10.000000000000000
                               -10.000000000000000  -8.251308029789307  -5.525293728262644  -2.680989078088968  -0.037147849162087   2.592905241565989   5.517182443528941   8.314119742916560  10.000000000000000
                               -10.000000000000000  -8.092770404970334  -4.867107153883556  -1.269063520210062   1.592205012553793   3.697156483150698   6.009472317085986   8.512363403574627  10.000000000000000
                               -10.000000000000000  -8.041084395499977  -4.433158022754155  -0.514173601456390   2.510907248669338   4.270961254735414   6.133068702624179   8.654174991018625  10.000000000000000];
expSolCPNitsche{2}(:, :, 2) = [ 10.000000000000000   7.978730647247588   4.429391763965572   0.432886642536501  -2.452896720666225  -4.619139547239231  -6.047442386523848  -8.478547467310431 -10.000000000000000
                                10.000000000000000   8.085064988182198   4.835879243397914   1.206397532934784  -1.630576481993684  -3.918546594896684  -5.935303123245863  -8.394336828190131 -10.000000000000000
                                10.000000000000000   8.314119742916477   5.517182443528875   2.592905241565907  -0.037147849162239  -2.680989078089041  -5.525293728262695  -8.251308029789344 -10.000000000000000
                                10.000000000000000   8.512363403574630   6.009472317085700   3.697156483150454   1.592205012553654  -1.269063520210173  -4.867107153883675  -8.092770404970329 -10.000000000000000
                                10.000000000000000   8.654174991018555   6.133068702624017   4.270961254735210   2.510907248669164  -0.514173601456524  -4.433158022754217  -8.041084395500004 -10.000000000000000];
expSolCPNitsche{2}(:, :, 3) = [ 10.000000000000000   8.664016704310619   6.512056156504556   5.016760733834547   4.479545269423578   5.016760733834574   6.512056156504615   8.664016704310651  10.000000000000000
                                10.000000000000000   8.672513530490036   6.600245867911538   5.172659429791882   4.752257870464942   5.172659429791913   6.600245867911612   8.672513530490061  10.000000000000000
                                10.000000000000000   8.700370537287416   6.687480188342523   5.321782231431345   4.977568750452381   5.321782231431389   6.687480188342556   8.700370537287462  10.000000000000000
                                10.000000000000000   8.698798237028601   6.658409217260806   5.135581326879664   4.790030276688952   5.135581326879712   6.658409217260948   8.698798237028601  10.000000000000000
                                10.000000000000000   8.719520564626549   6.582880782105378   4.951140579145116   4.547777294296350   4.951140579145181   6.582880782105453   8.719520564626588  10.000000000000000];
expSolCPNitsche{2}(:, :, 4) = [1     1     1     1     1     1     1     1     1
                               1     1     1     1     1     1     1     1     1
                               1     1     1     1     1     1     1     1     1
                               1     1     1     1     1     1     1     1     1
                               1     1     1     1     1     1     1     1     1];
expSolCPNitsche{3}(:, :, 1) = [ -10.000000000000000  -4.059342039207384   2.755225951939460   5.812411065524540  10.000000000000000
                                -6.619095825470098  -1.847135168617869   3.775895761899217   5.618728953162550   7.617545846841324
                                 0.142289767787962   2.852302809574366   5.825049697640672   6.320480316910356   6.394815800506212
                                 6.755378919349129   7.548238972869966   8.294973382275110   8.156566543432946   7.692591589942994
                                10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000];
expSolCPNitsche{3}(:, :, 2) = [ 10.000000000000000   5.812411065524343   2.755225951939271  -4.059342039207450 -10.000000000000000
                                 7.617545846841311   5.618728953162511   3.775895761899088  -1.847135168617986  -6.619095825470152
                                 6.394815800506210   6.320480316910327   5.825049697640653   2.852302809574325   0.142289767787908
                                 7.692591589942997   8.156566543432968   8.294973382275117   7.548238972869965   6.755378919349121
                                10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000  10.000000000000000];
expSolCPNitsche{3}(:, :, 3) = [ 10.000000000000000   6.119747631714021   4.193704173673707   6.119747631714124  10.000000000000000
                                 7.865564409313001   5.430447104342821   4.078128848254451   5.430447104342876   7.865564409313016
                                 4.951793322175211   3.887883899026013   3.064200803538159   3.887883899026030   4.951793322175235
                                 2.039396307163653   1.596769038103244   1.257587295578011   1.596769038103251   2.039396307163658
                                               0                   0                   0                   0                   0];
expSolCPNitsche{3}(:, :, 4) = [1     1     1     1     1
                               1     1     1     1     1
                               1     1     1     1     1
                               1     1     1     1     1
                               1     1     1     1     1];
expSolResHistoryNitsche = [   14.870644980804045
                              11.120624155568400
                               2.431007845632684
                               0.858957517106702
                               0.389632725400186
                               0.197584568640394
                               0.107195891000782
                               0.060874456265265
                               0.035619930257437
                               0.021264247162000];
expSolAutomaticStabilization = 1.0e+06 *[0.216583389266793   0.333487588358564   0.458880826279685   0.579478929297393   0.685937881215083   0.772252341097829   0.837163003095889   0.883136793136238   0.795729458477417   0.815276510748167
                                         0.194743365955183   0.301190765460125   0.452721343332373   0.636070337708309   0.823971213911730   0.978241305322404   1.081630224258979   1.142966731250184   1.177841012889437   1.197656911313093];

%% 9. Define the form-finding analysis properties
propFormFinding.tolerance = 1e-6;
propFormFinding.maxNoIter = 10;
propFormFinding.minNoIter = 1;

%% 10. Solve the multipatch problem with Penalty

% Initialize the B-Spline patch array for the penalty method
BSplinePatchesPenalty = BSplinePatches;

% Assign the parameters for multipatch coupling
% ---------------------------------------------

propCouplingPenalty.method = 'penalty';
propCouplingPenalty.intC = intC;
propCouplingPenalty.alphaD = zeros(connections.No, 1);
propCouplingPenalty.alphaR = zeros(connections.No, 1);
for iConnections = 1:connections.No
    % Get the id's of the patches
    IDPatchI = connections.xiEtaCoup(iConnections, 1);
    IDPatchJ = connections.xiEtaCoup(iConnections, 2);

    % Get the mean polynomial order between the patches
    isOnXiI = false;
    if connections.xiEtaCoup(iConnections, 5) == connections.xiEtaCoup(iConnections, 6)
        isOnXiI = true;
    end
    if isOnXiI
        polOrderI = BSplinePatchesPenalty{IDPatchI}.p;
    else
        polOrderI = BSplinePatchesPenalty{IDPatchI}.q;
    end
    isOnXiJ = false;
    if connections.xiEtaCoup(iConnections, 9) == connections.xiEtaCoup(iConnections, 10)
        isOnXiJ = true;
    end
    if isOnXiJ
        polOrderJ = BSplinePatchesPenalty{IDPatchJ}.p;
    else
        polOrderJ = BSplinePatchesPenalty{IDPatchJ}.q;
    end
    polOrderMean = mean([polOrderI polOrderJ]);

    % Assign the penalty parameters
    propCouplingPenalty.alphaD(iConnections,1) = ...
        max([norm(eval(['Dm' num2str(IDPatchI)])) ...
        norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
        min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea]);
    propCouplingPenalty.alphaR(iConnections, 1) = 0;
end

% Solve the problem using the Penalty method
% ------------------------------------------

[BSplinePatchesPenalty, ~, ~, resHistoryPenalty, ~, ~] = ...
    solve_DDMFormFindingIGAMembrane...
    (BSplinePatchesPenalty, connections, propCouplingPenalty, ...
    propFormFinding, solve_LinearSystem, '');

%% 11. Solve the multipatch problem with Lagrange Multipliers

% Initialize the B-Spline patch array for the penalty method
BSplinePatchesLM = BSplinePatches;

% Assign the parameters for multipatch coupling
% ---------------------------------------------

propCouplingLM.method = 'lagrangeMultipliers';
propCouplingLM.alphaD = zeros(connections.No, 1);
propCouplingLM.alphaR = zeros(connections.No, 1);
propCouplingLM.intC = intC;

% Solve the problem using the Lagrange Multipliers method
% -------------------------------------------------------

[BSplinePatchesLM, ~, ~, resHistoryLM, ~, ~] = ...
    solve_DDMFormFindingIGAMembrane ...
    (BSplinePatchesLM, connectionsLM, propCouplingLM, propFormFinding, ...
    solve_LinearSystem, '');

%% 12. Solve the multipatch problem with Mortar

% % Initialize the B-Spline patch array for the penalty method
% BSplinePatchesMortar = BSplinePatches;
% 
% % Assign the parameters for multipatch coupling
% % ---------------------------------------------
% 
% propCouplingMortar.method = 'mortar';
% propCouplingMortar.isSlaveSideCoarser = false;
% propCouplingMortar.computeRearrangedProblemMtrcs = @computeRearrangedProblemMtrcs4MortarIGAMembrane;
% propCouplingMortar.intC = intC;
% 
% % Solve the problem using the mortar method
% % -------------------------------------------------------
% 
% [BSplinePatchesMortar,~,~,resHistoryMortar,~,~] = ...
%     solve_DDMFormFindingIGAMembrane ...
%     (BSplinePatchesMortar, connections, propCouplingMortar, ...
%     propFormFinding, solve_LinearSystem, '');

%% 13. Solve the multipatch problem with Nitsche

% Initialize the B-Spline patch array for the penalty method
BSplinePatchesNitsche = BSplinePatches;

% Assign the parameters for multipatch coupling
% ---------------------------------------------

propCouplingNitsche.method = 'nitsche';
propCouplingNitsche.estimationStabilPrm = true;
propCouplingNitsche.gammaTilde = .5;
propCouplingNitsche.intC = intC;

% Solve the problem using the Penalty method
% ------------------------------------------

[BSplinePatchesNitsche, ~, propCouplingNitsche, resHistoryNitsche, ~, ~] = ...
    solve_DDMFormFindingIGAMembrane ...
    (BSplinePatchesNitsche, connections, propCouplingNitsche, ...
    propFormFinding, solve_LinearSystem, '');

%% 14. Verify the solution
for iPatches = 1:noPatches
    testCase.verifyEqual(BSplinePatchesPenalty{iPatches}.CP, expSolCPPenalty{iPatches}, 'AbsTol', absTolDisp*1e3);
    testCase.verifyEqual(BSplinePatchesLM{iPatches}.CP,expSolCPLM{iPatches}, 'AbsTol', absTolDisp);
%     testCase.verifyEqual(BSplinePatchesMortar{iPatches}.CP,expSolCPMortar{iPatches}, 'AbsTol', absTolDisp*1e1);
    testCase.verifyEqual(BSplinePatchesNitsche{iPatches}.CP, expSolCPNitsche{iPatches}, 'AbsTol', absTolDisp);
end
testCase.verifyEqual(resHistoryPenalty, expSolResHistoryPenalty, 'AbsTol', absTolResHistory*1e2);
testCase.verifyEqual(resHistoryLM, expSolResHistoryLM, 'AbsTol', absTolResHistory);
% testCase.verifyEqual(resHistoryMortar, expSolResHistoryMortar, 'AbsTol', absTolResHistory*1e1);
testCase.verifyEqual(resHistoryNitsche, expSolResHistoryNitsche, 'AbsTol', absTolResHistory);
testCase.verifyEqual(propCouplingNitsche.automaticStabilization, expSolAutomaticStabilization, 'AbsTol', absTolStabil);

end
