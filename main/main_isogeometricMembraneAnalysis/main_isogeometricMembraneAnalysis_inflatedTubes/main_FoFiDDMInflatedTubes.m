%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Form-finding for the torus structure.
%
% Date : 01.02.2018
%
%% Preamble
clear;
clc;

%% Includes 

% Add general math functions
addpath('../../../generalMath/');

% Add general auxiliary functions
addpath('../../../auxiliary/');

% Add system solvers
addpath('../../../equationSystemSolvers/');

% Add efficient computation functions
addpath('../../../efficientComputation/');

% Add transient analysis solvers
addpath('../../../transientAnalysis/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../../CAGDKernel/CAGDKernel_graphics/',...
        '../../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the isogeometric Kirchhoff-Love shell formulation
addpath('../../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../../isogeometricThinStructureAnalysis/graphicsMultipatches/',...
        '../../../isogeometricThinStructureAnalysis/loads/',...
        '../../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../../isogeometricThinStructureAnalysis/solvers/',...
        '../../../isogeometricThinStructureAnalysis/metrics/',...
        '../../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../../isogeometricThinStructureAnalysis/penaltyDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/penaltyDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../../isogeometricThinStructureAnalysis/output/',...
        '../../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/');
    
% Add functions related to the visualization of the multipatch geometry
% related to the isogeometric mortar-based mapping
addpath('../../../isogeometricMortarBasedMappingAnalysis/graphics/');

%% Define the geometry in terms of NURBS multipatches

% General parameters:
radius = 125e-2;
Radius = 1e1 + 125e-2;
arcLengthXi1 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta1 = 2*pi*radius/4;
arcLengthXi2 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta2 = 2*pi*radius/4;
arcLengthXi3 = 2*pi*Radius/4;
arcLengthEta3 = 2*pi*radius/4;
arcLengthXi4 = 2*pi*Radius/4;
arcLengthEta4 = 2*pi*radius/4;
arcLengthXi5 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta5 = 2*pi*radius/4;
arcLengthXi6 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta6 = 2*pi*radius/4;
arcLengthXi7 = 2*pi*Radius/4;
arcLengthEta7 = 2*pi*radius/4;
arcLengthXi8 = 2*pi*Radius/4;
arcLengthEta8 = 2*pi*radius/4;
arcLengthXi9 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta9 = 2*pi*radius/4;
arcLengthXi10 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta10 = 2*pi*radius/4;
arcLengthXi11 = 2*pi*Radius/4;
arcLengthEta11 = 2*pi*radius/4;
arcLengthXi12 = 2*pi*Radius/4;
arcLengthEta12 = 2*pi*radius/4;
arcLengthXi13 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta13 = 2*pi*radius/4;
arcLengthXi14 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta14 = 2*pi*radius/4;
arcLengthXi15 = 2*pi*Radius/4;
arcLengthEta15 = 2*pi*radius/4;
arcLengthXi16 = 2*pi*Radius/4;
arcLengthEta16 = 2*pi*radius/4;
arcLengthXi17 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta17 = 2*pi*radius/4;
arcLengthXi18 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta18 = 2*pi*radius/4;
arcLengthXi19 = 2*pi*Radius/4;
arcLengthEta19 = 2*pi*radius/4;
arcLengthXi20 = 2*pi*Radius/4;
arcLengthEta20 = 2*pi*radius/4;
arcLengthXi21 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta21 = 2*pi*radius/4;
arcLengthXi22 = 2*pi*(Radius + 2*radius)/4;
arcLengthEta22 = 2*pi*radius/4;
arcLengthXi23 = 2*pi*Radius/4;
arcLengthEta23 = 2*pi*radius/4;
arcLengthXi24 = 2*pi*Radius/4;
arcLengthEta24 = 2*pi*radius/4;

% Preassumed number of patches
noPatches = 8; % 24

% Define the patches of the multipatch geometry
for iPatches = 1
    % Patch 1 :
    % _________

    % Polynomial orders
    p1 = 2;
    q1 = 2;

    % Knot vectors
    Xi1 = [0 0 0 1 1 1];
    Eta1 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP1(:,:,1) = [Radius + radius Radius + radius Radius
                  Radius + radius Radius + radius Radius
                  0               0               0];

    % z-coordinates
    CP1(:,:,2) = [0 radius radius
                  0 radius radius
                  0 radius radius];

    % y-coordinates
    CP1(:,:,3) = [0               0               0
                  Radius + radius Radius + radius Radius
                  Radius + radius Radius + radius Radius];

    % weights
    CP1(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS1 = 0;
    nxi1 = length(CP1(:,1,1));
    neta1 = length(CP1(1,:,1));
    for i= 1:nxi1
        for j = 1:neta1
            if CP1(i,j,4)~=1
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

    % Polynomial orders
    p2 = 2;
    q2 = 2;

    % Knot vectors
    Xi2 = [0 0 0 1 1 1];
    Eta2 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP2(:,:,1) = [Radius + radius Radius + radius Radius
                  Radius + radius Radius + radius Radius
                  0               0               0];

    % y-coordinates
    CP2(:,:,2) = -[0 radius radius
                   0 radius radius
                   0 radius radius];

    % z-coordinates
    CP2(:,:,3) = [0               0               0
                  Radius + radius Radius + radius Radius
                  Radius + radius Radius + radius Radius];

    % weights
    CP2(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS2 = 0;
    nxi2 = length(CP2(:,1,1));
    neta2 = length(CP2(1,:,1));
    for i = 1:nxi2
        for j = 1:neta2
            if CP2(i,j,4)~=1
                isNURBS2 = 1;
                break;
            end
        end
        if isNURBS2
            break;
        end
    end

    % Patch 3 :
    % _________

    % Polynomial orders
    p3 = 2;
    q3 = 2;

    % Knot vectors
    Xi3 = [0 0 0 1 1 1];
    Eta3 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP3(:,:,1) = [Radius - radius Radius - radius Radius
                  Radius - radius Radius - radius Radius
                  0               0               0];

    % y-coordinates
    CP3(:,:,2) = -[0 radius radius
                   0 radius radius
                   0 radius radius];

    % z-coordinates
    CP3(:,:,3) = [0               0               0
                  Radius - radius Radius - radius Radius
                  Radius - radius Radius - radius Radius];

    % weights
    CP3(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS3 = 0;
    nxi3 = length(CP3(:,1,1));
    neta3 = length(CP3(1,:,1));
    for i = 1:nxi3
        for j = 1:neta3
            if CP3(i,j,4)~=1
                isNURBS3 = 1;
                break;
            end
        end
        if isNURBS3
            break;
        end
    end

    % Patch 4 :
    % _________

    % Polynomial orders
    p4 = 2;
    q4 = 2;

    % Knot vectors
    Xi4 = [0 0 0 1 1 1];
    Eta4 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP4(:,:,1) = [Radius - radius Radius - radius Radius
                  Radius - radius Radius - radius Radius
                  0               0               0];

    % y-coordinates
    CP4(:,:,2) = [0 radius radius
                  0 radius radius
                  0 radius radius];

    % z-coordinates
    CP4(:,:,3) = [0               0               0
                  Radius - radius Radius - radius Radius
                  Radius - radius Radius - radius Radius];

    % weights
    CP4(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS4 = 0;
    nxi4 = length(CP4(:,1,1));
    neta4 = length(CP4(1,:,1));
    for i = 1:nxi4
        for j = 1:neta4
            if CP4(i,j,4)~=1
                isNURBS4 = 1;
                break;
            end
        end
        if isNURBS4
            break;
        end
    end

    % Patch 5 :
    % _________

    % Polynomial orders
    p5 = 2;
    q5 = 2;

    % Knot vectors
    Xi5 = [0 0 0 1 1 1];
    Eta5 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP5(:,:,1) = -[Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius
                   0               0               0];

    % y-coordinates
    CP5(:,:,2) = [0 radius radius
                  0 radius radius
                  0 radius radius];

    % z-coordinates
    CP5(:,:,3) = [0               0               0
                  Radius + radius Radius + radius Radius
                  Radius + radius Radius + radius Radius];

    % weights
    CP5(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS5 = 0;
    nxi5 = length(CP5(:,1,1));
    neta5 = length(CP5(1,:,1));
    for i= 1:nxi5
        for j = 1:neta5
            if CP5(i,j,4)~=1
                isNURBS5 = 1;
                break;
            end
        end
        if isNURBS5
            break;
        end
    end

    % Patch 6 :
    % _________

    % Polynomial orders
    p6 = 2;
    q6 = 2;

    % Knot vectors
    Xi6 = [0 0 0 1 1 1];
    Eta6 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP6(:,:,1) = -[Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius
                   0               0               0];

    % y-coordinates
    CP6(:,:,2) = -[0 radius radius
                   0 radius radius
                   0 radius radius];

    % z-coordinates
    CP6(:,:,3) = [0               0               0
                  Radius + radius Radius + radius Radius
                  Radius + radius Radius + radius Radius];

    % weights
    CP6(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS6 = 0;
    nxi6 = length(CP6(:,1,1));
    neta6 = length(CP6(1,:,1));
    for i = 1:nxi6
        for j = 1:neta6
            if CP6(i,j,4)~=1
                isNURBS6 = 1;
                break;
            end
        end
        if isNURBS6
            break;
        end
    end

    % Patch 7 :
    % _________

    % Polynomial orders
    p7 = 2;
    q7 = 2;

    % Knot vectors
    Xi7 = [0 0 0 1 1 1];
    Eta7 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP7(:,:,1) = -[Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius
                   0               0               0];

    % y-coordinates
    CP7(:,:,2) = -[0 radius radius
                   0 radius radius
                   0 radius radius];

    % z-coordinates
    CP7(:,:,3) = [0               0               0
                  Radius - radius Radius - radius Radius
                  Radius - radius Radius - radius Radius];

    % weights
    CP7(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS7 = 0;
    nxi7 = length(CP7(:,1,1));
    neta7 = length(CP7(1,:,1));
    for i = 1:nxi7
        for j = 1:neta7
            if CP7(i,j,4)~=1
                isNURBS7 = 1;
                break;
            end
        end
        if isNURBS7
            break;
        end
    end

    % Patch 8 :
    % _________

    % Polynomial orders
    p8 = 2;
    q8 = 2;

    % Knot vectors
    Xi8 = [0 0 0 1 1 1];
    Eta8 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP8(:,:,1) = -[Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius
                   0               0               0];

    % y-coordinates
    CP8(:,:,2) = [0 radius radius
                  0 radius radius
                  0 radius radius];

    % z-coordinates
    CP8(:,:,3) = [0               0               0
                  Radius - radius Radius - radius Radius
                  Radius - radius Radius - radius Radius];

    % weights
    CP8(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS8 = 0;
    nxi8 = length(CP8(:,1,1));
    neta8 = length(CP8(1,:,1));
    for i = 1:nxi8
        for j = 1:neta8
            if CP8(i,j,4)~=1
                isNURBS8 = 1;
                break;
            end
        end
        if isNURBS8
            break;
        end
    end
end

%% %% Material constants

% Membrane :
% ----------

% Young's modulus
parameters.E = 1e9;

% Poisson ratio
parameters.nue = 0.0;

% Thickness
parameters.t = 1.0;

% Define inner pressure and second normal component of the stress
FAmp = 10.0; % 1e3
sigma0 = abs(FAmp)*radius/(2*parameters.t); % abs(FAmp)*radius/2/parameters.t;

% Prestress for the membrane
% parameters.prestress.computeParametricCoordinates = @(X) [atan(X(3,1)/X(1,1))
%                                                           asin(X(2,1)/radius)];

% parameters.prestress.computeParametricCoordinates = @(X) [atan(X(3,1)/X(1,1))
%                                                           real(acos((sqrt(X(1,1)^2 + X(3,1)^2) - Radius)/radius))];
                                                      
parameters.prestress.computeParametricCoordinates = @(X) [atan(X(3,1)/X(1,1))
                                                          real(acos((sqrt(X(1,1)^2 + X(3,1)^2) - Radius)/radius)) + pi/2];

% baseVctTheta = @(theta1,theta2) [-(Radius + radius*cos(theta2))*sin(theta1)
%                                  0
%                                  (Radius + radius*cos(theta2))*cos(theta1)];
% baseVctPhi = @(theta1,theta2) [(Radius - radius*sin(theta2))*cos(theta1)
%                                radius*cos(theta2)
%                                (Radius - radius*sin(theta2))*sin(theta1)];
% parameters.prestress.computeBaseVectors = @(theta1,theta2) [baseVctTheta(theta1,theta2) baseVctPhi(theta1,theta2)];

parameters.prestress.voigtVector = @(theta)sigma0*[1
                                                   (2*Radius + radius*sin(theta(2,1)))/(Radius + radius*sin(theta(2,1)))
                                                   0];

%% GUI

% Case name
caseName = 'DDMInflatedTubes';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Choose a method for the application of weak Dirichlet boundary conditions and the multipatch coupling
method = 'Nitsche';
if ~strcmp(method,'Penalty') && ~strcmp(method,'LagrangeMultipliers') && ...
        ~strcmp(method,'Mortar') && ~strcmp(method,'AugmentedLagrangeMultipliers') && ...
        ~strcmp(method,'Nitsche')
    error('%s is not a valid method (Nitsche, Penalty, LagrangeMultipliers, AugmentedLagrangeMultipliers)',method);
end

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points
int.type = 'default';
if strcmp(int.type,'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.etaNGPForLoad = 6;
    int.nGPForLoad = 6;
    int.nGPError = 12;
end

% Loop over the patches
for iPatches = 1:noPatches
    assignin('base',['int' num2str(iPatches)],int);
end

% Interface integration :
% _______________________

intC.type = 'user';
intC.method = 'Nitsche';
if strcmp(intC.type,'user')
    if strcmp(intC.method,'lagrangeMultipliers')
        intC.nGP1 = 16;
        intC.nGP2 = 16;
    else
        intC.noGPs = 16;
    end
    intC.nGPError = 16;
end

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'current';

% Plot strain or stress field
% .resultant: 'displacement','strain','curvature','force','moment','shearForce'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x','y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

% Define the coupling properties

% Material matrices for the patches :
% ___________________________________

Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
     [1              parameters.nue 0
      parameters.nue 1              0
      0               0              (1-parameters.nue)/2];
Db = parameters.t^2/12*Dm;
for iPatches = 1:noPatches
    assignin('base',['Dm' num2str(iPatches)],Dm);
    assignin('base',['Db' num2str(iPatches)],Db);
end

% Assign the penalty factors

% Function handle to writing out the results
% writeOutput = @writeResults4Carat;
% writeOutput = @writeResults4GiD;
writeOutput = 'undefined';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

%% Plot the prestress value along phi
% dphi = (2*pi - 0)/100;
% phi = zeros(101,1);
% prestresses = zeros(101,1);
% f = @(phi) (2*Radius + radius*sin(phi))/(Radius + radius*sin(phi));
% for iPhi = 1:101
%     phi(iPhi,1) = (iPhi - 1)*dphi;
%     ps = parameters.prestress.voigtVector([0;phi(iPhi)]);
%     prestresses(iPhi,1) = ps(2,1);
% end
% figure(graph.index)
% plot(phi,prestresses);
% graph.index = graph.index + 1;

%% Refinement

% Mesh size
meshSize = 'fine'; % 'coarsest', 'coarse' % 'dummy', 'coarse', 'fine'
if ~strcmp(meshSize,'dummy') && ~strcmp(meshSize,'coarse') && ~strcmp(meshSize,'coarsest') && ~strcmp(meshSize,'fine')
   error('Choose between a dummy, a coarse, a coarsest and a fine setting of the isogeometric discretization');
end

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%
for iPRef = 1
    % Patch 1 :
    % _________

    if strcmp(meshSize,'coarse') || strcmp(meshSize,'coarsest')
        a1 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a1 = 0;
    elseif strcmp(meshSize,'dummy')
        a1 = 0;
    end
    tp1 = a1;
    tq1 = a1;

    % Patch 2 :
    % _________

    if strcmp(meshSize,'coarse') || strcmp(meshSize,'coarsest')
        a2 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a2 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a2 = 0;
    end
    tp2 = a2;
    tq2 = a2;

    % Patch 3 :
    % _________

    if strcmp(meshSize,'coarse') || strcmp(meshSize,'coarsest')
        a3 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a3 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a3 = 0;
    end
    tp3 = a3;
    tq3 = a3;

    % Patch 4 :
    % _________

    if strcmp(meshSize,'coarse') || strcmp(meshSize,'coarsest')
        a4 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a4 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a4 = 0;
    end
    tp4 = a4;
    tq4 = a4;

    % Patch 5 :
    % _________

    if strcmp(meshSize,'coarse') || strcmp(meshSize,'coarsest')
        a5 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a5 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a5 = 0;
    end
    tp5 = a5;
    tq5 = a5;

    % Patch 6 :
    % _________

    if strcmp(meshSize,'coarse') || strcmp(meshSize,'coarsest')
        a6 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a6 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a6 = 0;
    end
    tp6 = a6;
    tq6 = a6;

    % Patch 7 :
    % _________

    if strcmp(meshSize,'coarse') || strcmp(meshSize,'coarsest')
        a7 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a7 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a7 = 0;
    end
    tp7 = a7;
    tq7 = a7;

    % Patch 8 :
    % _________

    if strcmp(meshSize,'coarse') || strcmp(meshSize,'coarsest')
        a8 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a8 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a8 = 0;
    end
    tp8 = a8;
    tq8 = a8;
    
    % Loop over the patches
    for iPatches = 1:noPatches
        tp = eval(['tp' num2str(iPatches)]);
        tq = eval(['tq' num2str(iPatches)]);
        p = eval(['p' num2str(iPatches)]);
        q = eval(['q' num2str(iPatches)]);
        Xi = eval(['Xi' num2str(iPatches)]);
        Eta = eval(['Eta' num2str(iPatches)]);
        CP = eval(['CP' num2str(iPatches)]);
        [Xi,Eta,CP,p,q] = degreeElevateBSplineSurface...
            (p,q,Xi,Eta,CP,tp,tq,'');
        assignin('base',['Xi' num2str(iPatches)],Xi);
        assignin('base',['Eta' num2str(iPatches)],Eta);
        assignin('base',['CP' num2str(iPatches)],CP);
        assignin('base',['p' num2str(iPatches)],p);
        assignin('base',['q' num2str(iPatches)],q);
    end    
end

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%
for iHRef = 1
    % Patch 1 :
    % _________

    if strcmp(meshSize,'coarse')
        a1 = 8; % 13
    elseif strcmp(meshSize,'coarsest')
        a1 = 2; % 2
    elseif strcmp(meshSize,'fine')
        a1 = 27; % 27
    elseif strcmp(meshSize,'dummy')
        a1 = 0; 
    end

    % Patch 2 :
    % _________

    if strcmp(meshSize,'coarse')
        a2 = 8; % 8
    elseif strcmp(meshSize,'coarsest')
        a2 = 2; % 2
    elseif strcmp(meshSize,'fine')
        a2 = 16; % 16
    elseif strcmp(meshSize,'dummy')
        a2 = 0; 
    end

    % Patch 3 :
    % _________

    if strcmp(meshSize,'coarse')
        a3 = 8; % 12
	elseif strcmp(meshSize,'coarsest')
        a3 = 2; % 2
    elseif strcmp(meshSize,'fine')
        a3 = 25; % 25
    elseif strcmp(meshSize,'dummy')
        a3 = 0; 
    end

    % Patch 4 :
    % _________

    if strcmp(meshSize,'coarse')
        a4 = 8; % 9
	elseif strcmp(meshSize,'coarsest')
        a4 = 2; % 2
    elseif strcmp(meshSize,'fine')
        a4 = 18; % 18
    elseif strcmp(meshSize,'dummy')
        a4 = 0; 
    end

    % Patch 5 :
    % _________

    if strcmp(meshSize,'coarse')
        a5 = 8; % 8
	elseif strcmp(meshSize,'coarsest')
        a5 = 2; % 2
    elseif strcmp(meshSize,'fine')
        a5 = 16; % 16
    elseif strcmp(meshSize,'dummy')
        a5 = 0; 
    end

    % Patch 6 :
    % _________

    if strcmp(meshSize,'coarse')
        a6 = 8; % 12
	elseif strcmp(meshSize,'coarsest')
        a6 = 2; % 2
    elseif strcmp(meshSize,'fine')
        a6 = 25; % 25
    elseif strcmp(meshSize,'dummy')
        a6 = 0; 
    end

    % Patch 7 :
    % _________

    if strcmp(meshSize,'coarse')
        a7 = 8; % 9
    elseif strcmp(meshSize,'coarsest')
        a7 = 2; % 2
    elseif strcmp(meshSize,'fine')
        a7 = 18; % 18
    elseif strcmp(meshSize,'dummy')
        a7 = 0; 
    end

    % Patch 8 :
    % _________

    if strcmp(meshSize,'coarse')
        a8 = 8; % 13
    elseif strcmp(meshSize,'coarsest')
        a8 = 2; % 2
    elseif strcmp(meshSize,'fine')
        a8 = 27; % 27
    elseif strcmp(meshSize,'dummy')
        a8 = 0; 
    end
    
    % Loop over the patches
    for iPatches = 1:noPatches
        a = eval(['a' num2str(iPatches)]);
        arcLengthXi = eval(['arcLengthXi' num2str(iPatches)]);
        arcLengthEta = eval(['arcLengthEta' num2str(iPatches)]);
        noKnotsXi = a;
        noKnotsEta = ceil(a*(arcLengthEta/arcLengthXi));
        p = eval(['p' num2str(iPatches)]);
        Xi = eval(['Xi' num2str(iPatches)]);
        q = eval(['q' num2str(iPatches)]);
        Eta = eval(['Eta' num2str(iPatches)]);
        CP = eval(['CP' num2str(iPatches)]);
        [Xi,Eta,CP] = knotRefineUniformlyBSplineSurface...
            (p,Xi,q,Eta,CP,noKnotsXi,noKnotsEta,'');
        assignin('base',['Xi' num2str(iPatches)],Xi);
        assignin('base',['Eta' num2str(iPatches)],Eta);
        assignin('base',['CP' num2str(iPatches)],CP);
    end
end

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Homogeneous Dirichlet boundary conditions
homDOFs = [];

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC.noCnd = 1; % 2
weakDBC.xiExtension = {[0 0]}; % [0 0] [1 1]
weakDBC.etaExtension = {[0 1]}; % [0 1] [0 1]
weakDBC.int.type = 'default';
weakDBC.int.noGPs = 16;

% Loop over the patches
for iPatches = 1:noPatches
    assignin('base',['homDOFs' num2str(iPatches)],homDOFs);
    assignin('base',['inhomDOFs' num2str(iPatches)],inhomDOFs);
    assignin('base',['valuesInhomDOFs' num2str(iPatches)],valuesInhomDOFs);
    assignin('base',['weakDBC' num2str(iPatches)],weakDBC);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cables                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cable
cables.No = 0;

% Loop over all patches
for iPatches = 1:noPatches
    assignin('base',['cables' num2str(iPatches)],cables);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iNBC = 1
    % General parameter
    scaling = 1.0;
    loadAmplitude = - scaling*FAmp;
    scaling = 0; % 2e1, 4e1
    wind = scaling*1e1;

    % Patch 1 :
    % _________

    FAmp1 = loadAmplitude;
    NBC1.noCnd = 1;
    xib1 = [0 1];   etab1 = [0 1];   dirForce1 = 'normal';
    NBC1.xiLoadExtension = {xib1};
    NBC1.etaLoadExtension = {etab1};
    NBC1.loadAmplitude = {FAmp1};
    NBC1.loadDirection = {dirForce1};
    NBC1.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC1.isFollower = true;
    NBC1.isTimeDependent = false;

    % Patch 2 :
    % _________

    FAmp2 = - loadAmplitude;
    NBC2.noCnd = 1;
    xib2 = [0 1];   etab2 = [0 1];   dirForce2 = 'normal';
    NBC2.xiLoadExtension = {xib2};
    NBC2.etaLoadExtension = {etab2};
    NBC2.loadAmplitude = {FAmp2};
    NBC2.loadDirection = {dirForce2};
    NBC2.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC2.isFollower = true;
    NBC2.isTimeDependent = false;

    % Patch 3 :
    % _________

    FAmp3 = loadAmplitude;
    NBC3.noCnd = 1;
    xib3 = [0 1];   etab3 = [0 1];   dirForce3 = 'normal';
    NBC3.xiLoadExtension = {xib3};
    NBC3.etaLoadExtension = {etab3};
    NBC3.loadAmplitude = {FAmp3};
    NBC3.loadDirection = {dirForce3};
    NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC3.isFollower = true;
    NBC3.isTimeDependent = false;

    % Patch 4 :
    % _________

    FAmp4 = - loadAmplitude;
    NBC4.noCnd = 1;
    xib4 = [0 1];   etab4 = [0 1];   dirForce4 = 'normal';
    NBC4.xiLoadExtension = {xib4};
    NBC4.etaLoadExtension = {etab4};
    NBC4.loadAmplitude = {FAmp4};
    NBC4.loadDirection = {dirForce4};
    NBC4.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC4.isFollower = true;
    NBC4.isTimeDependent = false;

    % Patch 5 :
    % _________

    FAmp5 = - loadAmplitude;
    NBC5.noCnd = 1;
    xib5 = [0 1];   etab5 = [0 1];   dirForce5 = 'normal';
    NBC5.xiLoadExtension = {xib5};
    NBC5.etaLoadExtension = {etab5};
    NBC5.loadAmplitude = {FAmp5};
    NBC5.loadDirection = {dirForce5};
    NBC5.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC5.isFollower = true;
    NBC5.isTimeDependent = false;

    % Patch 6 :
    % _________

    FAmp6 = loadAmplitude;
    NBC6.noCnd = 1;
    xib6 = [0 1];   etab6 = [0 1];   dirForce6 = 'normal';
    NBC6.xiLoadExtension = {xib6};
    NBC6.etaLoadExtension = {etab6};
    NBC6.loadAmplitude = {FAmp6 wind};
    NBC6.loadDirection = {dirForce6};
    NBC6.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC6.isFollower = true;
    NBC6.isTimeDependent = false;

    % Patch 7 :
    % _________

    FAmp7 = - loadAmplitude;
    NBC7.noCnd = 1;
    xib7 = [0 1];   etab7 = [0 1];   dirForce7 = 'normal';
    NBC7.xiLoadExtension = {xib7};
    NBC7.etaLoadExtension = {etab7};
    NBC7.loadAmplitude = {FAmp7};
    NBC7.loadDirection = {dirForce7};
    NBC7.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC7.isFollower = true;
    NBC7.isTimeDependent = false;

    % Patch 8 :
    % _________

    FAmp8 = loadAmplitude;
    NBC8.noCnd = 1;
    xib8 = [0 1];   etab8 = [0 1];   dirForce8 = 'normal';
    NBC8.xiLoadExtension = {xib8};
    NBC8.etaLoadExtension = {etab8};
    NBC8.loadAmplitude = {FAmp8};
    NBC8.loadDirection = {dirForce8};
    NBC8.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC8.isFollower = true;
    NBC8.isTimeDependent = false;

    % Collect all the Neumann boundary conditions into an arra<y
    NBC = {NBC1 NBC2 NBC3 NBC4 NBC5 NBC6 NBC7 NBC8};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface parametrizations     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iGamma = 1
    % Patch 1 :
    % _________

    % Connection with patch 2:
    xicoup1_2 = [0 1];
    etacoup1_2 = [0 0];

    % Connection with patch 4:
    xicoup1_4 = [0 1];
    etacoup1_4 = [1 1];

    % Connection with patch 5:
    xicoup1_5 = [1 1];
    etacoup1_5 = [0 1];

    % Connection with patch 10:
    xicoup1_10 = [0 1];
    etacoup1_10 = [1 1];

    % Connection with patch 11:
    xicoup1_11 = [0 1];
    etacoup1_11 = [1 1];

    % Collect all interfaces into arrays:
    xicoup1 = [xicoup1_2 xicoup1_4 xicoup1_5 xicoup1_10 xicoup1_11];
    etacoup1 = [etacoup1_2 etacoup1_4 etacoup1_5 etacoup1_10 etacoup1_11];

    % Patch 2 :
    % _________

    % Connection with patch 1:
    xicoup2_1 = [0 1];
    etacoup2_1 = [0 0];

    % Connection with patch 3:
    xicoup2_3 = [0 1];
    etacoup2_3 = [1 1];

    % Connection with patch 6:
    xicoup2_6 = [1 1];
    etacoup2_6 = [0 1];

    % Collect all interfaces into arrays:
    xicoup2 = [xicoup2_1 xicoup2_3 xicoup2_6];
    etacoup2 = [etacoup2_1 etacoup2_3 etacoup2_6];

    % Patch 3 :
    % _________

    % Connection with patch 2:
    xicoup3_2 = [0 1];
    etacoup3_2 = [1 1];

    % Connection with patch 4:
    xicoup3_4 = [0 1];
    etacoup3_4 = [0 0];

    % Connection with patch 7:
    xicoup3_7 = [1 1];
    etacoup3_7 = [0 1];

    % Collect all interfaces into arrays:
    xicoup3 = [xicoup3_2 xicoup3_4 xicoup3_7];
    etacoup3 = [etacoup3_2 etacoup3_4 etacoup3_7];

    % Patch 4 :
    % _________

    % Connection with patch 1:
    xicoup4_1 = [0 1];
    etacoup4_1 = [1 1];

    % Connection with patch 3:
    xicoup4_3 = [0 1];
    etacoup4_3 = [0 0];

    % Connection with patch 8:
    xicoup4_8 = [1 1];
    etacoup4_8 = [0 1];

    % Connection with patch 10:
    xicoup4_10 = [0 1];
    etacoup4_10 = [1 1];

    % Connection with patch 11:
    xicoup4_11 = [0 1];
    etacoup4_11 = [1 1];

    % Collect all interfaces into arrays:
    xicoup4 = [xicoup4_1 xicoup4_3 xicoup4_8 xicoup4_10 xicoup4_11];
    etacoup4 = [etacoup4_1 etacoup4_3 etacoup4_8 etacoup4_10 etacoup4_11];

    % Patch 5 :
    % _________

    % Connection with patch 1:
    xicoup5_1 = [1 1];
    etacoup5_1 = [0 1];

    % Connection with patch 6:
    xicoup5_6 = [0 1];
    etacoup5_6 = [0 0];

    % Connection with patch 8:
    xicoup5_8 = [0 1];
    etacoup5_8 = [1 1];

    % Connection with patch 14:
    xicoup5_14 = [0 1];
    etacoup5_14 = [1 1];

    % Connection with patch 15:
    xicoup5_15 = [0 1];
    etacoup5_15 = [1 1];

    % Collect all interfaces into arrays:
    xicoup5 = [xicoup5_1 xicoup5_6 xicoup5_8 xicoup5_14 xicoup5_15];
    etacoup5 = [etacoup5_1 etacoup5_6 etacoup5_8 etacoup5_14 etacoup5_15];

    % Patch 6 :
    % _________

    % Connection with patch 2:
    xicoup6_2 = [1 1];
    etacoup6_2 = [0 1];

    % Connection with patch 5:
    xicoup6_5 = [0 1];
    etacoup6_5 = [0 0];

    % Connection with patch 7:
    xicoup6_7 = [0 1];
    etacoup6_7 = [1 1];

    % Collect all interfaces into arrays:
    xicoup6 = [xicoup2_6 xicoup6_5 xicoup6_7];
    etacoup6 = [etacoup2_6 etacoup6_5 etacoup6_7];

    % Patch 7 :
    % _________

    % Connection with patch 3:
    xicoup7_3 = [1 1];
    etacoup7_3 = [0 1];

    % Connection with patch 6:
    xicoup7_6 = [0 1];
    etacoup7_6 = [1 1];

    % Connection with patch 8:
    xicoup7_8 = [0 1];
    etacoup7_8 = [0 0];

    % Collect all interfaces into arrays:
    xicoup7 = [xicoup7_3 xicoup7_6 xicoup7_8];
    etacoup7 = [etacoup7_3 etacoup7_6 etacoup7_8];

    % Patch 8 :
    % _________

    % Connection with patch 4:
    xicoup8_4 = [1 1];
    etacoup8_4 = [0 1];

    % Connection with patch 5:
    xicoup8_5 = [0 1];
    etacoup8_5 = [1 1];

    % Connection with patch 7:
    xicoup8_7 = [0 1];
    etacoup8_7 = [0 0];

    % Connection with patch 14:
    xicoup8_14 = [0 1];
    etacoup8_14 = [1 1];

    % Connection with patch 15:
    xicoup8_15 = [0 1];
    etacoup8_15 = [1 1];

    % Collect all interfaces into arrays:
    xicoup8 = [xicoup8_4 xicoup8_5 xicoup8_7 xicoup8_14 xicoup8_15];
    etacoup8 = [etacoup8_4 etacoup8_5 etacoup8_7 etacoup8_14 etacoup8_15];

    % Patch 9 :
    % _________

    % Connection with patch 10:
    xicoup9_10 = [0 1];
    etacoup9_10 = [0 0];

    % Connection with patch 12:
    xicoup9_12 = [0 1];
    etacoup9_12 = [1 1];

    % Connection with patch 13:
    xicoup9_13 = [1 1];
    etacoup9_13 = [0 1];

    % Connection with patch 18:
    xicoup9_18 = [0 1];
    etacoup9_18 = [1 1];

    % Connection with patch 19:
    xicoup9_19 = [0 1];
    etacoup9_19 = [1 1];

    % Collect all interfaces into arrays:
    xicoup9 = [xicoup9_10 xicoup9_12 xicoup9_13 xicoup9_19];
    etacoup9 = [etacoup9_10 etacoup9_12 etacoup9_13 etacoup9_19];

    % Patch 10 :
    % _________

    % Connection with patch 9:
    xicoup10_9 = [0 1];
    etacoup10_9 = [0 0];

    % Connection with patch 11:
    xicoup10_11 = [0 1];
    etacoup10_11 = [1 1];

    % Connection with patch 14:
    xicoup10_14 = [1 1];
    etacoup10_14 = [0 1];

    % Connection with patch 1:
    xicoup10_1 = [0 1];
    etacoup10_1 = [1 1];

    % Connection with patch 4:
    xicoup10_4 = [0 1];
    etacoup10_4 = [1 1];

    % Collect all interfaces into arrays:
    xicoup10 = [xicoup10_9 xicoup10_11 xicoup10_14 xicoup10_1 xicoup10_4];
    etacoup10 = [etacoup10_9 etacoup10_11 etacoup10_14 etacoup10_1 etacoup10_4];

    % Patch 11 :
    % _________

    % Connection with patch 10:
    xicoup11_10 = [0 1];
    etacoup11_10 = [1 1];

    % Connection with patch 12:
    xicoup11_12 = [0 1];
    etacoup11_12 = [0 0];

    % Connection with patch 15:
    xicoup11_15 = [1 1];
    etacoup11_15 = [0 1];

    % Connection with patch 1:
    xicoup11_1 = [0 1];
    etacoup11_1 = [1 1];

    % Connection with patch 4:
    xicoup11_4 = [0 1];
    etacoup11_4 = [1 1];

    % Collect all interfaces into arrays:
    xicoup11 = [xicoup11_10 xicoup11_12 xicoup11_15 xicoup11_1 xicoup11_4];
    etacoup11 = [etacoup11_10 etacoup11_12 etacoup11_15 etacoup11_1 etacoup11_4];

    % Patch 12 :
    % _________

    % Connection with patch 9:
    xicoup12_9 = [0 1];
    etacoup12_9 = [1 1];

    % Connection with patch 11:
    xicoup12_11 = [0 1];
    etacoup12_11 = [0 0];

    % Connection with patch 16:
    xicoup12_16 = [1 1];
    etacoup12_16 = [0 1];

    % Connection with patch 18:
    xicoup12_18 = [0 1];
    etacoup12_18 = [1 1];

    % Connection with patch 18:
    xicoup12_19 = [0 1];
    etacoup12_19 = [1 1];

    % Collect all interfaces into arrays:
    xicoup12 = [xicoup12_9 xicoup12_11 xicoup12_16 xicoup12_18 xicoup12_19];
    etacoup12 = [etacoup12_9 etacoup12_11 etacoup12_16 etacoup12_18 etacoup12_19];

    % Patch 13 :
    % _________

    % Connection with patch 9:
    xicoup13_9 = [1 1];
    etacoup13_9 = [0 1];

    % Connection with patch 14:
    xicoup13_14 = [0 1];
    etacoup13_14 = [0 0];

    % Connection with patch 16:
    xicoup13_16 = [0 1];
    etacoup13_16 = [1 1];

    % Connection with path 22:
    xicoup13_22 = [0 1];
    etacoup13_22 = [1 1];

    % Connection with path 23:
    xicoup13_23 = [0 1];
    etacoup13_23 = [1 1];

    % Collect all interfaces into arrays:
    xicoup13 = [xicoup13_9 xicoup13_14 xicoup13_16 xicoup13_22 xicoup13_23];
    etacoup13 = [etacoup13_9 etacoup13_14 etacoup13_16 etacoup13_22 etacoup13_23];

    % Patch 14 :
    % _________

    % Connection with patch 10:
    xicoup14_10 = [1 1];
    etacoup14_10 = [0 1];

    % Connection with patch 13:
    xicoup14_13 = [0 1];
    etacoup14_13 = [0 0];

    % Connection with patch 15:
    xicoup14_15 = [0 1];
    etacoup14_15 = [1 1];

    % Connection with patch 5:
    xicoup14_5 = [0 1];
    etacoup14_5 = [1 1];

    % Connection with patch 8:
    xicoup14_8 = [0 1];
    etacoup14_8 = [1 1];

    % Collect all interfaces into arrays:
    xicoup14 = [xicoup14_10 xicoup14_13 xicoup14_15 xicoup14_5 xicoup14_8];
    etacoup14 = [etacoup14_10 etacoup14_13 etacoup14_15 etacoup14_5 etacoup14_8];

    % Patch 15 :
    % _________

    % Connection with patch 11:
    xicoup15_11 = [1 1];
    etacoup15_11 = [0 1];

    % Connection with patch 14:
    xicoup15_14 = [0 1];
    etacoup15_14 = [1 1];

    % Connection with patch 16:
    xicoup15_16 = [0 1];
    etacoup15_16 = [0 0];

    % Connection with patch 5:
    xicoup15_5 = [0 1];
    etacoup15_5 = [1 1];

    % Connection with patch 8:
    xicoup15_8 = [0 1];
    etacoup15_8 = [1 1];

    % Collect all interfaces into arrays:
    xicoup15 = [xicoup15_11 xicoup15_14 xicoup15_16 xicoup15_5 xicoup15_8];
    etacoup15 = [etacoup15_11 etacoup15_14 etacoup15_16 etacoup15_5 etacoup15_8];

    % Patch 16 :
    % _________

    % Connection with patch 12:
    xicoup16_12 = [1 1];
    etacoup16_12 = [0 1];

    % Connection with patch 13:
    xicoup16_13 = [0 1];
    etacoup16_13 = [1 1];

    % Connection with patch 15:
    xicoup16_15 = [0 1];
    etacoup16_15 = [0 0];

    % Connection with patch 22:
    xicoup16_22 = [0 1];
    etacoup16_22 = [1 1];

    % Connection with patch 23:
    xicoup16_23 = [0 1];
    etacoup16_23 = [1 1];

    % Collect all interfaces into arrays:
    xicoup16 = [xicoup16_12 xicoup16_13 xicoup16_15 xicoup16_22 xicoup16_23];
    etacoup16 = [etacoup16_12 etacoup16_13 etacoup16_15 etacoup16_22 etacoup16_23];

    % Patch 17 :
    % _________

    % Connection with patch 18:
    xicoup17_18 = [0 1];
    etacoup17_18 = [0 0];

    % Connection with patch 20:
    xicoup17_20 = [0 1];
    etacoup17_20 = [1 1];

    % Connection with patch 21:
    xicoup17_21 = [1 1];
    etacoup17_21 = [0 1];

    % Collect all interfaces into arrays:
    xicoup17 = [xicoup17_18 xicoup17_20 xicoup17_21];
    etacoup17 = [etacoup17_18 etacoup17_20 etacoup17_21];

    % Patch 18 :
    % _________

    % Connection with patch 17:
    xicoup18_17 = [0 1];
    etacoup18_17 = [0 0];

    % Connection with patch 19:
    xicoup18_19 = [0 1];
    etacoup18_19 = [1 1];

    % Connection with patch 22:
    xicoup18_22 = [1 1];
    etacoup18_22 = [0 1];

    % Connection with patch 9:
    xicoup18_9 = [0 1];
    etacoup18_9 = [1 1];

    % Connection with patch 12:
    xicoup18_12 = [0 1];
    etacoup18_12 = [1 1];

    % Collect all interfaces into arrays:
    xicoup18 = [xicoup18_17 xicoup18_19 xicoup18_22 xicoup18_9 xicoup18_12];
    etacoup18 = [etacoup18_17 etacoup18_19 etacoup18_22 etacoup18_9 etacoup18_12];

    % Patch 19 :
    % _________

    % Connection with patch 18:
    xicoup19_18 = [0 1];
    etacoup19_18 = [1 1];

    % Connection with patch 20:
    xicoup19_20 = [0 1];
    etacoup19_20 = [0 0];

    % Connection with patch 23:
    xicoup19_23 = [1 1];
    etacoup19_23 = [0 1];

    % Connection with patch 9:
    xicoup19_9 = [0 1];
    etacoup19_9 = [1 1];

    % Connection with patch 12:
    xicoup19_12 = [0 1];
    etacoup19_12 = [1 1];

    % Collect all interfaces into arrays:
    xicoup19 = [xicoup19_18 xicoup19_20 xicoup19_23 xicoup19_9 xicoup19_12];
    etacoup19 = [etacoup19_18 etacoup19_20 etacoup19_23 etacoup19_9 etacoup19_12];

    % Patch 20 :
    % _________

    % Connection with patch 17:
    xicoup20_17 = [0 1];
    etacoup20_17 = [1 1];

    % Connection with patch 19:
    xicoup20_19 = [0 1];
    etacoup20_19 = [0 0];

    % Connection with patch 24:
    xicoup20_24 = [1 1];
    etacoup20_24 = [0 1];

    % Collect all interfaces into arrays:
    xicoup20 = [xicoup20_17 xicoup20_19 xicoup20_24];
    etacoup20 = [etacoup20_17 etacoup20_19 etacoup20_24];

    % Patch 21 :
    % _________

    % Connection with patch 17:
    xicoup21_17 = [1 1];
    etacoup21_17 = [0 1];

    % Connection with patch 22:
    xicoup21_22 = [0 1];
    etacoup21_22 = [0 0];

    % Connection with patch 24:
    xicoup21_24 = [0 1];
    etacoup21_24 = [1 1];

    % Collect all interfaces into arrays:
    xicoup21 = [xicoup21_17 xicoup21_22 xicoup21_24];
    etacoup21 = [etacoup21_17 etacoup21_22 etacoup21_24];

    % Patch 22 :
    % _________

    % Connection with patch 18:
    xicoup22_18 = [1 1];
    etacoup22_18 = [0 1];

    % Connection with patch 21:
    xicoup22_21 = [0 1];
    etacoup22_21 = [0 0];

    % Connection with patch 23:
    xicoup22_23 = [0 1];
    etacoup22_23 = [1 1];

    % Connection with patch 13:
    xicoup22_13 = [0 1];
    etacoup22_13 = [1 1];

    % Connection with patch 16:
    xicoup22_16 = [0 1];
    etacoup22_16 = [1 1];

    % Collect all interfaces into arrays:
    xicoup22 = [xicoup22_18 xicoup22_21 xicoup22_23 xicoup22_13 xicoup22_16];
    etacoup22 = [etacoup22_18 etacoup22_21 etacoup22_23 etacoup22_13 etacoup22_16];

    % Patch 23 :
    % _________

    % Connection with patch 19:
    xicoup23_19 = [1 1];
    etacoup23_19 = [0 1];

    % Connection with patch 22:
    xicoup23_22 = [0 1];
    etacoup23_22 = [1 1];

    % Connection with patch 24:
    xicoup23_24 = [0 1];
    etacoup23_24 = [0 0];

    % Connection with patch 13:
    xicoup23_13 = [0 1];
    etacoup23_13 = [1 1];

    % Connection with patch 16:
    xicoup23_16 = [0 1];
    etacoup23_16 = [1 1];

    % Collect all interfaces into arrays:
    xicoup23 = [xicoup23_19 xicoup23_22 xicoup23_24 xicoup23_13 xicoup23_16];
    etacoup23 = [etacoup23_19 etacoup23_22 etacoup23_24 etacoup23_13 etacoup23_16];

    % Patch 24 :
    % _________

    % Connection with patch 20:
    xicoup24_20 = [1 1];
    etacoup24_20 = [0 1];

    % Connection with patch 21:
    xicoup24_21 = [0 1];
    etacoup24_21 = [1 1];

    % Connection with patch 23:
    xicoup24_23 = [0 1];
    etacoup24_23 = [0 0];

    % Collect all interfaces into arrays:
    xicoup24 = [xicoup24_20 xicoup24_21 xicoup24_23];
    etacoup24 = [etacoup24_20 etacoup24_21 etacoup24_23];

    % Define connections :
    % ____________________

    % Number of connections
    noConnections = 12; % 12 % 24 % 44 % 52

    % Define connections by patch numbers
    connections.No = noConnections;
    connections.xiEtaCoup = zeros(noConnections,10);
    connections.xiEtaCoup(:,:) = [1 2 xicoup1_2 etacoup1_2 xicoup2_1 etacoup2_1 % -> First tube
                                  2 3 xicoup2_3 etacoup2_3 xicoup3_2 etacoup3_2
                                  3 4 xicoup3_4 etacoup3_4 xicoup4_3 etacoup4_3
                                  1 4 xicoup1_4 etacoup1_4 xicoup4_1 etacoup4_1
                                  1 5 xicoup1_5 etacoup1_5 xicoup5_1 etacoup5_1
                                  2 6 xicoup2_6 etacoup2_6 xicoup6_2 etacoup6_2
                                  3 7 xicoup3_7 etacoup3_7 xicoup7_3 etacoup7_3
                                  4 8 xicoup4_8 etacoup4_8 xicoup8_4 etacoup8_4
                                  5 6 xicoup5_6 etacoup5_6 xicoup6_5 etacoup6_5
                                  6 7 xicoup6_7 etacoup6_7 xicoup7_6 etacoup7_6
                                  7 8 xicoup7_8 etacoup7_8 xicoup8_7 etacoup8_7
                                  5 8 xicoup5_8 etacoup5_8 xicoup8_5 etacoup8_5]; % <- First tube
%                                   9 10 xicoup9_10 etacoup9_10 xicoup10_9 etacoup10_9 % <- Second tube
%                                   10 11 xicoup10_11 etacoup10_11 xicoup11_10 etacoup11_10
%                                   11 12 xicoup11_12 etacoup11_12 xicoup12_11 etacoup12_11
%                                   9 12 xicoup9_12 etacoup9_12 xicoup12_9 etacoup12_9
%                                   9 13 xicoup9_13 etacoup9_13 xicoup13_9 etacoup13_9
%                                   10 14 xicoup10_14 etacoup10_14 xicoup14_10 etacoup14_10
%                                   11 15 xicoup11_15 etacoup11_15 xicoup15_11 etacoup15_11
%                                   12 16 xicoup12_16 etacoup12_16 xicoup16_12 etacoup16_12
%                                   13 14 xicoup13_14 etacoup13_14 xicoup14_13 etacoup14_13
%                                   14 15 xicoup14_15 etacoup14_15 xicoup15_14 etacoup15_14
%                                   15 16 xicoup15_16 etacoup15_16 xicoup16_15 etacoup16_15
%                                   13 16 xicoup13_16 etacoup13_16 xicoup16_13 etacoup16_13 % <- Second tube
%                                   17 18 xicoup17_18 etacoup17_18 xicoup18_17 etacoup18_17 % -> Third Tube
%                                   18 19 xicoup18_19 etacoup18_19 xicoup19_18 etacoup19_18
%                                   19 20 xicoup19_20 etacoup19_20 xicoup20_19 etacoup20_19
%                                   17 20 xicoup17_20 etacoup17_20 xicoup20_17 etacoup20_17
%                                   17 21 xicoup17_21 etacoup17_21 xicoup21_17 etacoup21_17
%                                   18 22 xicoup18_22 etacoup18_22 xicoup22_18 etacoup22_18
%                                   19 23 xicoup19_23 etacoup19_23 xicoup23_19 etacoup23_19
%                                   20 24 xicoup20_24 etacoup20_24 xicoup24_20 etacoup24_20
%                                   21 22 xicoup21_22 etacoup21_22 xicoup22_21 etacoup22_21
%                                   22 23 xicoup22_23 etacoup22_23 xicoup23_22 etacoup23_22
%                                   23 24 xicoup23_24 etacoup23_24 xicoup24_23 etacoup24_23
%                                   21 24 xicoup21_24 etacoup21_24 xicoup24_21 etacoup24_21 % <- Third Tube
%                                   5 14 xicoup5_14 etacoup5_14 xicoup14_5 etacoup14_5 % -> Coupling betwenn the first and the second tube
%                                   5 15 xicoup5_15 etacoup5_15 xicoup15_5 etacoup15_5
%                                   8 14 xicoup8_14 etacoup8_14 xicoup14_8 etacoup14_8
%                                   8 15 xicoup8_15 etacoup8_15 xicoup15_8 etacoup15_8
%                                   1 10 xicoup1_10 etacoup1_10 xicoup10_1 etacoup10_1
%                                   1 11 xicoup1_11 etacoup1_11 xicoup11_1 etacoup11_1
%                                   4 10 xicoup4_10 etacoup4_10 xicoup10_4 etacoup10_4
%                                   4 11 xicoup4_11 etacoup4_11 xicoup11_4 etacoup11_4 % <- Coupling betwenn the first and the second tube
%                                   13 22 xicoup13_22 etacoup13_22 xicoup22_13 etacoup22_13 % -> Coupling betwenn the second and the third tube
%                                   13 23 xicoup13_23 etacoup13_23 xicoup23_13 etacoup23_13
%                                   16 22 xicoup16_22 etacoup16_22 xicoup22_16 etacoup22_16
%                                   16 23 xicoup16_23 etacoup16_23 xicoup23_16 etacoup23_16
%                                   9 18 xicoup9_18 etacoup9_18 xicoup18_9 etacoup18_9
%                                   9 19 xicoup9_19 etacoup9_19 xicoup19_9 etacoup19_9
%                                   12 18 xicoup12_18 etacoup12_18 xicoup18_12 etacoup18_12
%                                   12 19 xicoup12_19 etacoup12_19 xicoup19_12 etacoup19_12]; % <- Coupling betwenn the second and the third tube
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connections                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
connections.masterSlave = true(12,1);
connectionsLM = connections;
connectionsALM = connections;

%% Fill up the arrays for the patches
BSplinePatches = struct([]);
for iPatches = 1:noPatches
    p = eval(['p' num2str(iPatches)]);
    Xi = eval(['Xi' num2str(iPatches)]);
    q = eval(['q' num2str(iPatches)]);
    Eta = eval(['Eta' num2str(iPatches)]);
    CP = eval(['CP' num2str(iPatches)]);
    isNURBS = eval(['isNURBS' num2str(iPatches)]);
    homDOFs = eval(['homDOFs' num2str(iPatches)]);
    inhomDOFs = eval(['inhomDOFs' num2str(iPatches)]);
    valuesInhomDOFs = eval(['valuesInhomDOFs' num2str(iPatches)]);
    weakDBC = eval(['weakDBC' num2str(iPatches)]);
    cables = eval(['cables' num2str(iPatches)]);
    NBCPatch = eval(['NBC' num2str(iPatches)]);
    xicoup = [];
    etacoup = [];
    int = eval(['int' num2str(iPatches)]);
    patch = fillUpPatch...
        (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,...
        inhomDOFs,valuesInhomDOFs,weakDBC,cables,NBCPatch,[],[],[],...
        xicoup,etacoup,int);
    assignin('base',['patch' num2str(iPatches)],patch);
    BSplinePatches{iPatches} = patch;
end

%% Plot the multipatch geometry together with the parametric coordinates
% color = [217 218 219]/255;
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,color,graph);

%% Compute the load vector for the visualization of the reference configuration
for counterPatches = 1:noPatches
    BSplinePatches{counterPatches}.FGamma = ...
        zeros(3*BSplinePatches{counterPatches}.noCPs,1);
%     for iNBC = 1:NBC{counterPatches}.noCnd
%         funcHandle = str2func(NBC{counterPatches}.computeLoadVct{iNBC});
%         BSplinePatches{counterPatches}.FGamma = funcHandle...
%             (BSplinePatches{counterPatches}.FGamma,...
%             BSplinePatches{counterPatches},...
%             NBC{counterPatches}.xiLoadExtension{iNBC},...
%             NBC{counterPatches}.etaLoadExtension{iNBC},...
%             NBC{counterPatches}.loadAmplitude{iNBC},...
%             NBC{counterPatches}.loadDirection{iNBC},...
%             NBC{counterPatches}.isFollower(iNBC,1),0,...
%             BSplinePatches{counterPatches}.int,'outputEnabled');
%     end
end

%% Plot the reference configuration for the multipatch geometry
% color = [217 218 219]/255;
% graph.isPrestressEnabled = true;
% graph.compPrestress = '2';
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');
% az = -37.500000000000000;
% el = 30;
% view(az,el);

%% Plot the torus
% figure(graph.index)
% 
% % Plot first tube
% color = [217 218 219]/255;
% center = [0; 0; 0];
% Radius = Radius;
% radius = radius;
% thetaInterval = [-pi/2 pi/2];
% phiInterval = [0 +2*pi]; % [0 2*pi]
% azimuthallySymmetryAxis = 'y';
% compPrestress = 2;
% plot_torus3d(center,Radius,radius,thetaInterval,phiInterval,...
%     azimuthallySymmetryAxis,parameters.prestress,compPrestress,color);
% axis equal;
% hold on;
% % hold off;

%% Create Lagrange Multipiers fields for all interfaces
if strcmp(method,'LagrangeMultipliers') || strcmp(method,'AugmentedLagrangeMultipliers')
    fprintf('**************************************************\n');
    fprintf('* Creating interface Lagrange Multipliers fields *\n');
    fprintf('**************************************************\n\n');
    for iConnections = 1:connections.No
        %% Get the IDs of the patches involved
        idI = connections.xiEtaCoup(iConnections,1);
        idJ = connections.xiEtaCoup(iConnections,2);
        fprintf(['\t' 'Coupling between patches %d and %d \n'],idI,idJ);
        fprintf(['\t' '---------------------------------- \n\n']);

        %% Create a basic Lagrange Multipliers field
        pLambda = 0;
        XiLambda = [0 1];
        CPLambda(:,4) = [1];
        isNURBSLambda = 0;
        nxiLambda = length(CPLambda(:,1,1));
        for i = 1:nxiLambda
            if CPLambda(i,4)~=1
                isNURBSLambda = 1;
                break;
            end
        end

        %% Find the interface parametrization for the involved patches
        xiCoupI = connections.xiEtaCoup(iConnections,3:4);
        etaCoupI = connections.xiEtaCoup(iConnections,5:6);
        if xiCoupI(1,1) ~= xiCoupI(1,2) && etaCoupI(1,1) == etaCoupI(1,2)
            isOnXiI = true;
            fprintf(['\t' 'The coupling interface of patch %d is along xi\n'],idI);
        elseif xiCoupI(1,1) == xiCoupI(1,2) && etaCoupI(1,1) ~= etaCoupI(1,2)
            isOnXiI = false;
            fprintf(['\t' 'The coupling interface of patch %d is along eta\n'],idI);
        else
            error('Either the xi or the eta direction has to be fixed along the interface for patch %d',idI);
        end
        xiCoupJ = connections.xiEtaCoup(iConnections,7:8);
        etaCoupJ = connections.xiEtaCoup(iConnections,9:10);
        if xiCoupJ(1,1) ~= xiCoupJ(1,2) && etaCoupJ(1,1) == etaCoupJ(1,2)
            isOnXiJ = true;
            fprintf(['\t' 'The coupling interface of patch %d is along xi\n'],idJ);
        elseif xiCoupJ(1,1) == xiCoupJ(1,2) && etaCoupJ(1,1) ~= etaCoupJ(1,2)
            isOnXiJ = false;
            fprintf(['\t' 'The coupling interface of patch %d is along eta\n'],idJ);
        else
            error('Either the xi or the eta direction has to be fixed along the interface for patch %d',idJ);
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
        fprintf(['\t' 'Degree elevating the Lagrange Multipliers field to %d\n'],pLM);
        if pLM > 0
            clear pLambda XiLambda CPLambda;
            pLambda = 1;
            XiLambda = [0 0 1 1];
            CPLambda(:,4) = [1 1];
            isNURBSLambda = 0;
            nxiLambda = length(CPLambda(:,1,1));
            for i = 1:nxiLambda
                if CPLambda(i,4) ~= 1
                    isNURBSLambda = 1;
                    break;
                end
            end

            % Perform accordingly a p-refinement
            tpLambda = pLM;
            [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
                (pLambda,XiLambda,CPLambda,tpLambda,'');
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
        scaleLM = 1.1;
        
        %% HACK %%
        %% HACK %%
        %% HACK %%
        noLambda = ceil(max([noKnotsI noKnotsJ])*scaleLM);
        %% HACK %%
        %% HACK %%
        %% HACK %%
        
        fprintf(['\t' 'Uniformly inserting %d knots in the Lagrange Multipliers field\n'],noLambda);
        [XiLambdaLM,CPLambdaLM] = knotRefineUniformlyBSplineCurve...
            (noLambda,pLambda,XiLambda,CPLambda,'');
        scaleALM = 0.8;
        noLambda = ceil(max([noKnotsI noKnotsJ])*scaleALM);
        fprintf(['\t' 'Uniformly inserting %d knots in the augmented Lagrange Multipliers field\n'],noLambda);
        [XiLambdaALM,CPLambdaALM] = knotRefineUniformlyBSplineCurve...
            (noLambda,pLambda,XiLambda,CPLambda,'');

        %% Fill up the Lagrange Multipliers patch and add it to the array
        lambdaLM = fillUpLagrangeMultipliers...
            (pLambda,XiLambdaLM,CPLambdaLM,isNURBSLambda);
        connectionsLM.lambda{iConnections} = lambdaLM;
        lambdaALM = fillUpLagrangeMultipliers...
            (pLambda,XiLambdaALM,CPLambdaALM,isNURBSLambda);
        connectionsALM.lambda{iConnections} = lambdaALM;
        fprintf('\n');
    end
end

%% Set up the parameters and properties for each method
if strcmp(method,'Penalty')
    % General parameters
    penaltyPrmScale = 1e0;

    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    for iPatches = 1:noPatches
        % Check if weak boundary conditions are to be applied
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
        % Assign the method name
        BSplinePatches{iPatches}.weakDBC.method = 'penalty';

        % Get the polynomial order along the Dirichlet boundary
        isOnXi = false;
        if BSplinePatches{iPatches}.weakDBC.etaExtension{1}(1,1) == ...
                BSplinePatches{iPatches}.weakDBC.etaExtension{1}(1,2)
            isOnXi = true;
        end
        if isOnXi
            polOrder = BSplinePatches{iPatches}.p;
        else
            polOrder = BSplinePatches{iPatches}.q;
        end

        % Assign the penalty parameter
        BSplinePatches{iPatches}.weakDBC.alpha = ...
            norm(eval(['Dm' num2str(iPatches)]))*polOrder/...
            sqrt(BSplinePatches{iPatches}.minElArea)*...
            penaltyPrmScale;
    end

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'penalty';
    propCoupling.intC = intC;
    propCoupling.alphaD = zeros(connections.No,1);
    propCoupling.alphaR = zeros(connections.No,1);
    for iConnections = 1:connections.No
        % Get the id's of the patches
        IDPatchI = connections.xiEtaCoup(iConnections,1);
        IDPatchJ = connections.xiEtaCoup(iConnections,2);

        % Get the mean polynomial order between the patches
        isOnXiI = false;
        if connections.xiEtaCoup(iConnections,5) == ...
                connections.xiEtaCoup(iConnections,6)
            isOnXiI = true;
        end
        if isOnXiI
            polOrderI = BSplinePatches{IDPatchI}.p;
        else
            polOrderI = BSplinePatches{IDPatchI}.q;
        end
        isOnXiJ = false;
        if connections.xiEtaCoup(iConnections,9) == ...
                connections.xiEtaCoup(iConnections,10)
            isOnXiJ = true;
        end
        if isOnXiJ
            polOrderJ = BSplinePatches{IDPatchJ}.p;
        else
            polOrderJ = BSplinePatches{IDPatchJ}.q;
        end
        polOrderMean = mean([polOrderI polOrderJ]);

        % Assign the penalty parameters
        propCoupling.alphaD(iConnections,1) = ...
            max([norm(eval(['Dm' num2str(IDPatchI)])) ...
            norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
            sqrt(min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea]))*...
            penaltyPrmScale;
    end
elseif strcmp(method,'LagrangeMultipliers')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    % Properties for the weak Dirichlet boundary conditions
    for iPatches = 1:noPatches
        % Check if weak boundary conditions are to be applied
        
        %% HACK %%
        %% HACK %%
        %% HACK %%
%         BSplinePatches{iPatches}.weakDBC.noCnd = 0;
        %% HACK %%
        %% HACK %%
        %% HACK %%
        
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
        % Assign the method name
        BSplinePatches{iPatches}.weakDBC.method = 'lagrangeMultipliers';
        BSplinePatches{iPatches}.weakDBC.alpha = 0;

        % Find along which parametric line the weak Dirichlet 
        % condition is to be imposed
        isOnXi = false;
        if BSplinePatches{iPatches}.weakDBC.xiExtension{1}(1,1) == ...
             BSplinePatches{iPatches}.weakDBC.xiExtension{1}(1,2)
            isOnXi = true;
        end

        % Make a Lagrange Multipliers discretization
        clear pLambda XiLambda CPLambda isNURBSLambda; 

        % Find out up to which polynomial degree the Lagrange
        % Multipliers discretization needs to be increased
        if isOnXi
            polOrderPatch =  BSplinePatches{iPatches}.p;
        else
            polOrderPatch =  BSplinePatches{iPatches}.q;
        end
        pLM = polOrderPatch - 2;

        if pLM <= 0
            pLambda = 0;
            XiLambda = [0 1];
            CPLambda = zeros(1,4);
        else
            pLambda = 1;
            XiLambda = [0 0 1 1];
            CPLambda = zeros(2,4);

            % Peform a p-refinement
            tpLambda = polOrderPatch - pLM;
            [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
                (pLambda,XiLambda,CPLambda,tpLambda,'');
        end
        isNURBSLambda = 0;
        nxiLambda = length(CPLambda(:,1));
        for i = 1:nxiLambda
            if CPLambda(i,4) ~= 1
                isNURBSLambda = 1;
                break;
            end
        end

        % Perform an h-refinement
        percentage = 1.0;
        if isOnXi
            Rxi = unique(BSplinePatches{iPatches}.Xi);
        else
            Rxi = unique(BSplinePatches{iPatches}.Eta);
        end
        noXi = ceil(percentage*(length(Rxi) - 2));
        [XiLambda,CPLambda] = knotRefineUniformlyBSplineCurve...
            (noXi,pLambda,XiLambda,CPLambda,'');

        % Create the corresponding Lagrange Multipliers structure
        BSplinePatches{iPatches}.weakDBC.lambda{1} = fillUpLagrangeMultipliers...
         (pLambda,XiLambda,CPLambda,isNURBSLambda);
    end

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'lagrangeMultipliers';
    connections = connectionsLM;
    propCoupling.alphaD = zeros(connections.No,1);
    propCoupling.alphaR = zeros(connections.No,1);
    propCoupling.intC = intC;
elseif strcmp(method,'Mortar')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    % Properties for the weak Dirichlet boundary conditions
    for iPatches = 1:noPatches
        % Check if weak boundary conditions are to be applied
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
        % Get the boundary extensions
        for iCnd = 1:BSplinePatches{iPatches}.weakDBC.noCnd
            xisup = [BSplinePatches{iPatches}.weakDBC.xiExtension{iCnd}(1,1) BSplinePatches{iPatches}.weakDBC.xiExtension{iCnd}(1,2)];
            etasup = [BSplinePatches{iPatches}.weakDBC.etaExtension{iCnd}(1,1) BSplinePatches{iPatches}.weakDBC.etaExtension{iCnd}(1,2)];

            for dir = 1:3
                BSplinePatches{iPatches}.homDOFs = findDofs3D...
                    (BSplinePatches{iPatches}.homDOFs,xisup,etasup,dir,BSplinePatches{iPatches}.CP);
            end
        end
        
        % Remove the application of weak Dirichlet boundary conditions
        BSplinePatches{iPatches}.weakDBC.noCnd = 0;
    end

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'mortar';
    propCoupling.isSlaveSideCoarser = false;
    propCoupling.computeRearrangedProblemMtrcs = @computeRearrangedProblemMtrcs4MortarIGAMembrane;
    propCoupling.intC = intC;
elseif strcmp(method,'AugmentedLagrangeMultipliers')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------
    
    % Scaling factor for the augmented Lagrange Multipliers method
    scalePenaltyFactorALM = 1e-2;
    
    for iPatches = 1:noPatches
        % Check if weak boundary conditions are to be applied
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
        % Assign the method name
        BSplinePatches{iPatches}.weakDBC.method = 'lagrangeMultipliers';

        % Find along which parametric line the weak Dirichlet
        % condition is to be imposed
        isOnXi = false;
        if BSplinePatches{iPatches}.weakDBC.xiExtension{1}(1,1) == ...
             BSplinePatches{iPatches}.weakDBC.xiExtension{1}(1,2)
            isOnXi = true;
        end

        % Find out up to which polynomial degree the Lagrange
        % Multipliers discretization needs to be increased
        if isOnXi
            polOrderPatch = BSplinePatches{iPatches}.p;
        else
            polOrderPatch = BSplinePatches{iPatches}.q;
        end
        pLM = polOrderPatch - 2; %polOrderPatch - 2 has worked well

        % Assign the penalty parameter
        BSplinePatches{iPatches}.weakDBC.alpha = ...
            norm(eval(['Dm' num2str(iPatches)]))*polOrderPatch/...
            sqrt(BSplinePatches{iPatches}.minElArea)*scalePenaltyFactorALM;

        % Make a Lagrange Multipliers discretization
        clear pLambda XiLambda CPLambda isNURBSLambda; 
        if pLM <= 0
            pLambda = 0;
            XiLambda = [0 1];
            CPLambda = zeros(1,4);
        else
            pLambda = 1;
            XiLambda = [0 0 1 1];
            CPLambda = zeros(2,4);

            % Peform a p-refinement
            tpLambda = polOrderPatch - pLM;
            [XiLambda,CPLambda,pLambda] = degreeElevateBSplineCurve...
                (pLambda,XiLambda,CPLambda,tpLambda,'');
        end
        isNURBSLambda = 0;
        nxiLambda = length(CPLambda(:,1));
        for i = 1:nxiLambda
            if CPLambda(i,4) ~= 1
                isNURBSLambda = 1;
                break;
            end
        end

        % Perform an h-refinement
        percentage = .5;
        if isOnXi
            Rxi = unique(BSplinePatches{iPatches}.Xi);
        else
            Rxi = unique(BSplinePatches{iPatches}.Eta);
        end
        noXi = ceil(percentage*(length(Rxi) - 2));
        [XiLambda,CPLambda] = knotRefineUniformlyBSplineCurve...
            (noXi,pLambda,XiLambda,CPLambda,'');

        % Create the corresponding Lagrange Multipliers structure
        BSplinePatches{iPatches}.weakDBC.lambda{1} = fillUpLagrangeMultipliers...
            (pLambda,XiLambda,CPLambda,isNURBSLambda);
    end
    
    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'lagrangeMultipliers';
    connections = connectionsALM;
    propCoupling.intC = intC;
    propCoupling.alphaD = zeros(connections.No,1);
    propCoupling.alphaR = zeros(connections.No,1);
    for iConnections = 1:connections.No
        % Get the Patch IDs
        IDPatchI = connections.xiEtaCoup(iConnections,1);
        IDPatchJ = connections.xiEtaCoup(iConnections,2);

        % Get the mean polynomial order between the patches
        isOnXiI = false;
        if connections.xiEtaCoup(iConnections,5) == ...
                connections.xiEtaCoup(iConnections,6)
            isOnXiI = true;
        end
        if isOnXiI
            polOrderI = BSplinePatches{IDPatchI}.p;
        else
            polOrderI = BSplinePatches{IDPatchI}.q;
        end
        isOnXiJ = false;
        if connections.xiEtaCoup(iConnections,9) == ...
                connections.xiEtaCoup(iConnections,10)
            isOnXiJ = true;
        end
        if isOnXiJ
            polOrderJ = BSplinePatches{IDPatchJ}.p;
        else
            polOrderJ = BSplinePatches{IDPatchJ}.q;
        end
        polOrderMean = mean([polOrderI polOrderJ]);

        % Assign the penalty parameter
        propCoupling.alphaD(iConnections,1) = ...
            max([norm(eval(['Dm' num2str(IDPatchI)])) ...
            norm(eval(['Dm' num2str(IDPatchJ)]))])*polOrderMean/...
            sqrt(min([BSplinePatches{IDPatchI}.minElArea BSplinePatches{IDPatchJ}.minElArea]))*...
            scalePenaltyFactorALM;
    end
elseif strcmp(method,'Nitsche')
    % Assign the parameters for the application of weak DBC
    % -----------------------------------------------------

    % Properties for the weak Dirichlet boundary conditions
    for iPatches = 1:noPatches
        % Check if weak boundary conditions are to be applied
        if BSplinePatches{iPatches}.weakDBC.noCnd == 0
            continue;
        end
        
        % Assign the method
        BSplinePatches{iPatches}.weakDBC.method = 'nitsche';
        
        % Assign the estimation of the stabilization parameter
        BSplinePatches{iPatches}.weakDBC.estimationStabilPrm = true;
    end

    % Assign the parameters for multipatch coupling
    % ---------------------------------------------

    propCoupling.method = 'nitsche';
    propCoupling.estimationStabilPrm = true;
    propCoupling.gammaTilde = .5;
    propCoupling.intC = intC;
end

%% Form-finding properties
propFormFinding.tolerance = 1e-4;
propFormFinding.maxNoIter = 6;
propFormFinding.minNoIter = 2;

%% Solve the form-finding problem
[BSplinePatchesFoFi,CPHistoryMultipatch,propCouplingNitsche,resHistory,...
    hasConverged,noIter] = ...
    solve_DDMFormFindingIGAMembrane(BSplinePatches,connections,...
    propCoupling,propFormFinding,solve_LinearSystem,'outputEnabled');
% save(['data_FoFiDDMMiddleSailOlympiadach' '_' meshSize '_' method]);

%% Postprocessing
graph.postprocConfig = 'reference';
graph.resultant = 'displacement';
graph.component = '2norm';
color = [217 218 219]/255;
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatchesFoFi,connections,color,graph,'outputEnabled');
az = -37.500000000000000;
el = 30;
view(az,el);
axis off;
title('');

%% END