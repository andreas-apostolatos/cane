%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : A rectangular 2-patch geometry representing a membrane is clamped
%        at both its edges and is subject to a surface load. For the
%        coupling of the multipatches a Penalty, a Lagrange Multipliers and
%        a Nitsche method are employed. The analytical solution for this
%        problem must be 0.880339146912915 for the tip displacement
%
% Date : 28.10.2015
%
%% Preamble
clear;
clc;

%% Includes 

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Add system solvers
addpath('../../equationSystemSolvers/');

% Add functions related to the efficient computation
addpath('../../efficientComputation/');

% Add transient analysis solvers
addpath('../../transientAnalysis/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the isogeometric Kirchhoff-Love shell formulation
addpath('../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../isogeometricThinStructureAnalysis/graphicsMultipatches/',...
        '../../isogeometricThinStructureAnalysis/loads/',...
        '../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../isogeometricThinStructureAnalysis/solvers/',...
        '../../isogeometricThinStructureAnalysis/metrics/',...
        '../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../isogeometricThinStructureAnalysis/penaltyDecompositionKLShell/',...
        '../../isogeometricThinStructureAnalysis/penaltyDecompositionMembrane/',...
        '../../isogeometricThinStructureAnalysis/lagrangeMultipliersDecompositionKLShell/',...
        '../../isogeometricThinStructureAnalysis/nitscheDecompositionKLShell/',...
        '../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../isogeometricThinStructureAnalysis/precomputedData/',...
        '../../isogeometricThinStructureAnalysis/output/',...
        '../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../isogeometricThinStructureAnalysis/weakDBCMembrane/');

% Add functions related to the visualization of the multipatch geometry
% related to the isogeometric mortar-based mapping
addpath('../../isogeometricMortarBasedMappingAnalysis/graphics/');

%% CAD modelling via NURBS

% Global variables:
Length = 10;
Width = 1;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

% Control Point coordinates and weights

% x-coordinates
CP(:,:,1) = [-Length/2 -Length/2
             Length/2  Length/2];

% y-coordinates
CP(:,:,2) = [-Width/2 Width/2
             -Width/2 Width/2];

% z-coordinates
CP(:,:,3) = [0 0
             0 0];
          
% Weights
CP(:,:,4) = [1 1
             1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i= 1:nxi
    for j = 1:neta
        if CP(i,j,4)~=1
            isNURBS = 1;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% Material constants

% Parameters
EYoung = 1e2;
poissonRatio = 0.0;
thickness = 1;
density = 7810;
prestress.voigtVector  = [0
                          0
                          0];
                      
% Young's modulus
parameters.E = EYoung;

% Poisson ratio
parameters.nue = poissonRatio;

% Thickness of the plate
parameters.t = thickness;

% Density of the membrane (used only for dynamics)
parameters.rho = density;

% Prestress for the membrane
parameters.prestress = prestress;

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Case name
caseName = 'rectangularPlate';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points

int.type = 'default';
if strcmp(int.type,'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.etaNGPForLoad = 6;
    int.nGPForLoad = 6;
    int.nGPError = 12;
end

% On the nonlinear equation system solver
newtonRaphson.eps = 1e-14;
newtonRaphson.maxIt = 10;

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','force'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x','y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

%% Refinement

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

tp = 0;
tq = 0;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

a = 1;
nKnotsXi = a;
nKnotsEta = ceil(a*Width/Length);
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,nKnotsXi,nKnotsEta,'');

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Homogeneous DOFs
homDOFs = [];
xisup = [0 0];   
etasup = [0 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xisup,etasup,dirSupp,CP);
end
xisup = [0 1];
etasup = [0 1];
for dirSupp = [2 3]
    homDOFs = findDofs3D(homDOFs,xisup,etasup,dirSupp,CP);
end

% Inhomogeneous DOFs
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;
% weakDBC.method = 'Nitsche';
% weakDBC.scalingStabilPrm = 0;
% weakDBC.alpha = 0.0;
% weakDBC.xiExtension = {[0 0]};
% weakDBC.etaExtension = {[0 1]};
% weakDBC.computeConstMtx = 'computeWeakDBCMtxPenaltyIGAMembrane';
% % weakDBC1.computeConstMtx = 'computeWeakDBCConstStabilMtxNitscheIGAMembrane';
% weakDBC.computeTangMtxResVct = 'computeWeakDBCTangMtxResVctNitscheIGAMembrane';
% weakDBC.int.type = 'default';
% weakDBC.int.noGPs = 16;

% Cables
cables.No = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load amplitude
FAmp = + 1e1;

% Patch 1 :
% _________

FAmp = FAmp;
NBC.noCnd = 1;
xib = 1;   etab = [0 1];   dirForce = 'x';
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude{1} = FAmp;
NBC.loadDirection = {dirForce};
NBC.computeLoadVct{1} = 'computeLoadVctLineIGAThinStructure';
NBC.isFollower(1,1) = false;
NBC.isTimeDependent(1,1) = false;

%% Fill up the patch array
patch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,...
    inhomDOFs,valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch};
noPatches = length(BSplinePatches);

connections.No = 0;

NBC = {NBC};

%% Compute the load vector for the visualization of the reference configuration
for counterPatches = 1:noPatches
    BSplinePatches{counterPatches}.FGamma = ...
        zeros(3*BSplinePatches{counterPatches}.noCPs,1);
    for counterNBC = 1:NBC{counterPatches}.noCnd
        funcHandle = str2func(NBC{counterPatches}.computeLoadVct{counterNBC});
        BSplinePatches{counterPatches}.FGamma = funcHandle...
            (BSplinePatches{counterPatches}.FGamma,...
            BSplinePatches{counterPatches},...
            NBC{counterPatches}.xiLoadExtension{counterNBC},...
            NBC{counterPatches}.etaLoadExtension{counterNBC},...
            NBC{counterPatches}.loadAmplitude{counterNBC},...
            NBC{counterPatches}.loadDirection{counterNBC},...
            NBC{counterPatches}.isFollower(counterNBC,1),0,...
            BSplinePatches{counterPatches}.int,'outputEnabled');
    end
end

%% Output the initial geometry to be read by GiD
% pathToOutput = '../../outputGiD/isogeometricMembraneAnalysis/';
% writeOutMultipatchBSplineSurface4GiD(BSplinePatches,pathToOutput,caseName);

%% Plot the multipatch geometry together with the parametric coordinates
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,graph);

%% Plot the reference configuration for the multipatch geometry
color = [.4 .4 .4];
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 10e-9;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 1;

% On the load conservativeness for each patch
noPatches = length(BSplinePatches);
propNLinearAnalysis.conservativeLoad = true(noPatches,1);

%% Solve the steady-state nonlinear problem using the Nitsche method
plot_IGANLinear = 'undefined';
[dHat,CPHistory,resHistory,hasConverged,BSplinePatches{1},minElASize] = ...
    solve_IGAMembraneNLinear...
    (BSplinePatches{1},propNLinearAnalysis,solve_LinearSystem,plot_IGANLinear,...
    graph,'outputEnabled');

%% Solve the steady-state nonlinear problem using the Lagrange Multipliers method
% propCoupling.alphaD = max([norm(Dm1) norm(Dm2)])*(1/0.2143)*0;
% propCoupling.alphaR = max([norm(Db1) norm(Db2)])*0;
% propCoupling.intC = intC;
% plot_IGANLinear = 'undefined';
% [dHatLM,CPHistoryLM,resHistoryLM,minElASize] = ...
%     solve_DDMLagrangeMultipliersIGAMembraneMultipatchesNLinear...
%     (BSplinePatches,connections,propCoupling,propNLinearAnalysis,...
%     solve_LinearSystem,plot_IGANLinear,graph,'outputEnabled');

%% Postprocessing
scaling = 1;
graph.index = plot_postprocIGAMembraneNLinear...
    (BSplinePatches{1}{1},dHat*scaling,graph,'outputEnabled');
% graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
%     (BSplinePatches,dHatLM*scaling,graph,'outputEnabled');
% title('Lagrange Multipliers');

%% END OF SCRIPT
