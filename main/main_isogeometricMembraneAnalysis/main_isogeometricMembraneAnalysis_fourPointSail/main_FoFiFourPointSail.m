%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Form-finding analysis for a four-point sail
%
% Date : 08.02.2020
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

%% CAD modelling via NURBS

% Global variables:
Length = 20;
Width = Length;
Height = Length/2;

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [-Length/2 Length/2
             -Length/2 Length/2];
         
% y-coordinates
CP(:,:,2) = [-Width/2 -Width/2
             Width/2  Width/2];
         
% z-coordinates
CP(:,:,3) = [Height 0    
             0      Height];
       
% Weights
CP(:,:,4) = [1 1
             1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i = 1:nxi
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

% general parameters
EYoung = 8e+8;
nue = .4;
thickness = 1e-3;
sigma0 = 3e+3;
prestress.voigtVector = [sigma0/thickness
                         sigma0/thickness
                         0];
density = 8050;

% Young's modulus
parameters.E = EYoung;

% Poisson ratio
parameters.nue = nue;

% Thickness of the membrane
parameters.t = thickness;

% Density of the membrane
parameters.rho = density;

% Prestress of the membrane
parameters.prestress = prestress;

% Cable :
% _______

parametersCable.E = 1.6e+11;
parametersCable.radiusCS = 12e-3/2;
parametersCable.areaCS = pi*parametersCable.radiusCS^2;
parametersCable.rho = 8050;
parametersCable.prestress = 6e+4/parametersCable.areaCS;

%% GUI

% Case name
caseName = 'SinglePatchFourPointSail';

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

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
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','curvature','force','moment','shearForce'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x','y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

% Define the coupling properties

% Patch 1 :
% _________

Dm = parameters.E*parameters.t/(1 - parameters.nue^2)*...
      [1              parameters.nue 0
       parameters.nue 1              0
       0               0              (1-parameters.nue)/2];

% Assign the penalty factors

% Function handle to writing out the results
% writeOutput = @writeResults4Carat;
% writeOutput = @writeResults4GiD;
writeOutput = 'undefined';

% Postprocessing
propPostproc.resultant = {'displacement'};
propPostproc.computeResultant = {'computeDisplacement'};

%% Refinement

% Select p-refinement level
iPRef = 1; % Coarse : 1 || Fine : 2

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

a = 1;
tp = a;
tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'');

%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

noKnotsXi = 10;
noKnotsEta = 10;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface...
    (p,Xi,q,Eta,CP,noKnotsXi,noKnotsEta,'');

%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs = [];
xisup = [0 0];   etasup = [0 0];
for dir = 1:3
    homDOFs = findDofs3D...
        (homDOFs,xisup,etasup,dir,CP);
end
xisup = [1 1];   etasup = [0 0];
for dir = 1:3
    homDOFs = findDofs3D...
        (homDOFs,xisup,etasup,dir,CP);
end
xisup = [0 0];   etasup = [1 1];
for dir = 1:3
    homDOFs = findDofs3D...
        (homDOFs,xisup,etasup,dir,CP);
end
xisup = [1 1];   etasup = [1 1];
for dir = 1:3
    homDOFs = findDofs3D...
        (homDOFs,xisup,etasup,dir,CP);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables.No = 4;
cables.xiExtension = {[0 1] [0 1] [0 0] [1 1]};
cables.etaExtension = {[0 0] [1 1] [0 1] [0 1]};
cables.parameters = {parametersCable parametersCable parametersCable parametersCable};
cables.int.type = 'default';
% cables1.int.type = 'user';
cables.int.noGPs = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parameter
loadAmplitude = 0;

% Patch 1 :
% _________

FAmp = loadAmplitude;
NBC.noCnd = 1;
xib = [0 1];   etab = [0 1];   dirForce = 'z';
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC.isFollower = false;
NBC.isTimeDependent = true;

%% Fill up the arrays for the patches
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,...
    inhomDOFs,valuesInhomDOFs,weakDBC,cables,NBC,[],[],[],[],[],int);
BSplinePatches = {BSplinePatch};
noPatches = length(BSplinePatches);
connections.No = 0;

%% Plot the distribution of the determinant of the geometrical Jacobian for the multipatch geometry
% figure(graph.index)
% for iPatches = 1:noPatches
%     plot_postprocBSplineSurfaceGeometricJacobian(BSplinePatches{iPatches}.p,...
%         BSplinePatches{iPatches}.q,BSplinePatches{iPatches}.Xi,...
%         BSplinePatches{iPatches}.Eta,BSplinePatches{iPatches}.CP,...
%         BSplinePatches{iPatches}.isNURBS);
%     hold on;
% end
% hold off;
% shading interp;
% colormap('jet');
% camlight left;
% lighting phong;
% colorbar;
% axis equal;
% graph.index = graph.index + 1;

%% Plot the multipatch geometry with the patch numbering
% color = [217 218 219]/255;
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,color,graph);

%% Compute the load vector for the visualization of the reference configuration
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.FGamma = ...
        zeros(3*BSplinePatches{iPatches}.noCPs,1);
%     for counterNBC = 1:NBC{counterPatches}.noCnd
%         funcHandle = str2func(NBC{counterPatches}.computeLoadVct{counterNBC});
%         BSplinePatches{counterPatches}.FGamma = funcHandle...
%             (BSplinePatches{counterPatches}.FGamma,...
%             BSplinePatches{counterPatches},...
%             NBC{counterPatches}.xiLoadExtension{counterNBC},...
%             NBC{counterPatches}.etaLoadExtension{counterNBC},...
%             NBC{counterPatches}.loadAmplitude{counterNBC},...
%             NBC{counterPatches}.loadDirection{counterNBC},...
%             NBC{counterPatches}.isFollower(counterNBC,1),0,...
%             BSplinePatches{counterPatches}.int,'outputEnabled');
%     end
end

%% Plot the reference configuration for the multipatch geometry before the form-finding analysis
color = [217 218 219]/255;
% color = 'none';
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = 260;
el = 30;
view(az,el);
camlight(30,70);
lighting phong;

%% Load a Finite Element solution
% nodes = importdata('../../../preComputedData/isogeometricMembraneAnalysis/4ptSail/referenceFEMSolution/nodes');
% elements = importdata('../../../preComputedData/isogeometricMembraneAnalysis/4ptSail/referenceFEMSolution/elements');
% displacement = importdata('../../../preComputedData/isogeometricMembraneAnalysis/4ptSail/referenceFEMSolution/displacements');
% mesh.elements = elements(:,1:3);
% for iNodes = 1:length(nodes(:,1))
%     [index,~] = find(displacement(iNodes,1) == nodes(:,1));
%     nodes(index,2:4) = nodes(index,2:4) + displacement(iNodes,2:4);
% end
% mesh.nodes = nodes;

%% Plot the form-found geometry using classical Finite Elements
% color = [217 218 219]/255;
% labelsEnabled = false;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,mesh,labelsEnabled,color,graph);
% az = 30;
% el = 45;
% view([az el]);
% camlight(20,40);
% % camlight left;
% lighting phong;

%% Perform a form finding analysis to find the shape of the multipatch membrane with the Nitsche method
propFormFinding.tolerance = 1e-5;
propFormFinding.maxNoIter = 1e2;
[BSplinePatch,CPHistory,resHistory,hasConverged,noIter] = ...
    solve_formFindingIGAMembrane(BSplinePatch,propFormFinding,...
    solve_LinearSystem,'outputEnabled');
BSplinePatches = {BSplinePatch};
% save data_FoFiFourPointSailCoarse.mat;
% return;

%% Compute the relative error from the reference results
% referenceData = importdata('../../../preComputedData/isogeometricMembraneAnalysis/referenceSolution_fourPointSailp3q3Xi50Eta50.mat');
% propReferenceSolution.referenceBSplinePatch = referenceData.BSplinePatch;
% propNewtonRapshon.eps = 1e-9;
% propNewtonRapshon.maxIt = 10;
% propError.noSamplingPoints = 10;
% propError.tolClose = 1e-4;
% propInt.type = 'default';
% [relGeoDomainErrL2,relGeoInterfaceErrL2] = ....
%     computeDomainAndInterfaceErrorInL2NormMembraneFormFiding...
%     (BSplinePatches,connections,propReferenceSolution,propNewtonRapshon,...
%     propError,propInt,'outputEnabled');

%% Plot the reference configuration for the multipatch geometry after the form-finding analysis
color = [217 218 219]/255;
% color = 'none';
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
az = 260;
el = 30;
view(az,el);
camlight(30,70);
lighting phong;

%% END
