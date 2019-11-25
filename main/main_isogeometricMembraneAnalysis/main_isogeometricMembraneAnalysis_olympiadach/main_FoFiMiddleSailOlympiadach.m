%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Form-finding analysis over the middle tent of the Olympiadach
%        modelled with a single patch
%
% Date : 13.02.2016
%
%% Preamble
clear;
clc;

%% Includes

% Add general math functions
addpath('../../../generalMath/');

% Add general auxiliary functions
addpath('../../../auxiliary/');

% Add all linear equation system solvers
addpath('../../../equationSystemSolvers/');

% Add all efficient computation functions
addpath('../../../efficientComputation/');

% Add graphics-related functions for the isogeometric mortar-based mapping
addpath('../../../isogeometricMortarBasedMappingAnalysis/graphics/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../../CAGDKernel/CAGDKernel_graphics/',...
        '../../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the Isogeometric Kirchhoff-Love shell formulation
addpath('../../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../../isogeometricThinStructureAnalysis/graphicsMultipatches/',...
        '../../../isogeometricThinStructureAnalysis/loads/',...
        '../../../isogeometricThinStructureAnalysis/solutionMatricesAndVectors/',...
        '../../../isogeometricThinStructureAnalysis/solvers/',...
        '../../../isogeometricThinStructureAnalysis/metrics/',...
        '../../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../../isogeometricThinStructureAnalysis/BOperatorMatrices/',...
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/');

%% NURBS parameters

% Global variables:
noPnts = 10;
P = zeros(3,noPnts);
P(:,1) = [14.947
          64.456
          39.79];
P(:,2) = [40.82
          59.937
          39.757];
P(:,3) = [57.515
          77.156
          45.661];
P(:,4) = [70.56
          104.72
          47.657];
P(:,5) = [74.976
          126.2
          31.367];
P(:,6) = [54.05
          135.65
          25.293];
P(:,7) = [43.504
          138.44
          25.723];
P(:,8) = [20.445
          140.89
          33.204];
P(:,9) = [12.832
          120.58
          49.194];
P(:,10) = [10.41
           90.228
           46.113];

% Polynomial degrees
p = 1;
q = 1;

% Knot vectors
Xi = [0 0 1 1];
Eta = [0 0 .25 .5 .75 1 1];

% Control Point coordinates and weights
noCPXi = noPnts/5;
noCPEta = noPnts/2;
CP = zeros(noCPXi,noCPEta,4);
noCoord = 3;
for iCoord = 1:noCoord
    CP(:,:,iCoord) = [P(iCoord,1) P(iCoord,10) P(iCoord,9) P(iCoord,8) P(iCoord,7)
                      P(iCoord,2) P(iCoord,3)  P(iCoord,4) P(iCoord,5) P(iCoord,6)];
end
CP(:,:,4) = [1 1 1 1 1
             1 1 1 1 1];

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

%% Material constants

% Membrane :
% ----------

% Young's modulus
parameters.E = 0.0;

% Poisson ratio
parameters.nue = 0.0;

% Thickness
parameters.t = 1e-2;

% Density (used only for dynamics)
parameters.rho = 0.0;

% Prestress for the membrane
sigma0 = 1e2;
parameters.prestress.voigtVector = [sigma0
                                    sigma0
                                    0];

% Cable 1 :
% ---------

% Young's modulus
parametersCable1.E = 0.0;

% Cross sectional area
parametersCable1.areaCS = 5e-2;

% Density
parametersCable1.rho = 0.0;

% Prestress
parametersCable1.prestress = 5e3;

% Cable 2 :
% ---------

% Young's modulus
parametersCable2.E = 0.0;

% Cross sectional area
parametersCable2.areaCS = 5e-2;

% Density
parametersCable2.rho = 0.0;

% Prestress
parametersCable2.prestress = 2e3;

%% GUI

% Analysis type
analysis.type = 'isogeometricMembraneAnalysis';

% Define linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'default';
if strcmp(int.type,'user')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.xetaNGPForLoad = 6;
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
% .component: 'x', 'y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

%% Refinement

% Degree by which to elevate
a = 0; % 2
tp = a;
tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Number of knots to exist in both directions
refXi = 0; % 40
refEta = 2*refXi;
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

% supports (Dirichlet boundary conditions)
homDOFs = [];
xiSup = {[0 0] [0 0] [1 1] [1 1]};
etaSup = {[0 0] [1 1] [0 0] [1 1]};
noSupPnts = length(xiSup);
for iSupPnts = 1:noSupPnts
    for dirSupp = 1:3
        homDOFs = findDofs3D(homDOFs,xiSup{iSupPnts},etaSup{iSupPnts},dirSupp,CP);
    end
end
xiSup = {[1 1] [1 1] [1 1] [0 0] [0 0] [0 0]};
etaSup = {[.25 .25] [.5 .5] [.75 .75] [.25 .25] [.5 .5] [.75 .75]};
noSupPnts = length(xiSup);
for iCoord = 1:noCoord
    for iSupPnts = 1:noSupPnts
        xi = mean(xiSup{iSupPnts});
        eta = mean(etaSup{iSupPnts});
        if xi ~= mean(xiSup{iSupPnts}) || eta ~= mean(etaSup{iSupPnts})
            error('For point constraints the support span vectors must have equal components');
        end
        nxi = length(CP(:,1,1));
        neta = length(CP(1,:,1));
        xiSpan = findKnotSpan(xi,Xi,nxi);
        etaSpan = findKnotSpan(eta,Eta,neta);
        R = computeIGABasisFunctionsAndDerivativesForSurface...
            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,0);
        X = computeCartesianCoordinatesOfAPointOnBSplineSurface...
            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,R);
        indexX = find(X(1,1) == CP(:,:,1));
        indexY = find(X(2,1) == CP(:,:,2));
        indexZ = find(X(3,1) == CP(:,:,3));
        [sX,mX] = size(indexX);
        [sY,mY] = size(indexY);
        [sZ,mZ] = size(indexZ);
        if mX > 1 || mY > 1 || mZ > 1
            error('The second dimension of the array must be 1, check why this is not the case');
        end
        if sX == sY && sX == sZ && sY == sZ
            noMultipleCPs = sX;
        elseif sX ~= sY || sX ~= sZ || sY ~= sZ
            error('No interpolation point has been found')
        end
        for iMultipleCPs = 1:noMultipleCPs
            indexXi = indexX(iMultipleCPs,1);
            r = length(homDOFs) + 1;
            homDOFs(r) = noCoord*indexXi - (noCoord - iCoord);
        end
    end
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables.No = 10;
cables.xiExtension = {[0 1] [1 1] [1 1] [1 1] [1 1] [0 1] [0 0] [0 0] [0 0] [0 0]};
cables.etaExtension = {[0 0] [0 .25] [.25 .5] [.5 .75] [.75 1] [1 1] [.75 1] [.5 .75] [.25 .5] [0 .25]};
cables.parameters = {parametersCable1 parametersCable1 parametersCable2 parametersCable1 ...
    parametersCable1 parametersCable1 parametersCable1 parametersCable1 parametersCable2 ...
    parametersCable1};
cables.int.type = 'default'; % 'user'
cables.int.noGPs = 16;

% load (Neuman boundary conditions)
NBC.noCnd = 0;

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,valuesInhomDOFs,...
    weakDBC,cables,NBC,[],[],[],[],[],int);
BSplinePatches = {BSplinePatch};
noPatches = length(BSplinePatches);
connections = [];

%% Plot the distribution of the determinant of the geometrical Jacobian
% figure(graph.index)
% plot_postprocBSplineSurfaceGeometricJacobian(p,q,Xi,Eta,CP,isNURBS);
% shading interp;
% colormap('jet');
% camlight left;
% lighting phong;
% colorbar;
% axis equal;
% graph.index = graph.index + 1;

%% Plot the geometry with the parametrization
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%    (BSplinePatches,msh,labelsEnabled,graph);

%% Compute the load vector for the visualization of the reference configuration
% for iPatches = 1:noPatches
%     BSplinePatches{iPatches}.FGamma = ...
%         zeros(3*BSplinePatches{iPatches}.noCPs,1);
% %     for counterNBC = 1:NBC{counterPatches}.noCnd
% %         funcHandle = str2func(NBC{counterPatches}.computeLoadVct{counterNBC});
% %         BSplinePatches{counterPatches}.FGamma = funcHandle...
% %             (BSplinePatches{counterPatches}.FGamma,...
% %             BSplinePatches{counterPatches},...
% %             NBC{counterPatches}.xiLoadExtension{counterNBC},...
% %             NBC{counterPatches}.etaLoadExtension{counterNBC},...
% %             NBC{counterPatches}.loadAmplitude{counterNBC},...
% %             NBC{counterPatches}.loadDirection{counterNBC},...
% %             NBC{counterPatches}.isFollower(counterNBC,1),0,...
% %             BSplinePatches{counterPatches}.int,'outputEnabled');
% %     end
% end

%% Plot the reference configuration for the multipatch geometry before the form-finding analysis
% color = [217 218 219]/255;
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');
% % az = -40;
% % el = 40;
% % view(az,el);
% % axis off;
% % title('');
% % camlight left;
% % lighting phong;

%% Perform a form finding analysis to find the shape of the membrane
propFormFinding.tolerance = 1e-3;
propFormFinding.maxNoIter = 16;
[BSplinePatch,CPHistory,resHistory,hasConverged,noIter] = ...
    solve_formFindingIGAMembrane(BSplinePatch,propFormFinding,...
    solve_LinearSystem,'outputEnabled');
BSplinePatches = {BSplinePatch};
save data_FoFiMiddleSailOlympiadachp0q0Xi40Eta80.mat;
return;

%% Compute the load vector for the visualization of the reference configuration
FGamma = [];
for iNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{iNBC});
    FGamma = funcHandle(FGamma,BSplinePatch,NBC.xiLoadExtension{iNBC},...
                    NBC.etaLoadExtension{iNBC},...
                    NBC.loadAmplitude{iNBC},...
                    NBC.loadDirection{iNBC},...
                    NBC.isFollower(iNBC,1),0,int,'outputEnabled');
end

%% Compute the load vector for the visualization of the reference configuration
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.FGamma = ...
        zeros(3*BSplinePatches{iPatches}.noCPs,1);
%     for iNBC = 1:NBC{iPatches}.noCnd
%         funcHandle = str2func(NBC{iPatches}.computeLoadVct{iNBC});
%         BSplinePatches{iPatches}.FGamma = funcHandle...
%             (BSplinePatches{iPatches}.FGamma,...
%             BSplinePatches{iPatches},...
%             NBC{iPatches}.xiLoadExtension{iNBC},...
%             NBC{iPatches}.etaLoadExtension{iNBC},...
%             NBC{iPatches}.loadAmplitude{iNBC},...
%             NBC{iPatches}.loadDirection{iNBC},...
%             NBC{iPatches}.isFollower(iNBC,1),0,...
%             BSplinePatches{iPatches}.int,'outputEnabled');
%     end
end

%% Plot the reference configuration for the multipatch geometry after the form-finding analysis
color = [217 218 219]/255;
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
% az = -40;
% el = 40;
% view(az,el);
% axis off;
% title('');
% camlight left;
% lighting phong;

%% End
