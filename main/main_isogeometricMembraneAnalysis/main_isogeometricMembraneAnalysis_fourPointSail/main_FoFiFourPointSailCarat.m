%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Reading and analysing the results of form-finding over the
%        four-point sail example from Carat++
%
% Date : 30.01.2017
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
    
%% Load the NURBS parameters from a FoFi file

% Patch 1 :
% _________

% Polynomial orders
p = 1;
q = 1;

% Knot vectors
Xi = [0.0000000000e+00, 0.0000000000e+00, 2.0000000000e+00, 4.0000000000e+00, 6.0000000000e+00, 8.0000000000e+00, 1.0000000000e+01, 1.0000000000e+01];
Eta = [0.0000000000e+00, 0.0000000000e+00, 2.0000000000e+00, 4.0000000000e+00, 6.0000000000e+00, 8.0000000000e+00, 1.0000000000e+01, 1.0000000000e+01];

% Control Point coordinates and weights 
CP(:,:,1) = [-1.0000000000e+01 -7.9172001006e+00 -6.8529900697e+00 -6.8529900697e+00 -7.9172001006e+00 -1.0000000000e+01
             -6.2738927771e+00 -5.2686674089e+00 -4.9065346315e+00 -4.9065346315e+00 -5.2686674089e+00 -6.2738927771e+00
             -2.1382330878e+00 -2.0095961303e+00 -1.8581603860e+00 -1.8581603860e+00 -2.0095961303e+00 -2.1382330878e+00
             2.1382330878e+00 2.0095961303e+00 1.8581603860e+00 1.8581603860e+00 2.0095961303e+00 2.1382330878e+00
             6.2738927771e+00 5.2686674089e+00 4.9065346315e+00 4.9065346315e+00 5.2686674089e+00 6.2738927771e+00
             1.0000000000e+01 7.9172001006e+00 6.8529900697e+00 6.8529900697e+00 7.9172001006e+00 1.0000000000e+01];
          
CP(:,:,2) = [-1.0000000000e+01 -6.2738927771e+00 -2.1382330878e+00 2.1382330878e+00 6.2738927771e+00 1.0000000000e+01
             -7.9172001006e+00 -5.2686674089e+00 -2.0095961303e+00 2.0095961303e+00 5.2686674089e+00 7.9172001006e+00
             -6.8529900697e+00 -4.9065346315e+00 -1.8581603860e+00 1.8581603860e+00 4.9065346315e+00 6.8529900697e+00
             -6.8529900697e+00 -4.9065346315e+00 -1.8581603860e+00 1.8581603860e+00 4.9065346315e+00 6.8529900697e+00 
             -7.9172001006e+00 -5.2686674089e+00 -2.0095961303e+00 2.0095961303e+00 5.2686674089e+00 7.9172001006e+00 
             -1.0000000000e+01 -6.2738927771e+00 -2.1382330878e+00 2.1382330878e+00 6.2738927771e+00 1.0000000000e+01];

CP(:,:,3) = [1.0000000000e+01 7.8656166362e+00 5.9206942888e+00 4.0793057112e+00 2.1343833638e+00 0.0000000000e+00
             7.8656166362e+00 6.6212599721e+00 5.5812750947e+00 4.4187249053e+00 3.3787400279e+00 2.1343833638e+00
             5.9206942888e+00 5.5812750947e+00 5.1959341872e+00 4.8040658128e+00 4.4187249053e+00 4.0793057112e+00
             4.0793057112e+00 4.4187249053e+00 4.8040658128e+00 5.1959341872e+00 5.5812750947e+00 5.9206942888e+00
             2.1343833638e+00 3.3787400279e+00 4.4187249053e+00 5.5812750947e+00 6.6212599721e+00 7.8656166362e+00
             0.0000000000e+00 2.1343833638e+00 4.0793057112e+00 5.9206942888e+00 7.8656166362e+00 1.0000000000e+01];
         
CP(:,:,4) = [1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00
             1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00
             1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00
             1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00
             1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00
             1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00 1.0000000000e+00];
          
% Flag on whether the basis is a B-Spline or a NURBS
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

% Young's modulus
parameters.E = 2.1e9;

% Poisson ratio
parameters.nue = 0.4;

% Thickness
parameters.t = 1e-3;

% Density (used only for dynamics)
parameters.rho = 8050;

% Prestress for the membrane
sigma0 = 1e3;
parameters.prestress.voigtVector = [sigma0
                                    sigma0
                                    0];
                                
% Cable parameters
Radius = 10;
parametersCable.E = parameters.E;
parametersCable.radiusCS = 1e-2;
parametersCable.areaCS = pi*parametersCable.radiusCS^2;
parametersCable.rho = parameters.rho;
parametersCable.prestress = sigma0*Radius;

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
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement','strain','curvature','force','moment','shearForce'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x', 'y','z','2norm','1','2','12','1Principal','2Principal'
graph.component = '2norm';

%% Refinement

% % Degree by which to elevate
% a = 0;
% tp = a;
% tq = a;
% [Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');
% 
% % Number of knots to exist in both directions
% scaling = 0;
% refXi = ceil((Length/Width)*scaling);
% refEta = scaling;
% [Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

% supports (Dirichlet boundary conditions)
homDOFs = [];
xiSup = [0 0];   etaSup = [0 0];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [0 0];   etaSup = [1 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [0 0];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end
xiSup = [1 1];   etaSup = [1 1];
for dirSupp = 1:3
    homDOFs = findDofs3D(homDOFs,xiSup,etaSup,dirSupp,CP);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak homogeneous Dirichlet boundary conditions
weakDBC.noCnd = 4;
% weakDBC.method = 'nitsche';
weakDBC.method = 'penalty';
weakDBC.estimationStabilPrm = true;
weakDBC.alpha = 1e3;
weakDBC.xiExtension = {[0 0] [0 0] [1 1] [1 1]};
weakDBC.etaExtension = {[0 0] [1 1] [0 0] [1 1]};
weakDBC.int.type = 'default';
weakDBC.int.noGPs = 16;

% Embedded cables
cables.No = 4;
cables.xiExtension = {[0 1] [0 1] [0 0] [1 1]};
cables.etaExtension = {[0 0] [1 1] [0 1] [0 1]};
cables.parameters = {parametersCable parametersCable parametersCable parametersCable};
cables.int.type = 'default';
% cables.int.type = 'user';
cables.int.noGPs = 16;

% load (Neuman boundary conditions)
scaling = 0.0;
FAmp = - 1e5*scaling;
xib = [0 1];   etab = [0 1];   dirForce = 'z';
NBC.noCnd = 1;
NBC.xiLoadExtension = {xib};
NBC.etaLoadExtension = {etab};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {dirForce};
NBC.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC.isFollower = false;
NBC.isTimeDependent = false;

%% Create the B-Spline patch array
BSplinePatch = fillUpPatch...
    (analysis,p,Xi,q,Eta,CP,isNURBS,parameters,homDOFs,inhomDOFs,valuesInhomDOFs,...
    weakDBC,cables,NBC,[],[],[],[],[],int);
BSplinePatches = {BSplinePatch};
noPatches = length(BSplinePatches);
connections = [];

%% Compute the load vector for the visualization of the reference configuration
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.FGamma = ...
        zeros(3*BSplinePatches{iPatches}.noCPs,1);
    NBC = BSplinePatches{iPatches}.NBC;
    for iNBC = 1:NBC.noCnd
        funcHandle = str2func(NBC.computeLoadVct{iNBC});
        BSplinePatches{iPatches}.FGamma = funcHandle...
            (BSplinePatches{iPatches}.FGamma,...
            BSplinePatches{iPatches},...
            NBC.xiLoadExtension{iNBC},...
            NBC.etaLoadExtension{iNBC},...
            NBC.loadAmplitude{iNBC},...
            NBC.loadDirection{iNBC},...
            NBC.isFollower(iNBC,1),0,...
            BSplinePatches{iPatches}.int,'outputEnabled');
    end
end

%% Plot the reference configuration for the multipatch geometry before the form-finding analysis
color = [217 218 219]/255;
graph.index = 3;
graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
    (BSplinePatches,connections,color,graph,'outputEnabled');
% az = -40;
% el = 40;
% view(az,el);
% axis off;
% title('');
% camlight left;
% lighting phong;

%% Nonlinear analysis properties

% % number of load steps for the non-linear analysis
% propNLinearAnalysis.method = 'newtonRapshon';
% 
% % number of load steps for the non-linear analysis
% propNLinearAnalysis.noLoadSteps = 1;
% 
% % Assign a tolerance for the Newton iterations
% propNLinearAnalysis.eps = 1e-5;
% 
% % Assign the maximum number of iterations
% propNLinearAnalysis.maxIter = 1000;

%% Solve the system applying nonlinear analysis
% % plot_IGANLinear = @plot_postprocIGAMembraneNLinear;
% plot_IGANLinear = '';
% [dHatNLinear,CPHistory,resHistory,hasConverged,minElArea] = solve_IGAMembraneNLinear...
%     (BSplinePatch,propNLinearAnalysis,solve_LinearSystem,plot_IGANLinear,...
%     graph,'outputEnabled');

%% Postprocessing
% scaling = 1;
% graph.index = plot_postprocIGAMembraneNLinear...
%     (BSplinePatch,scaling*dHatNLinear,graph,'outputEnabled');

%% END
