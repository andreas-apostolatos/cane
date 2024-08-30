%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : The benchmark example is a curved beam (modelled as a thin plate) 
%        which is subject into tip shear force. For this example there is 
%        closed form solution in terms of the stress field.
%
% Date : 15.11.2013
%        19.04.2020
%
%% Preamble
clc;
clear;

%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Add linear solvers
addpath('../../equationSystemSolvers/');

% Add functions related to the efficient computation
addpath('../../efficientComputation/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BSplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the Isogeometric Plane Stress formulation
addpath('../../isogeometricThinStructureAnalysis/stiffnessMatrices/',...
        '../../isogeometricThinStructureAnalysis/loads/',...
        '../../isogeometricThinStructureAnalysis/solvers/',...
        '../../isogeometricThinStructureAnalysis/graphicsSinglePatch/',...
        '../../isogeometricThinStructureAnalysis/postprocessing/',...
        '../../isogeometricThinStructureAnalysis/auxiliary/',...
        '../../isogeometricThinStructureAnalysis/errorComputation/', ...
        '../../isogeometricThinStructureAnalysis/graphics/');
    
%% NURBS parameters

% Geometrical parameters
internalRadius = 4;
externalRadius = 5;

% Polynomial degrees
p = 2;
q = 1;

% Knot vectors
Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [0 0
             internalRadius externalRadius;
             internalRadius externalRadius];
         
% y-coordinates
CP(:,:,2) = [internalRadius externalRadius
             internalRadius externalRadius
             0 0];
         
% z-coordinates
CP(:,:,3) = [0 0;
             0 0
             0 0];

% Weights
w = cos(45*pi/180);
CP(:,:,4) = [1 1;
              w w
              1 1];
          
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = false;
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));
for i = 1:numCPs_xi
    for j = 1:numCPs_eta
        if CP(i, j, 4)~=1
            isNURBS = true;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% Material constants 

% Young's modulus
propParameters.E = 1e5;

% Poisson ratio
propParameters.nue = 0.0;

% Thickness of the plate
propParameters.t = 1;

%% UI

% Analysis type
propAnalysis.type = "isogeometricPlateInMembraneActionAnalysis";

% Function handle for the solution of the linear equation system
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points
propInt.type = 'default';
if strcmp(propInt.type,'user')
    propInt.xiNGP = 6;
    propInt.etaNGP = 6;
    propInt.xiNGPForLoad = 6;
    propInt.xetaNGPForLoad = 6;
end

% On the graphics
propGraph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
propGraph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement', 'strain', 'stress'
propGraph.resultant = 'stress';

% Component of the resultant to plot
% .component: 'x', 'y', '2norm','xy', '1Principal', '2Principal'
propGraph.component = 'xy';

% Error computation for this benchmark
% .resultant: 'displacement', 'strain', 'stress'
propError.resultant = 'stress';
% .component: 'x', 'y', 'tensor','xy'
propError.component = 'tensor';
% .xiNGP .eta.NGP
propError.xiNGP = 16;
propError.etaNGP = 16;

% Function handle to the computation of the body force vector
computeBodyForces = @computeConstantVecrticalBodyForceVct;

%% Refinement 

% Degree elevation of the surface
a = 0;
tp = a;  tq = a + 1;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface ...
    (p, q, Xi, Eta, CP, tp, tq, 'outputEnabled');

% Knot insertion on the surface
a = 1/3; % 2
scaling_xi = 25;
scaling_eta = 7;
refXi = ceil(scaling_xi*a);
refEta = ceil(scaling_eta*a);
[Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, refXi, refEta, 'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

% Homogeneous Dirichlet boundary conditions
homDOFs = [];
xiSup = [0 0];
etaSup = [0 1];
dirSupp = 1;
homDOFs = findDofs2D(homDOFs, xiSup, etaSup, dirSupp, CP);
xiSup = [0 0];
etaSup = [0 0];
dirSupp = 2;
homDOFs = findDofs2D(homDOFs, xiSup, etaSup, dirSupp, CP);

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% Weak Dirichlet boundary conditions
weakDBC.noCnd = 0;

% Embedded cables
cables.No = 0;

% load (Neuman boundary conditions)
FAmp = -1;
NBC.noCnd = 1;
NBC.computeLoadVct = {'computeLoadVctLineIGAPlateInMembraneAction'};
NBC.xiLoadExtension = {1};
NBC.etaLoadExtension = {[0 1]};
NBC.loadAmplitude = {FAmp};
NBC.loadDirection = {'x'};
NBC.isFollower(1, 1) = false;
NBC.isTimeDependent(1, 1) = false;

%% Fill up the patch array
BSplinePatch = fillUpPatch ...
    (propAnalysis, p, Xi, q, Eta, CP, isNURBS, propParameters, homDOFs, ...
    inhomDOFs, valuesInhomDOFs, weakDBC, cables, NBC, [], [], [], [], [], ...
    propInt);
BSplinePatches = {BSplinePatch};
numPatches = length(BSplinePatches);

%% Compute the load vector for the visualization of the reference configuration
for iPatches = 1:numPatches
    BSplinePatches{iPatches}.FGamma = ...
        zeros(2*BSplinePatches{iPatches}.noCPs,1);
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

%% Plot the initial configuration
figure(propGraph.index)
plot_referenceConfigurationIGAPlateInMembraneAction(p, q, Xi, Eta, CP, ...
    isNURBS, homDOFs, BSplinePatches{1}.FGamma, 'outputEnabled');
title('Reference configuration for an isogeometric plate in membrane action');
propGraph.index = propGraph.index + 1;

%% Solve the system
[dHat, FComplete, rankD, condK, minEig, minElSize] = ...
    solve_IGAPlateInMembraneActionLinear...
    (propAnalysis, BSplinePatch, 'undefined', solve_LinearSystem, 'outputEnabled');

%% Postprocessing 

% Plot the current configuration and the selected resultant
propGraph.index = plot_postprocIGAPlateInMembraneAction ...
    (p, q, Xi, Eta, CP, isNURBS, homDOFs, propParameters, ...
    BSplinePatches{1}.FGamma, dHat, propGraph, 'outputEnabled');

% Compute the relative error in the L2-norm for the stress resultant
[errorCurvedBeamTipShear, minElArea] =  ...
    computeRelErrorL2CurvedBeamTipShearIGAPlateInMembraneAction ...
    (p, q, Xi, Eta, CP, isNURBS, propParameters, internalRadius, ...
    externalRadius, abs(FAmp), dHat, propError, 'outputEnabled');

%% End of the scrpit