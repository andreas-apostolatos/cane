%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% IGA curved beam under tip shear: convergence study
%
% This script runs a compact h-refinement study for the curved beam
% benchmark and compares the numerical stress field with the analytical
% solution used by main_curvedBeamTipShear.m.

clc;
clear;

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(repoRoot));

%% Geometry and material parameters
internalRadius = 4;
externalRadius = 5;

p0 = 2;
q0 = 1;
Xi0 = [0 0 0 1 1 1];
Eta0 = [0 0 1 1];

CP0(:, :, 1) = [0 0
                internalRadius externalRadius
                internalRadius externalRadius];
CP0(:, :, 2) = [internalRadius externalRadius
                internalRadius externalRadius
                0 0];
CP0(:, :, 3) = [0 0
                0 0
                0 0];
w = cos(45*pi/180);
CP0(:, :, 4) = [1 1
                w w
                1 1];

parameters.E = 1e5;
parameters.nue = 0.0;
parameters.t = 1;

analysis.type = "isogeometricPlateInMembraneActionAnalysis";
solveLinearSystem = @solve_LinearSystemMatlabBackslashSolver;

int.type = 'default';
error.resultant = 'stress';
error.component = 'tensor';
error.xiNGP = 16;
error.etaNGP = 16;

weakDBC.noCnd = 0;
cables.No = 0;

loadAmplitude = -1;

refXiValues = [9 12 15 18 21 24];
refEtaValues = max(1, ceil(refXiValues*7/25));

relativeStressError = zeros(numel(refXiValues), 1);
minimumElementArea = zeros(numel(refXiValues), 1);
numberOfDofs = zeros(numel(refXiValues), 1);

fprintf('IGA curved beam convergence benchmark\n');
fprintf('------------------------------------------------\n');

for iRef = 1:numel(refXiValues)
    p = p0;
    q = q0;
    Xi = Xi0;
    Eta = Eta0;
    CP = CP0;
    
    isNURBS = any(any(CP(:, :, 4) ~= 1));
    
    [Xi, Eta, CP, p, q] = degreeElevateBSplineSurface ...
        (p, q, Xi, Eta, CP, 0, 1, '');
    [Xi, Eta, CP] = knotRefineUniformlyBSplineSurface ...
        (p, Xi, q, Eta, CP, refXiValues(iRef), refEtaValues(iRef), '');
    
    homDOFs = [];
    homDOFs = findDofs2D(homDOFs, [0 0], [0 1], 1, CP);
    homDOFs = findDofs2D(homDOFs, [0 0], [0 0], 2, CP);
    
    NBC.noCnd = 1;
    NBC.computeLoadVct = {'computeLoadVctLineIGAPlateInMembraneAction'};
    NBC.xiLoadExtension = {1};
    NBC.etaLoadExtension = {[0 1]};
    NBC.loadAmplitude = {loadAmplitude};
    NBC.loadDirection = {'x'};
    NBC.isFollower(1, 1) = false;
    NBC.isTimeDependent(1, 1) = false;
    
    patch = fillUpPatch ...
        (analysis, p, Xi, q, Eta, CP, isNURBS, parameters, homDOFs, ...
        [], [], weakDBC, cables, NBC, [], [], [], [], [], int);
    
    [dHat, ~, ~, ~, ~, ~] = solve_IGAPlateInMembraneActionLinear ...
        (analysis, patch, 'undefined', solveLinearSystem, '');
    
    [relativeStressError(iRef), minimumElementArea(iRef)] = ...
        computeRelErrorL2CurvedBeamTipShearIGAPlateInMembraneAction ...
        (p, q, Xi, Eta, CP, isNURBS, parameters, internalRadius, ...
        externalRadius, abs(loadAmplitude), dHat, error, '');
    
    numberOfDofs(iRef) = 2*patch.noCPs;
    fprintf('refXi=%2d refEta=%2d DOFs=%4d minArea=%8.3e relError=%8.3e\n', ...
        refXiValues(iRef), refEtaValues(iRef), numberOfDofs(iRef), ...
        minimumElementArea(iRef), relativeStressError(iRef));
end

figure('Name', 'IGA curved beam stress convergence');
loglog(sqrt(minimumElementArea), relativeStressError, '-bo', ...
    'LineWidth', 2, 'MarkerSize', 7);
grid on;
xlabel('Minimum element size indicator');
ylabel('$\| \sigma - \sigma_h \|_{L^2}/\|\sigma\|_{L^2}$', ...
    'Interpreter', 'latex');
title('Curved beam stress convergence against analytical solution');

figure('Name', 'IGA curved beam stress convergence by DOFs');
loglog(numberOfDofs, relativeStressError, '-bo', ...
    'LineWidth', 2, 'MarkerSize', 7);
grid on;
xlabel('Number of DOFs');
ylabel('$\| \sigma - \sigma_h \|_{L^2}/\|\sigma\|_{L^2}$', ...
    'Interpreter', 'latex');
title('Curved beam stress convergence against analytical solution');

%% End of script
