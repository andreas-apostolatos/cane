%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
% 
% Task : Modal analysis for a form-found four-point sail, whose geometry is
%        directly read from a file.
%
% Date : 13.11.2016
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

addpath('../../../FEMPlateInMembraneActionAnalysis/graphics/');

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
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionKLShell/',...
        '../../../isogeometricThinStructureAnalysis/nitscheDecompositionMembrane/',...
        '../../../isogeometricThinStructureAnalysis/errorComputation/',...
        '../../../isogeometricThinStructureAnalysis/precomputedData/',...
        '../../../isogeometricThinStructureAnalysis/output/',...
        '../../../isogeometricThinStructureAnalysis/transientAnalysis/',...
        '../../../isogeometricThinStructureAnalysis/initialConditions/',...
        '../../../isogeometricThinStructureAnalysis/weakDBCMembrane/',...
        '../../../isogeometricThinStructureAnalysis/formFindingAnalysis/');
    
% Add functions related to the visualization of the multipatch geometry
% related to the isogeometric mortar-based mapping
addpath('../../../isogeometricMortarBasedMappingAnalysis/graphics/');

%% Read the geometry from the results of a form finding analysis

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
noPatches = 24; % 24

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

    % Patch 9 :
    % _________

    % Polynomial orders
    p9 = 2;
    q9 = 2;

    % Knot vectors
    Xi9 = [0 0 0 1 1 1];
    Eta9 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP9(:,:,1) = [Radius + radius Radius + radius Radius
                  Radius + radius Radius + radius Radius
                  0               0               0];

    % y-coordinates
    CP9(:,:,2) = 2*radius + [0 radius radius
                             0 radius radius
                             0 radius radius];

    % z-coordinates
    CP9(:,:,3) = [0               0               0
                  Radius + radius Radius + radius Radius
                  Radius + radius Radius + radius Radius];

    % weights
    CP9(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS9 = 0;
    nxi9 = length(CP9(:,1,1));
    neta9 = length(CP9(1,:,1));
    for i= 1:nxi9
        for j = 1:neta9
            if CP9(i,j,4)~=1
                isNURBS9 = 1;
                break;
            end
        end
        if isNURBS9
            break;
        end
    end

    % Patch 10 :
    % __________

    % Polynomial orders
    p10 = 2;
    q10 = 2;

    % Knot vectors
    Xi10 = [0 0 0 1 1 1];
    Eta10 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP10(:,:,1) = [Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius
                   0               0               0];

    % y-coordinates
    CP10(:,:,2) = 2*radius - [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP10(:,:,3) = [0               0               0
                   Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius];

    % weights
    CP10(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS10 = 0;
    nxi10 = length(CP10(:,1,1));
    neta10 = length(CP10(1,:,1));
    for i = 1:nxi10
        for j = 1:neta10
            if CP10(i,j,4)~=1
                isNURBS10 = 1;
                break;
            end
        end
        if isNURBS10
            break;
        end
    end

    % Patch 11 :
    % __________

    % Polynomial orders
    p11 = 2;
    q11 = 2;

    % Knot vectors
    Xi11 = [0 0 0 1 1 1];
    Eta11 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP11(:,:,1) = [Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius
                   0               0               0];

    % y-coordinates
    CP11(:,:,2) = 2*radius - [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP11(:,:,3) = [0               0               0
                   Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius];

    % weights
    CP11(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS11 = 0;
    nxi11 = length(CP11(:,1,1));
    neta11 = length(CP11(1,:,1));
    for i = 1:nxi11
        for j = 1:neta11
            if CP11(i,j,4)~=1
                isNURBS11 = 1;
                break;
            end
        end
        if isNURBS11
            break;
        end
    end

    % Patch 12 :
    % __________

    % Polynomial orders
    p12 = 2;
    q12 = 2;

    % Knot vectors
    Xi12 = [0 0 0 1 1 1];
    Eta12 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP12(:,:,1) = [Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius
                   0               0               0];

    % y-coordinates
    CP12(:,:,2) = 2*radius + [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP12(:,:,3) = [0               0               0
                   Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius];

    % weights
    CP12(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS12 = 0;
    nxi12 = length(CP12(:,1,1));
    neta12 = length(CP12(1,:,1));
    for i = 1:nxi12
        for j = 1:neta12
            if CP12(i,j,4)~=1
                isNURBS12 = 1;
                break;
            end
        end
        if isNURBS12
            break;
        end
    end

    % Patch 13 :
    % __________

    % Polynomial orders
    p13 = 2;
    q13 = 2;

    % Knot vectors
    Xi13 = [0 0 0 1 1 1];
    Eta13 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP13(:,:,1) = -[Radius + radius Radius + radius Radius
                    Radius + radius Radius + radius Radius
                    0               0               0];

    % y-coordinates
    CP13(:,:,2) = 2*radius + [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP13(:,:,3) = [0               0               0
                   Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius];

    % weights
    CP13(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS13 = 0;
    nxi13 = length(CP13(:,1,1));
    neta13 = length(CP13(1,:,1));
    for i= 1:nxi13
        for j = 1:neta13
            if CP13(i,j,4)~=1
                isNURBS13 = 1;
                break;
            end
        end
        if isNURBS13
            break;
        end
    end

    % Patch 14 :
    % __________

    % Polynomial orders
    p14 = 2;
    q14 = 2;

    % Knot vectors
    Xi14 = [0 0 0 1 1 1];
    Eta14 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP14(:,:,1) = -[Radius + radius Radius + radius Radius
                    Radius + radius Radius + radius Radius
                    0               0               0];

    % y-coordinates
    CP14(:,:,2) = 2*radius - [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP14(:,:,3) = [0               0               0
                   Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius];

    % weights
    CP14(:,:,4) = [1 1 2
                  1 1 2
                  2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS14 = 0;
    nxi14 = length(CP14(:,1,1));
    neta14 = length(CP14(1,:,1));
    for i = 1:nxi14
        for j = 1:neta14
            if CP14(i,j,4)~=1
                isNURBS14 = 1;
                break;
            end
        end
        if isNURBS14
            break;
        end
    end

    % Patch 15 :
    % __________

    % Polynomial orders
    p15 = 2;
    q15 = 2;

    % Knot vectors
    Xi15 = [0 0 0 1 1 1];
    Eta15 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP15(:,:,1) = -[Radius - radius Radius - radius Radius
                    Radius - radius Radius - radius Radius
                    0               0               0];

    % y-coordinates
    CP15(:,:,2) = 2*radius - [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP15(:,:,3) = [0               0               0
                   Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius];

    % weights
    CP15(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS15 = 0;
    nxi15 = length(CP15(:,1,1));
    neta15 = length(CP15(1,:,1));
    for i = 1:nxi15
        for j = 1:neta15
            if CP15(i,j,4)~=1
                isNURBS15 = 1;
                break;
            end
        end
        if isNURBS15
            break;
        end
    end

    % Patch 16 :
    % __________

    % Polynomial orders
    p16 = 2;
    q16 = 2;

    % Knot vectors
    Xi16 = [0 0 0 1 1 1];
    Eta16 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP16(:,:,1) = -[Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius
                   0               0               0];

    % y-coordinates
    CP16(:,:,2) = 2*radius + [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP16(:,:,3) = [0               0               0
                  Radius - radius Radius - radius Radius
                  Radius - radius Radius - radius Radius];

    % weights
    CP16(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS16 = 0;
    nxi16 = length(CP16(:,1,1));
    neta16 = length(CP16(1,:,1));
    for i = 1:nxi16
        for j = 1:neta16
            if CP16(i,j,4)~=1
                isNURBS16 = 1;
                break;
            end
        end
        if isNURBS16
            break;
        end
    end

    % Patch 17 :
    % __________

    % Polynomial orders
    p17 = 2;
    q17 = 2;

    % Knot vectors
    Xi17 = [0 0 0 1 1 1];
    Eta17 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP17(:,:,1) = [Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius
                   0               0               0];

    % y-coordinates
    CP17(:,:,2) = 4*radius + [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP17(:,:,3) = [0               0               0
                   Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius];

    % weights
    CP17(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS17 = 0;
    nxi17 = length(CP17(:,1,1));
    neta17 = length(CP17(1,:,1));
    for i= 1:nxi17
        for j = 1:neta17
            if CP17(i,j,4)~=1
                isNURBS17 = 1;
                break;
            end
        end
        if isNURBS17
            break;
        end
    end

    % Patch 18 :
    % __________

    % Polynomial orders
    p18 = 2;
    q18 = 2;

    % Knot vectors
    Xi18 = [0 0 0 1 1 1];
    Eta18 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP18(:,:,1) = [Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius
                   0               0               0];

    % y-coordinates
    CP18(:,:,2) = 4*radius - [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP18(:,:,3) = [0               0               0
                   Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius];

    % weights
    CP18(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS18 = 0;
    nxi18 = length(CP18(:,1,1));
    neta18 = length(CP18(1,:,1));
    for i = 1:nxi18
        for j = 1:neta18
            if CP18(i,j,4)~=1
                isNURBS18 = 1;
                break;
            end
        end
        if isNURBS18
            break;
        end
    end

    % Patch 19 :
    % __________

    % Polynomial orders
    p19 = 2;
    q19 = 2;

    % Knot vectors
    Xi19 = [0 0 0 1 1 1];
    Eta19 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP19(:,:,1) = [Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius
                   0               0               0];

    % y-coordinates
    CP19(:,:,2) = 4*radius - [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP19(:,:,3) = [0               0               0
                   Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius];

    % weights
    CP19(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS19 = 0;
    nxi19 = length(CP19(:,1,1));
    neta19 = length(CP19(1,:,1));
    for i = 1:nxi19
        for j = 1:neta19
            if CP19(i,j,4)~=1
                isNURBS19 = 1;
                break;
            end
        end
        if isNURBS19
            break;
        end
    end

    % Patch 20 :
    % __________

    % Polynomial orders
    p20 = 2;
    q20 = 2;

    % Knot vectors
    Xi20 = [0 0 0 1 1 1];
    Eta20 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP20(:,:,1) = [Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius
                   0               0               0];

    % y-coordinates
    CP20(:,:,2) = 4*radius + [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP20(:,:,3) = [0               0               0
                   Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius];

    % weights
    CP20(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS20 = 0;
    nxi20 = length(CP20(:,1,1));
    neta20 = length(CP20(1,:,1));
    for i = 1:nxi20
        for j = 1:neta20
            if CP20(i,j,4)~=1
                isNURBS20 = 1;
                break;
            end
        end
        if isNURBS20
            break;
        end
    end

    % Patch 21 :
    % __________

    % Polynomial orders
    p21 = 2;
    q21 = 2;

    % Knot vectors
    Xi21 = [0 0 0 1 1 1];
    Eta21 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP21(:,:,1) = -[Radius + radius Radius + radius Radius
                    Radius + radius Radius + radius Radius
                    0               0               0];

    % y-coordinates
    CP21(:,:,2) = 4*radius + [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP21(:,:,3) = [0               0               0
                   Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius];

    % weights
    CP21(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS21 = 0;
    nxi21 = length(CP21(:,1,1));
    neta21 = length(CP21(1,:,1));
    for i = 1:nxi21
        for j = 1:neta21
            if CP21(i,j,4)~=1
                isNURBS21 = 1;
                break;
            end
        end
        if isNURBS21
            break;
        end
    end

    % Patch 22 :
    % __________

    % Polynomial orders
    p22 = 2;
    q22 = 2;

    % Knot vectors
    Xi22 = [0 0 0 1 1 1];
    Eta22 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP22(:,:,1) = -[Radius + radius Radius + radius Radius
                    Radius + radius Radius + radius Radius
                    0               0               0];

    % y-coordinates
    CP22(:,:,2) = 4*radius - [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP22(:,:,3) = [0               0               0
                   Radius + radius Radius + radius Radius
                   Radius + radius Radius + radius Radius];

    % weights
    CP22(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS22 = 0;
    nxi22 = length(CP22(:,1,1));
    neta22 = length(CP22(1,:,1));
    for i = 1:nxi22
        for j = 1:neta22
            if CP22(i,j,4)~=1
                isNURBS22 = 1;
                break;
            end
        end
        if isNURBS22
            break;
        end
    end

    % Patch 23 :
    % __________

    % Polynomial orders
    p23 = 2;
    q23 = 2;

    % Knot vectors
    Xi23 = [0 0 0 1 1 1];
    Eta23 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP23(:,:,1) = -[Radius - radius Radius - radius Radius
                    Radius - radius Radius - radius Radius
                    0               0               0];

    % y-coordinates
    CP23(:,:,2) = 4*radius - [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP23(:,:,3) = [0               0               0
                   Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius];

    % weights
    CP23(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS23 = 0;
    nxi23 = length(CP23(:,1,1));
    neta23 = length(CP23(1,:,1));
    for i = 1:nxi23
        for j = 1:neta23
            if CP23(i,j,4)~=1
                isNURBS23 = 1;
                break;
            end
        end
        if isNURBS23
            break;
        end
    end

    % Patch 24 :
    % __________

    % Polynomial orders
    p24 = 2;
    q24 = 2;

    % Knot vectors
    Xi24 = [0 0 0 1 1 1];
    Eta24 = [0 0 0 1 1 1];

    % Control Point coordinates and weights

    % x-coordinates
    CP24(:,:,1) = -[Radius - radius Radius - radius Radius
                    Radius - radius Radius - radius Radius
                    0               0               0];

    % y-coordinates
    CP24(:,:,2) = 4*radius + [0 radius radius
                              0 radius radius
                              0 radius radius];

    % z-coordinates
    CP24(:,:,3) = [0               0               0
                   Radius - radius Radius - radius Radius
                   Radius - radius Radius - radius Radius];

    % weights
    CP24(:,:,4) = [1 1 2
                   1 1 2
                   2 2 4];

    % Flag on whether the basis is a B-Spline or a NURBS
    isNURBS24 = 0;
    nxi24 = length(CP24(:,1,1));
    neta24 = length(CP24(1,:,1));
    for i = 1:nxi24
        for j = 1:neta24
            if CP24(i,j,4)~=1
                isNURBS24 = 1;
                break;
            end
        end
        if isNURBS24
            break;
        end
    end
end

%% Material constants

% Membrane :
% ----------

% Young's modulus
parameters.E = 3.1e8;

% Poisson ratio
parameters.nue = 3e-1;

% Thickness
parameters.t = 6e-4;

% Density (used only for dynamics)
parameters.rho = 1250;

% Define inner pressure and second normal component of the stress
FAmp = 1e3;
sigma0 = abs(FAmp)*radius/(2*parameters.t);

% Prestress for the membrane
% parameters.prestress.computeParametricCoordinates = ...
%     @(X) [atan(X(3,1)/X(1,1))
%           real(acos((sqrt(X(1,1)^2 + X(3,1)^2) - Radius)/radius))];
parameters.prestress.computeParametricCoordinates = ...
    @(X) [atan(X(3,1)/X(1,1))
          real(acos((sqrt(X(1,1)^2 + X(3,1)^2) - Radius)/radius)) + pi/2];

% baseVctTheta = @(theta1,theta2) [-(Radius + radius*cos(theta2))*sin(theta1)
%                                  0
%                                  (Radius + radius*cos(theta2))*cos(theta1)];
% baseVctPhi = @(theta1,theta2) [(Radius - radius*sin(theta2))*cos(theta1)
%                                radius*cos(theta2)
%                                (Radius - radius*sin(theta2))*sin(theta1)];
% parameters.prestress.computeBaseVectors = @(theta1,theta2) [baseVctTheta(theta1,theta2) baseVctPhi(theta1,theta2)];

parameters.prestress.voigtVector = ...
    @(theta)sigma0*[1
                    ((2*Radius + radius*sin(theta(2,1)))/(Radius + radius*sin(theta(2,1))))
                    0];

%% GUI

% Case name
caseName = 'modalAnalysisDDMInflatedTubes';

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

%% Refinement

% Mesh size
meshSize = 'fine'; % 'dummy', 'coarse', 'fine'
if ~strcmp(meshSize,'dummy') && ~strcmp(meshSize,'coarse') && ~strcmp(meshSize,'fine')
   error('Choose between a dummy, a coarse and a fine setting of the isogeometric discretization');
end

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%
for iPRef = 1
    % Patch 1 :
    % _________

    if strcmp(meshSize,'coarse')
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

    if strcmp(meshSize,'coarse')
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

    if strcmp(meshSize,'coarse')
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

    if strcmp(meshSize,'coarse')
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

    if strcmp(meshSize,'coarse')
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

    if strcmp(meshSize,'coarse')
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

    if strcmp(meshSize,'coarse')
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

    if strcmp(meshSize,'coarse')
        a8 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a8 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a8 = 0;
    end
    tp8 = a8;
    tq8 = a8;

    % Patch 9 :
    % _________

    if strcmp(meshSize,'coarse')
        a9 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a9 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a9 = 0;
    end
    tp9 = a9;
    tq9 = a9;

    % Patch 10 :
    % __________

    if strcmp(meshSize,'coarse')
        a10 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a10 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a10 = 0;
    end
    tp10 = a10;
    tq10 = a10;

    % Patch 11 :
    % __________

    if strcmp(meshSize,'coarse')
        a11 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a11 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a11 = 0;
    end
    tp11 = a11;
    tq11 = a11;

    % Patch 12 :
    % __________

    if strcmp(meshSize,'coarse')
        a12 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a12 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a12 = 0;
    end
    tp12 = a12;
    tq12 = a12;

    % Patch 13 :
    % __________

    if strcmp(meshSize,'coarse')
        a13 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a13 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a13 = 0;
    end
    tp13 = a13;
    tq13 = a13;

    % Patch 14 :
    % __________

    if strcmp(meshSize,'coarse')
        a14 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a14 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a14 = 0;
    end
    tp14 = a14;
    tq14 = a14;

    % Patch 15 :
    % __________

    if strcmp(meshSize,'coarse')
        a15 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a15 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a15 = 0;
    end
    tp15 = a15;
    tq15 = a15;

    % Patch 16 :
    % __________

    if strcmp(meshSize,'coarse')
        a16 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a16 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a16 = 0;
    end
    tp16 = a16;
    tq16 = a16;

    % Patch 17 :
    % __________

    if strcmp(meshSize,'coarse')
        a17 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a17 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a17 = 0;
    end
    tp17 = a17;
    tq17 = a17;

    % Patch 18 :
    % __________

    if strcmp(meshSize,'coarse')
        a18 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a18 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a18 = 0;
    end
    tp18 = a18;
    tq18 = a18;

    % Patch 19 :
    % __________

    if strcmp(meshSize,'coarse')
        a19 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a19 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a19 = 0;
    end
    tp19 = a19;
    tq19 = a19;

    % Patch 20 :
    % __________

    if strcmp(meshSize,'coarse')
        a20 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a20 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a20 = 0;
    end
    tp20 = a20;
    tq20 = a20;

    % Patch 21 :
    % __________

    if strcmp(meshSize,'coarse')
        a21 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a21 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a21 = 0;
    end
    tp21 = a21;
    tq21 = a21;

    % Patch 22 :
    % __________

    if strcmp(meshSize,'coarse')
        a22 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a22 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a22 = 0;
    end
    tp22 = a22;
    tq22 = a22;

    % Patch 23 :
    % _________

    if strcmp(meshSize,'coarse')
        a23 = 0; % 1
    elseif strcmp(meshSize,'fine')
        a23 = 1; % 1
    elseif strcmp(meshSize,'dummy')
        a23 = 0;
    end
    tp23 = a23;
    tq23 = a23;

    % Patch 24 :
    % _________

    if strcmp(meshSize,'coarse')
        a24 = 0; % 0
    elseif strcmp(meshSize,'fine')
        a24 = 0; % 0
    elseif strcmp(meshSize,'dummy')
        a24 = 0;
    end
    tp24 = a24;
    tq24 = a24;

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
        a1 = 13; % 13
    elseif strcmp(meshSize,'fine')
        a1 = 27; % 27
    elseif strcmp(meshSize,'dummy')
        a1 = 0; 
    end

    % Patch 2 :
    % _________

    if strcmp(meshSize,'coarse')
        a2 = 8; % 8
    elseif strcmp(meshSize,'fine')
        a2 = 16; % 16
    elseif strcmp(meshSize,'dummy')
        a2 = 0; 
    end

    % Patch 3 :
    % _________

    if strcmp(meshSize,'coarse')
        a3 = 12; % 12
    elseif strcmp(meshSize,'fine')
        a3 = 25; % 25
    elseif strcmp(meshSize,'dummy')
        a3 = 0; 
    end

    % Patch 4 :
    % _________

    if strcmp(meshSize,'coarse')
        a4 = 9; % 9
    elseif strcmp(meshSize,'fine')
        a4 = 18; % 18
    elseif strcmp(meshSize,'dummy')
        a4 = 0; 
    end

    % Patch 5 :
    % _________

    if strcmp(meshSize,'coarse')
        a5 = 8; % 8
    elseif strcmp(meshSize,'fine')
        a5 = 16; % 16
    elseif strcmp(meshSize,'dummy')
        a5 = 0; 
    end

    % Patch 6 :
    % _________

    if strcmp(meshSize,'coarse')
        a6 = 12; % 12
    elseif strcmp(meshSize,'fine')
        a6 = 25; % 25
    elseif strcmp(meshSize,'dummy')
        a6 = 0; 
    end

    % Patch 7 :
    % _________

    if strcmp(meshSize,'coarse')
        a7 = 9; % 9
    elseif strcmp(meshSize,'fine')
        a7 = 18; % 18
    elseif strcmp(meshSize,'dummy')
        a7 = 0; 
    end

    % Patch 8 :
    % _________

    if strcmp(meshSize,'coarse')
        a8 = 13; % 13
    elseif strcmp(meshSize,'fine')
        a8 = 27; % 27
    elseif strcmp(meshSize,'dummy')
        a8 = 0; 
    end

    % Patch 9 :
    % _________

    if strcmp(meshSize,'coarse')
        a9 = 14; % 14
    elseif strcmp(meshSize,'fine')
        a9 = 28; % 28
    elseif strcmp(meshSize,'dummy')
        a9 = 0; 
    end

    % Patch 10 :
    % __________

    if strcmp(meshSize,'coarse')
        a10 = 8; % 8
    elseif strcmp(meshSize,'fine')
        a10 = 17; % 17
    elseif strcmp(meshSize,'dummy')
        a10 = 0; 
    end

    % Patch 11 :
    % __________

    if strcmp(meshSize,'coarse')
        a11 = 13; % 13
    elseif strcmp(meshSize,'fine')
        a11 = 26; % 26
    elseif strcmp(meshSize,'dummy')
        a11 = 0; 
    end

    % Patch 12 :
    % __________

    if strcmp(meshSize,'coarse')
        a12 = 7; % 7
    elseif strcmp(meshSize,'fine')
        a12 = 15; % 15
    elseif strcmp(meshSize,'dummy')
        a12 = 0;
    end

    % Patch 13 :
    % __________

    if strcmp(meshSize,'coarse')
        a13 = 8; % 8
    elseif strcmp(meshSize,'fine')
        a13 = 17; % 17
    elseif strcmp(meshSize,'dummy')
        a13 = 0;
    end

    % Patch 14 :
    % __________

    if strcmp(meshSize,'coarse')
        a14 = 14; % 14
    elseif strcmp(meshSize,'fine')
        a14 = 28; % 28
    elseif strcmp(meshSize,'dummy')
        a14 = 0; 
    end

    % Patch 15 :
    % __________

    if strcmp(meshSize,'coarse')
        a15 = 7; % 7
    elseif strcmp(meshSize,'fine')
        a15 = 15; % 15
    elseif strcmp(meshSize,'dummy')
        a15 = 0; 
    end

    % Patch 16 :
    % __________

    if strcmp(meshSize,'coarse')
        a16 = 13; % 13
    elseif strcmp(meshSize,'fine')
        a16 = 27; % 27
    elseif strcmp(meshSize,'dummy')
        a16 = 0; 
    end

    % Patch 17 :
    % __________

    if strcmp(meshSize,'coarse')
        a17 = 14; % 14
    elseif strcmp(meshSize,'fine')
        a17 = 28; % 28
    elseif strcmp(meshSize,'dummy')
        a17 = 0; 
    end

    % Patch 18 :
    % __________

    if strcmp(meshSize,'coarse')
        a18 = 7; % 7
    elseif strcmp(meshSize,'fine')
        a18 = 14; % 14
    elseif strcmp(meshSize,'dummy')
        a18 = 0; 
    end

    % Patch 19 :
    % __________

    if strcmp(meshSize,'coarse')
        a19 = 12; % 12
    elseif strcmp(meshSize,'fine')
        a19 = 24; % 24
    elseif strcmp(meshSize,'dummy')
        a19 = 0; 
    end

    % Patch 20 :
    % __________

    if strcmp(meshSize,'coarse')
        a20 = 9; % 9
    elseif strcmp(meshSize,'fine')
        a20 = 19; % 19
    elseif strcmp(meshSize,'dummy')
        a20 = 0; 
    end

    % Patch 21 :
    % __________

    if strcmp(meshSize,'coarse')
        a21 = 7; % 7
    elseif strcmp(meshSize,'fine')
        a21 = 14; % 14
    elseif strcmp(meshSize,'dummy')
        a21 = 0; 
    end

    % Patch 22 :
    % __________

    if strcmp(meshSize,'coarse')
        a22 = 12; % 12
    elseif strcmp(meshSize,'fine')
        a22 = 24; % 24
    elseif strcmp(meshSize,'dummy')
        a22 = 0; 
    end


    % Patch 23 :
    % __________

    if strcmp(meshSize,'coarse')
        a23 = 9; % 9
    elseif strcmp(meshSize,'fine')
        a23 = 19; % 19
    elseif strcmp(meshSize,'dummy')
        a23 = 0; 
    end

    % Patch 24 :
    % __________

    if strcmp(meshSize,'coarse')
        a24 = 14; % 14
    elseif strcmp(meshSize,'fine')
        a24 = 28; % 28
    elseif strcmp(meshSize,'dummy')
        a24 = 0; 
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

% Description of the Dirichlet boundary
xiExtension = [0 0];
etaExtension = [0 1];

% Weak homogeneous Dirichlet boundary conditions
if strcmp(method,'LagrangeMultipliers') || strcmp(method,'AugmentedLagrangeMultipliers')
    weakDBC.noCnd = 0;    
else
    weakDBC.noCnd = 1;
    weakDBC.xiExtension = {xiExtension};
    weakDBC.etaExtension = {etaExtension};
    weakDBC.imposedMotion = {@(x,y,z,t) [0; 0; 0]};
    weakDBC.int.type = 'default';
    weakDBC.int.noGPs = 16;
end

% Loop over the patches to assign the strong and weak Dirichlet boundary
% conditions
for iPatches = 1:noPatches
    if strcmp(method,'LagrangeMultipliers') || strcmp(method,'AugmentedLagrangeMultipliers')
        homDOFs = [];
        CP = eval(['CP' num2str(iPatches)]);
        for iCnd = 1:length(xiExtension)
            for iDOFs = 1:length(etaExtension)
                for iCoord = 1:3
                    homDOFs = findDofs3D(homDOFs,xiExtension,etaExtension,iCoord,CP);
                end
            end
        end
    end
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
    scaling = 0; % 2e1
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
    NBC2.noCnd = 2;
    xib2 = [0 1];   etab2 = [0 1];   dirForce2 = 'normal';
    NBC2.xiLoadExtension = {xib2 xib2};
    NBC2.etaLoadExtension = {etab2 etab2};
    NBC2.loadAmplitude = {FAmp2 wind};
    NBC2.loadDirection = {dirForce2 'y'};
    NBC2.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
    NBC2.isFollower = [true; false];
    NBC2.isTimeDependent = [false; false];

    % Patch 3 :
    % _________

    FAmp3 = loadAmplitude;
    NBC3.noCnd = 2;
    xib3 = [0 1];   etab3 = [0 1];   dirForce3 = 'normal';
    NBC3.xiLoadExtension = {xib3 xib3};
    NBC3.etaLoadExtension = {etab3 etab3};
    NBC3.loadAmplitude = {FAmp3 wind};
    NBC3.loadDirection = {dirForce3 'y'};
    NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
    NBC3.isFollower = [true; false];
    NBC3.isTimeDependent = [false; false];

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
    NBC6.noCnd = 2;
    xib6 = [0 1];   etab6 = [0 1];   dirForce6 = 'normal';
    NBC6.xiLoadExtension = {xib6 xib6};
    NBC6.etaLoadExtension = {etab6 etab6};
    NBC6.loadAmplitude = {FAmp6 wind};
    NBC6.loadDirection = {dirForce6 'y'};
    NBC6.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
    NBC6.isFollower = [true; false];
    NBC6.isTimeDependent = [false; false];

    % Patch 7 :
    % _________

    FAmp7 = - loadAmplitude;
    NBC7.noCnd = 2;
    xib7 = [0 1];   etab7 = [0 1];   dirForce7 = 'normal';
    NBC7.xiLoadExtension = {xib7 xib7};
    NBC7.etaLoadExtension = {etab7 etab7};
    NBC7.loadAmplitude = {FAmp7 wind};
    NBC7.loadDirection = {dirForce7 'y'};
    NBC7.computeLoadVct = {'computeLoadVctAreaIGAThinStructure' 'computeLoadVctAreaIGAThinStructure'};
    NBC7.isFollower = [true; false];
    NBC7.isTimeDependent = [false; false];

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

    % Patch 9 :
    % _________

    FAmp9 = loadAmplitude;
    NBC9.noCnd = 1;
    xib9 = [0 1];   etab9 = [0 1];   dirForce9 = 'normal';
    NBC9.xiLoadExtension = {xib9};
    NBC9.etaLoadExtension = {etab9};
    NBC9.loadAmplitude = {FAmp9};
    NBC9.loadDirection = {dirForce9};
    NBC9.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC9.isFollower = true;
    NBC9.isTimeDependent = false;

    % Patch 10 :
    % __________

    FAmp10 = - loadAmplitude;
    NBC10.noCnd = 1;
    xib10 = [0 1];   etab10 = [0 1];   dirForce10 = 'normal';
    NBC10.xiLoadExtension = {xib10};
    NBC10.etaLoadExtension = {etab10};
    NBC10.loadAmplitude = {FAmp10};
    NBC10.loadDirection = {dirForce10};
    NBC10.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC10.isFollower = true;
    NBC10.isTimeDependent = false;

    % Patch 11 :
    % __________

    FAmp11 = loadAmplitude;
    NBC11.noCnd = 1;
    xib11 = [0 1];   etab11 = [0 1];   dirForce11 = 'normal';
    NBC11.xiLoadExtension = {xib11};
    NBC11.etaLoadExtension = {etab11};
    NBC11.loadAmplitude = {FAmp11};
    NBC11.loadDirection = {dirForce11};
    NBC11.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC11.isFollower = true;
    NBC11.isTimeDependent = false;

    % Patch 12 :
    % __________

    FAmp12 = - loadAmplitude;
    NBC12.noCnd = 1;
    xib12 = [0 1];   etab12 = [0 1];   dirForce12 = 'normal';
    NBC12.xiLoadExtension = {xib12};
    NBC12.etaLoadExtension = {etab12};
    NBC12.loadAmplitude = {FAmp12};
    NBC12.loadDirection = {dirForce12};
    NBC12.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC12.isFollower = true;
    NBC12.isTimeDependent = false;

    % Patch 13 :
    % __________

    FAmp13 = - loadAmplitude;
    NBC13.noCnd = 1;
    xib13 = [0 1];   etab13 = [0 1];   dirForce13 = 'normal';
    NBC13.xiLoadExtension = {xib13};
    NBC13.etaLoadExtension = {etab13};
    NBC13.loadAmplitude = {FAmp13};
    NBC13.loadDirection = {dirForce13};
    NBC13.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC13.isFollower = true;
    NBC13.isTimeDependent = false;

    % Patch 14 :
    % __________

    FAmp14 = loadAmplitude;
    NBC14.noCnd = 1;
    xib14 = [0 1];   etab14 = [0 1];   dirForce14 = 'normal';
    NBC14.xiLoadExtension = {xib14};
    NBC14.etaLoadExtension = {etab14};
    NBC14.loadAmplitude = {FAmp14};
    NBC14.loadDirection = {dirForce14};
    NBC14.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC14.isFollower = true;
    NBC14.isTimeDependent = false;

    % Patch 15 :
    % __________

    FAmp15 = - loadAmplitude;
    NBC15.noCnd = 1;
    xib15 = [0 1];   etab15 = [0 1];   dirForce15 = 'normal';
    NBC15.xiLoadExtension = {xib15};
    NBC15.etaLoadExtension = {etab15};
    NBC15.loadAmplitude = {FAmp15};
    NBC15.loadDirection = {dirForce15};
    NBC15.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC15.isFollower = true;
    NBC15.isTimeDependent = false;

    % Patch 16 :
    % _________

    FAmp16 = loadAmplitude;
    NBC16.noCnd = 1;
    xib16 = [0 1];   etab16 = [0 1];   dirForce16 = 'normal';
    NBC16.xiLoadExtension = {xib16};
    NBC16.etaLoadExtension = {etab16};
    NBC16.loadAmplitude = {FAmp16};
    NBC16.loadDirection = {dirForce16};
    NBC16.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC16.isFollower = true;
    NBC16.isTimeDependent = false;

    % Patch 17 :
    % __________

    FAmp17 = loadAmplitude;
    NBC17.noCnd = 1;
    xib17 = [0 1];   etab17 = [0 1];   dirForce17 = 'normal';
    NBC17.xiLoadExtension = {xib17};
    NBC17.etaLoadExtension = {etab17};
    NBC17.loadAmplitude = {FAmp17};
    NBC17.loadDirection = {dirForce17};
    NBC17.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC17.isFollower = true;
    NBC17.isTimeDependent = false;

    % Patch 18 :
    % __________

    FAmp18 = - loadAmplitude;
    NBC18.noCnd = 1;
    xib18 = [0 1];   etab18 = [0 1];   dirForce18 = 'normal';
    NBC18.xiLoadExtension = {xib10};
    NBC18.etaLoadExtension = {etab10};
    NBC18.loadAmplitude = {FAmp10};
    NBC18.loadDirection = {dirForce10};
    NBC18.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC18.isFollower = true;
    NBC18.isTimeDependent = false;

    % Patch 19 :
    % __________

    FAmp19 = loadAmplitude;
    NBC19.noCnd = 1;
    xib19 = [0 1];   etab19 = [0 1];   dirForce19 = 'normal';
    NBC19.xiLoadExtension = {xib19};
    NBC19.etaLoadExtension = {etab19};
    NBC19.loadAmplitude = {FAmp19};
    NBC19.loadDirection = {dirForce19};
    NBC19.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC19.isFollower = true;
    NBC19.isTimeDependent = false;

    % Patch 20 :
    % __________

    FAmp20 = - loadAmplitude;
    NBC20.noCnd = 1;
    xib20 = [0 1];   etab20 = [0 1];   dirForce20 = 'normal';
    NBC20.xiLoadExtension = {xib20};
    NBC20.etaLoadExtension = {etab20};
    NBC20.loadAmplitude = {FAmp20};
    NBC20.loadDirection = {dirForce20};
    NBC20.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC20.isFollower = true;
    NBC20.isTimeDependent = false;

    % Patch 21 :
    % __________

    FAmp21 = - loadAmplitude;
    NBC21.noCnd = 1;
    xib21 = [0 1];   etab21 = [0 1];   dirForce21 = 'normal';
    NBC21.xiLoadExtension = {xib21};
    NBC21.etaLoadExtension = {etab21};
    NBC21.loadAmplitude = {FAmp21};
    NBC21.loadDirection = {dirForce21};
    NBC21.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC21.isFollower = true;
    NBC21.isTimeDependent = false;

    % Patch 22 :
    % __________

    FAmp22 = loadAmplitude;
    NBC22.noCnd = 1;
    xib22 = [0 1];   etab22 = [0 1];   dirForce22 = 'normal';
    NBC22.xiLoadExtension = {xib22};
    NBC22.etaLoadExtension = {etab22};
    NBC22.loadAmplitude = {FAmp22};
    NBC22.loadDirection = {dirForce22};
    NBC22.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC22.isFollower = true;
    NBC22.isTimeDependent = false;

    % Patch 23 :
    % __________

    FAmp23 = - loadAmplitude;
    NBC23.noCnd = 1;
    xib23 = [0 1];   etab23 = [0 1];   dirForce23 = 'normal';
    NBC23.xiLoadExtension = {xib23};
    NBC23.etaLoadExtension = {etab23};
    NBC23.loadAmplitude = {FAmp23};
    NBC23.loadDirection = {dirForce23};
    NBC23.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC23.isFollower = true;
    NBC23.isTimeDependent = false;

    % Patch 24 :
    % _________

    FAmp24 = loadAmplitude;
    NBC24.noCnd = 1;
    xib24 = [0 1];   etab24 = [0 1];   dirForce24 = 'normal';
    NBC24.xiLoadExtension = {xib24};
    NBC24.etaLoadExtension = {etab24};
    NBC24.loadAmplitude = {FAmp24};
    NBC24.loadDirection = {dirForce24};
    NBC24.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
    NBC24.isFollower = true;
    NBC24.isTimeDependent = false;

    % Collect all the Neumann boundary conditions into an arra<y
    NBC = {NBC1 NBC2 NBC3 NBC4 NBC5 NBC6 NBC7 NBC8 ...
        NBC9 NBC10 NBC11 NBC12 NBC13 NBC14 NBC15 NBC16 ...
        NBC17 NBC18 NBC19 NBC20 NBC21 NBC22 NBC23 NBC24};
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
    noConnections = 52; %12 % 24 % 44 % 52

    % Define connections by patch numbers
    connections.No = noConnections;
    connections.xiEtaCoup = zeros(noConnections,10);
    connections.xiEtaCoup(:,:) = [2 1 xicoup2_1 etacoup2_1 xicoup1_2 etacoup1_2 % -> First tube
                                  2 3 xicoup2_3 etacoup2_3 xicoup3_2 etacoup3_2
                                  4 3 xicoup4_3 etacoup4_3 xicoup3_4 etacoup3_4
                                  4 1 xicoup4_1 etacoup4_1 xicoup1_4 etacoup1_4
                                  1 5 xicoup1_5 etacoup1_5 xicoup5_1 etacoup5_1
                                  6 2 xicoup6_2 etacoup6_2 xicoup2_6 etacoup2_6
                                  3 7 xicoup3_7 etacoup3_7 xicoup7_3 etacoup7_3
                                  8 4 xicoup8_4 etacoup8_4 xicoup4_8 etacoup4_8
                                  5 6 xicoup5_6 etacoup5_6 xicoup6_5 etacoup6_5
                                  7 6 xicoup7_6 etacoup7_6 xicoup6_7 etacoup6_7
                                  7 8 xicoup7_8 etacoup7_8 xicoup8_7 etacoup8_7
                                  5 8 xicoup5_8 etacoup5_8 xicoup8_5 etacoup8_5 % <- First tube
                                  9 10 xicoup9_10 etacoup9_10 xicoup10_9 etacoup10_9 % <- Second tube
                                  10 11 xicoup10_11 etacoup10_11 xicoup11_10 etacoup11_10
                                  11 12 xicoup11_12 etacoup11_12 xicoup12_11 etacoup12_11
                                  9 12 xicoup9_12 etacoup9_12 xicoup12_9 etacoup12_9
                                  9 13 xicoup9_13 etacoup9_13 xicoup13_9 etacoup13_9
                                  10 14 xicoup10_14 etacoup10_14 xicoup14_10 etacoup14_10
                                  11 15 xicoup11_15 etacoup11_15 xicoup15_11 etacoup15_11
                                  12 16 xicoup12_16 etacoup12_16 xicoup16_12 etacoup16_12
                                  13 14 xicoup13_14 etacoup13_14 xicoup14_13 etacoup14_13
                                  14 15 xicoup14_15 etacoup14_15 xicoup15_14 etacoup15_14
                                  15 16 xicoup15_16 etacoup15_16 xicoup16_15 etacoup16_15
                                  13 16 xicoup13_16 etacoup13_16 xicoup16_13 etacoup16_13 % <- Second tube
                                  17 18 xicoup17_18 etacoup17_18 xicoup18_17 etacoup18_17 % -> Third Tube
                                  18 19 xicoup18_19 etacoup18_19 xicoup19_18 etacoup19_18
                                  19 20 xicoup19_20 etacoup19_20 xicoup20_19 etacoup20_19
                                  17 20 xicoup17_20 etacoup17_20 xicoup20_17 etacoup20_17
                                  17 21 xicoup17_21 etacoup17_21 xicoup21_17 etacoup21_17
                                  18 22 xicoup18_22 etacoup18_22 xicoup22_18 etacoup22_18
                                  19 23 xicoup19_23 etacoup19_23 xicoup23_19 etacoup23_19
                                  20 24 xicoup20_24 etacoup20_24 xicoup24_20 etacoup24_20
                                  21 22 xicoup21_22 etacoup21_22 xicoup22_21 etacoup22_21
                                  22 23 xicoup22_23 etacoup22_23 xicoup23_22 etacoup23_22
                                  23 24 xicoup23_24 etacoup23_24 xicoup24_23 etacoup24_23
                                  21 24 xicoup21_24 etacoup21_24 xicoup24_21 etacoup24_21 % <- Third Tube
                                  5 14 xicoup5_14 etacoup5_14 xicoup14_5 etacoup14_5 % -> Coupling betwenn the first and the second tube
                                  5 15 xicoup5_15 etacoup5_15 xicoup15_5 etacoup15_5
                                  8 14 xicoup8_14 etacoup8_14 xicoup14_8 etacoup14_8
                                  8 15 xicoup8_15 etacoup8_15 xicoup15_8 etacoup15_8
                                  1 10 xicoup1_10 etacoup1_10 xicoup10_1 etacoup10_1
                                  1 11 xicoup1_11 etacoup1_11 xicoup11_1 etacoup11_1
                                  4 10 xicoup4_10 etacoup4_10 xicoup10_4 etacoup10_4
                                  4 11 xicoup4_11 etacoup4_11 xicoup11_4 etacoup11_4 % <- Coupling betwenn the first and the second tube
                                  13 22 xicoup13_22 etacoup13_22 xicoup22_13 etacoup22_13 % -> Coupling betwenn the second and the third tube
                                  13 23 xicoup13_23 etacoup13_23 xicoup23_13 etacoup23_13
                                  16 22 xicoup16_22 etacoup16_22 xicoup22_16 etacoup22_16
                                  16 23 xicoup16_23 etacoup16_23 xicoup23_16 etacoup23_16
                                  9 18 xicoup9_18 etacoup9_18 xicoup18_9 etacoup18_9
                                  9 19 xicoup9_19 etacoup9_19 xicoup19_9 etacoup19_9
                                  12 18 xicoup12_18 etacoup12_18 xicoup18_12 etacoup18_12
                                  12 19 xicoup12_19 etacoup12_19 xicoup19_12 etacoup19_12]; % <- Coupling betwenn the second and the third tube
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
% msh.nodes = [];
% msh.elements = [];
% labelsEnabled = true;
% graph.index = plot_referenceConfigIGAMortarMapping...
%     (BSplinePatches,msh,labelsEnabled,graph);

%% Compute the load vector for the visualization of the reference configuration
% for counterPatches = 1:noPatches
%     BSplinePatches{counterPatches}.FGamma = ...
%         zeros(3*BSplinePatches{counterPatches}.noCPs,1);
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
% end

%% Plot the reference configuration for the multipatch geometry
% graph.index = graph.index + 1;
% graph.isPrestressEnabled = true;
% graph.compPrestress = '2';
% color = [217 218 219]/255;
% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatches,connections,color,graph,'outputEnabled');
% [az,el] = view;

%% Create Lagrange Multipiers fields for all interfaces
fprintf('**************************************************\n');
fprintf('* Creating interface Lagrange Multipliers fields *\n');
fprintf('**************************************************\n\n');
for iConnections = 1:connections.No
    %% Check method for weak enforcement of multipatch coupling
    if ~strcmp(method,'LagrangeMultipliers') && ~strcmp(method,'AugmentedLagrangeMultipliers')
        break;
    end
    
    %% Get the IDs of the patches involved
    idI = connections.xiEtaCoup(iConnections,1);
    idJ = connections.xiEtaCoup(iConnections,2);
    fprintf('Coupling between patches %d and %d \n',idI,idJ);
    fprintf('---------------------------------- \n\n');
    
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
        fprintf('The coupling interface of patch %d is along xi\n',idI);
    elseif xiCoupI(1,1) == xiCoupI(1,2) && etaCoupI(1,1) ~= etaCoupI(1,2)
        isOnXiI = false;
        fprintf('The coupling interface of patch %d is along eta\n',idI);
    else
        error('Either the xi or the eta direction has to be fixed along the interface for patch %d',idI);
    end
    xiCoupJ = connections.xiEtaCoup(iConnections,7:8);
    etaCoupJ = connections.xiEtaCoup(iConnections,9:10);
    if xiCoupJ(1,1) ~= xiCoupJ(1,2) && etaCoupJ(1,1) == etaCoupJ(1,2)
        isOnXiJ = true;
        fprintf('The coupling interface of patch %d is along xi\n',idJ);
    elseif xiCoupJ(1,1) == xiCoupJ(1,2) && etaCoupJ(1,1) ~= etaCoupJ(1,2)
        isOnXiJ = false;
        fprintf('The coupling interface of patch %d is along eta\n',idJ);
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
    fprintf('Degree elevating the Lagrange Multipliers field to %d\n',pLM);
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
        if strcmp(method,'LagrangeMultipliers')
            tpLambda = pLM - pLambda + 1;
        elseif strcmp(method,'AugmentedLagrangeMultipliers')
            tpLambda = pLM - pLambda;
        end
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
    scaleLM = .5; % .5
    noLambda = ceil(min([noKnotsI noKnotsJ])*scaleLM);    
    fprintf('Uniformly inserting %d knots in the Lagrange Multipliers field\n',noLambda);
    [XiLambdaLM,CPLambdaLM] = knotRefineUniformlyBSplineCurve...
        (noLambda,pLambda,XiLambda,CPLambda,'');
    scaleALM = .2;
    noLambda = ceil(min([noKnotsI noKnotsJ])*scaleALM);
    fprintf('Uniformly inserting %d knots in the augmented Lagrange Multipliers field\n',noLambda);
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
        pLM = polOrderPatch + 1; % polOrderPatch - 2

        if pLM <= 0
            pLambda = 0;
            XiLambda = [0 1];
            CPLambda = zeros(1,4);
        else
            pLambda = 1;
            XiLambda = [0 0 1 1];
            CPLambda = zeros(2,4);

            % Peform a p-refinement
            tpLambda = pLM - pLambda; % polOrderPatch - pLM;
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
        percentage = 1.4; % 1.0
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
            tpLambda = pLM - pLambda; % polOrderPatch - pLM;
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
        percentage = .25;
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

%% Nonlinear analysis properties

% number of load steps for the non-linear analysis
propNLinearAnalysis.method = 'newtonRapshon';

% number of load steps for the non-linear analysis
propNLinearAnalysis.noLoadSteps = 1;

% Assign a tolerance for the Newton iterations
propNLinearAnalysis.eps = 1e-4;

% Assign the maximum number of iterations
propNLinearAnalysis.maxIter = 35;

%% Perform modal analysis and visualize the chosen eigenmode shape

% Solve the eigenvalue problem
noEig = 100; % 100, 200, 400
[eigenmodeShapes,naturalFrequencies,BSplinePatches,dHat] = ...
    solve_DDMEigenmodeAnalysisIGAMembrane...
    (BSplinePatches,connections,propCoupling,solve_LinearSystem,noEig,...
    propNLinearAnalysis,'outputEnabled');
save(['./data_ModalAnalysisPool/' 'data_ModalAnalysis' '_' meshSize '_' method]);
return;

%% Visualize selected mode shapes

% Find the numbering of the DOFs which are purely associated with
% displacements (that is no Lagrange Multiplier DOFs)
noDOFsWTLM = 0;
for iPatches = 1:noPatches
    noDOFsWTLMPatch = 3*BSplinePatches{iPatches}.noCPs;
    noDOFsWTLM = noDOFsWTLM + noDOFsWTLMPatch;
end
EFTPatchesWTLM = zeros(1,noDOFsWTLM);
noDOFsSaved = 0;
for iPatches = 1:noPatches
    noDOFsWTLMPatch = 3*BSplinePatches{iPatches}.noCPs;
    EFTPatchesWTLM(1,noDOFsSaved + 1:noDOFsSaved + noDOFsWTLMPatch) = ...
        BSplinePatches{iPatches}.EFTPatches(1,1:noDOFsWTLMPatch);
    noDOFsSaved = noDOFsSaved + noDOFsWTLMPatch;
end

% Create a new B-Spline array with displaced Control Points according to
% the static equilibrium configuration computed at the static step
BSplinePatchesEq = BSplinePatches;
for iPatches = 1:noPatches
    BSplinePatchesEq{iPatches}.CP = computeDisplacedControlPointsForIGAKirchhoffLoveShell...
    (BSplinePatches{iPatches}.CP,dHat(EFTPatchesWTLM));
end

% Eigenfrequency 
idEig = 1; % 1, 2, 4

% Visualization of the eigenmode shapes
for i = idEig
    scalingFct = -40;
    noEigen = i;
    graph.index = plot_postprocIGAMembraneMultipatchesNLinear...
        (BSplinePatchesEq,scalingFct*eigenmodeShapes(EFTPatchesWTLM,noEigen),graph,'');
    az = -37.500000000000000;
    el = 30;
    view(az,el);
%     camlightHandle = camlight('left');
%     camlightHandle = camlight(-290,-320); % (0,40)
%     camlightHandle = camlight(0,40);
%     camlightHandle = camlight(-37.500000000000007,69.996738700909546);
    camlightHandle = camlight(-0.312793885148904e2 + 180,-0.370537693954654e2 + 180); % Mode shape 1
    lightHandle = findobj(gcf,'Type','Light')
    lightHandle.Position
    camHandle = camlight;
    lighting phong;
    xLimits = xlim;
    yLimits = ylim;
    zLimits = zlim;
    axis off;
    title('');
end

% Write the geometry for GiD
BSplinePatchesNitsche = computeUpdatedGeometryIGAThinStructureMultipatches...
    (BSplinePatches,scalingFct*eigenmodeShapes(EFTPatchesWTLM,idEig));
for iPatches = 1:noPatches
    BSplinePatchesNitsche{iPatches}.CP = BSplinePatchesNitsche{iPatches}.CPd;
end

% graph.index = plot_referenceConfigurationIGAThinStructureMultipatches...
%     (BSplinePatchesNitsche,connections,[217 218 219]/255,graph,'');

pathToOutput = '../../../outputGiD/isogeometricMembraneAnalysis/';
writeOutMultipatchBSplineSurface4GiD(BSplinePatchesNitsche,pathToOutput,[caseName '_' meshSize '_' method]);

%% Load a reference solution for the modal analysis
mesh.nodes = importdata('../../../preComputedData/isogeometricMembraneAnalysis/inflatedTubes/referenceFEMSolution_ModalAnalysis/nodes');
elements = importdata('../../../preComputedData/isogeometricMembraneAnalysis/inflatedTubes/referenceFEMSolution_ModalAnalysis/elements');
mesh.elements = elements(:,1:4);
fileNameFEM = '../../../preComputedData/isogeometricMembraneAnalysis/inflatedTubes/referenceFEMSolution_ModalAnalysis/displacementsModeShapes';
naturalFrequenciesFEM = importdata('../../../preComputedData/isogeometricMembraneAnalysis/inflatedTubes/referenceFEMSolution_ModalAnalysis/eigenfrequencies100');
naturalFrequenciesFEM = naturalFrequenciesFEM.data(:,2);
fstring = fileread(fileNameFEM);
stringSeparator = ['Result "Eigenvector" "Eigenvector nmb." ' num2str(idEig) ' Vector OnNodes'];
block = regexp(fstring,stringSeparator,'split');
block = regexp(block{2},'Values','split');
block(1) = [];
block(2:end) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
displacements = cell2mat(out);

%% Make the numbering of the nodes sequential
for iNodes = 1:length(mesh.nodes(:,1))
    if mesh.nodes(iNodes,1) ~= iNodes
        [idI,idJ] = find(mesh.nodes(iNodes,1) == mesh.elements(:,:));
        for k = 1:length(idI)
            mesh.elements(idI(k),idJ(k)) = iNodes;
        end
    end
end

%% Visualize the chosen reference modal shape
scaling = -40.0;
displacementsOrdered = zeros(3*length(displacements(:,1)),1);
for iDisp = 1:length(displacements(:,1))
%     [index,~] = find(displacements(iNodes,1) == mesh.nodes(:,1));
    displacementsOrdered(3*iDisp - 2,1) = scaling*displacements(iDisp,2);
    displacementsOrdered(3*iDisp - 1,1) = scaling*displacements(iDisp,3);
    displacementsOrdered(3*iDisp,1) = scaling*displacements(iDisp,4);
end
resultant = 'displacement';
component = '2norm';
graph.visualization.geometry = 'current'; % 'reference', 'current', 'reference_and_current'
mesh.nodes = mesh.nodes(:,2:4);
graph.index = plot_currentConfigurationAndResultants...
    (mesh,[],displacementsOrdered,[],[],resultant,component,graph);
title('');
axis on;
grid on;
epsY = 5e-2;
epsZ = 5e-2;
xlim(xLimits);
ylim([yLimits(1,1) - epsY yLimits(1,2)]);
zlim([zLimits(1,1) zLimits(1,2) + epsZ]);
az = -37.500000000000000;
el = 30;
view(az,el);
camlightHandle = camlight(-0.312793885148904e2 + 180,-0.370537693954654e2 + 180);
% axis on;
% colormap('jet');
% h = colorbar;
% limitsColorbar = get(h,'Limits');

lgt = light('Position',lightHandle(1).Position,'Style','local');
[azLight,elLight] = lightangle(lgt)

% camlight(0,230);
% camlight(180,230);

%% END
