function [stiffMtx, F, BSplinePatch, propCoupling, minElAreaSize] = ...
    computeStiffMtxAndLoadVctIGAPlateInMembraneActionLinearPagewise...
    (stiffMtx, tanMtxLoad, dHat, dHatSaved, dHatDot, dHatDotSaved, ...
    BSplinePatch, connections, propCoupling, loadFactor, noPatch, ...
    noTimeStep, noNonlinearIteration, noWeakDBCCnd, t, ...
    propTransientAnalysis, isReferenceUpdated, tab, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the linear stiffness matrix and force vector for the plate in
% membrane action analysis.
%
%                 Input :
%              constMtx : Constant precomputed part of the stiffness matrix
%            tanMtxLoad : Tangent stiffness matrix resulting from the
%                         application of followe loads
%                  dHat : The displacement field of the previous iteration
%                         step (dummy variable for this function)
%             dHatSaved : The displacement field of the previous time step
%                         (dummy variable for this function)
%               dHatDot : The velocity field of the previous iteration step
%                         (dummy variable for this function)
%          dHatDotSaved : The velocity field of the previous time step 
%                        (dummy variable for this function)
%          BSplinePatch : B-Spline patch with the following parameters:
%                                 .p,.q: Polynomial degrees
%                              .Xi,.Eta: knot vectors
%                                   .CP: Control Points coordinates and 
%                                        weights
%                              .isNURBS: Flag on whether the basis is a 
%                                        NURBS or a B-Spline
%                             .homDOFs : The global numbering of the
%                                        DOFs where homogeneous Dirichlet
%                                        boundary conditions are applied
%                           .inhomDOFs : The global numbering of the
%                                        DOFs where homogeneous Dirichlet
%                                        boundary conditions are applied
%                     .valuesInhomDOFs : Prescribed values to the DOFs 
%                                        where homogeneous Dirichlet
%                                        boundary conditions are applied
%                               FGamma : The boundary applied force vector
%                                        over the B-Spline patch
%                           bodyForces : Function handle to the computation 
%                                        of the body forces
%                        .DOFNumbering : Numbering of the DOFs sorted into
%                                        a 3D array
%                          .parameters : material parameters of the 
%                                        membrane
%                                 .int : On the numerical integration
%                                         .type : 'default' or 'user'
%                                        .xiNGP : No. of GPs along xi-
%                                                 direction for stiffness 
%                                                 entries
%                                       .etaNGP : No. of GPs along eta-
%                                                 direction for stiffness 
%                                                 entries
%                                 .xiNGPForLoad : No. of GPs along xi-
%                                                 direction for load 
%                                                 entries
%                                .etaNGPForLoad : No. of GPs along eta-
%                                                 direction for load 
%                                                 entries
%                                   .nGPForLoad : No. of GPs along boundary
%           connections : Dummy variable for this function
%          propCoupling : Dummy variable for this function
%            loadFactor : The load factor of the current time step
%            noTimeStep : Number of time step
%  noNonlinearIteration : Number of nonlinear iteration step
%          noWeakDBCCnd : Number of weak Dirichlet boundary conditions
%                     t : The time instance
% propTransientAnalysis : Structure on the transient analysis :
%                            .timeDependence : 'true' or 'false'
%                               .noTimeSteps : Number of time steps
%    isReferenceUpdated : Flag on whether the reference configuration is
%                         updated (mainly used in form-finding analysis)
%                   tab : Tabulation related to outputting information on
%                         the command window
%                outMsg : Enables outputting information onto the command
%                         window when chosen as 'outputEnabled'
%
%                Output :
%              stiffMtx : The master stiffness matrix
%                     F : The global force vector
%          BSplinePatch : The updated with the stabilization parameters
%                         BSpline patch array when weak enforcement of the
%                         Dirichlet boundary conditions is chosen
%          propCoupling : The updated with the stabilization terms coupling 
%                         properties array
%         minElAreaSize : Dummy output for this function
%             
% Function layout :
%
% 1. Loop over the elements in the B-Spline patch and reshape necessary 3D arrays if the stiffness matrix is not precomputed
%
% 2. Loop over all Gauss points if the stiffness matrix is not precomputed
% ->
%    2i. Compute the basis functions matrices pagewise
%
%   2ii. Get the covariant basis of the reference configuration pagewise
%
%  2iii. Compute the surface normal base vector pagewise
%
%   2iv. Compute the covariant basis of the reference configuration pagewise
%
%    2v. Compute element tangent stiffness matrix
%
%   2vi. Add the contributions from the Gauss Point
% <-
%
% 3. Assemble to the global tangent matrix
%
% 4. Compute the exernally applied load vector
%
% 5. Check output
%
%% Function main body

%% 0. Read input

% Check the NURBS geometry input
isBSplinePatchCell = false;
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
    isBSplinePatchCell = true;
end

% Get the patch parameters
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
parameters = BSplinePatch.parameters;

% Number of knots in xi-,eta-direction
numKnots_xi = length(Xi);
numKnots_eta = length(Eta);
numCPs_xi = length(CP(:, 1, 1));
numCPs_eta = length(CP(1, :, 1));

% Check input
checkInputForBSplineSurface ...
    (p, numKnots_xi, numCPs_xi, q, numKnots_eta, numCPs_eta);

% Number of DOFs per Control Point
numDOFsPerCP = 2;

% Number of degrees of freedom for the whole structure
numDOFs = length(dHat);

% Check the precomputed part of the stiffness matrix
isPrecomputed = false;
if isnumeric(stiffMtx)
    numDOFs_precomputed =  length(stiffMtx);
    if numDOFs_precomputed ~= numDOFs
        error('The precomputed part of the stiffness matrix has %d DOFs but the system has %d DOFs', ...
            numDOFs_precomputed, numDOFs);
    else
        isPrecomputed = true;
    end
end
if ~isfield(BSplinePatch, 'isStiffMtxPrecom')
    error('Structure BSplinePatch must define field isStiffMtxPrecom')
else
    if ~islogical(BSplinePatch.isStiffMtxPrecom)
        error('BSplinePatch.isStiffMtxPrecom must be a boolean');
    else
        isStiffMtxPrecom = BSplinePatch.isStiffMtxPrecom;
    end
end
if isStiffMtxPrecom && ~isPrecomputed
    error('Flag BSplinePatch.isStiffMtxPrecom regarding the existence of a precomputed stiffness matrix is true but provided constant part of the stiffness matrix is not provided');
end

% Number of Gauss Points in the element level
numGPsEl = BSplinePatch.noGPsEl;

% Number of Control Points and DOFs in the element level
numCPsEl = (p + 1)*(q + 1);
numDOFsEl = numDOFsPerCP*numCPsEl;

% Initialize the pagewise material matrix
DmPage = zeros(BSplinePatch.noElmnts, 3, 3);

% Initialize the pagewise displaced Control Point array
CPVctPage = zeros(BSplinePatch.noElmnts, numDOFsEl);

% Minimum element edge size is returned as a dummy array
minElAreaSize = 'undefined';

% Initialize pagewise basis function matrices
dRdxiMatrixPage = zeros(BSplinePatch.noElmnts, 3, numDOFsEl);
dRdetaMatrixPage = zeros(BSplinePatch.noElmnts, 3, numDOFsEl);

% Initialize pagewise element stiffness matrices and residual vectors
if ~isStiffMtxPrecom
    stiffMtxElPage = zeros(BSplinePatch.noElmnts, numDOFsEl, numDOFsEl);
end

%% 1. Loop over the elements in the B-Spline patch and reshape necessary 3D arrays if the stiffness matrix is not precomputed
if ~isStiffMtxPrecom
    for iElmnts = 1:BSplinePatch.noElmnts
        %% 1i. Assign the material matrices into a pagewise array
        DmPage(iElmnts, :, :) = parameters.E*parameters.t/(1 - parameters.nue^2)*...
            [1              parameters.nue 0
             parameters.nue 1              0
             0              0              (1 - parameters.nue)/2];

        %% 1iii. Assign the 3D array for the displaced Control Point array
        CPPage = BSplinePatch.CP(BSplinePatch.xiIndexCP(iElmnts, :), ...
            BSplinePatch.etaIndexCP(iElmnts, :), 1:numDOFsPerCP);
        CPVctPage(iElmnts, :) = reshape(reshape(CPPage, [(p + 1)*(q + 1), numDOFsPerCP])', ...
            [numDOFsPerCP*(p + 1)*(q + 1), 1]);
    end
end

%% 2. Loop over all Gauss points if the stiffness matrix is not precomputed
if ~isStiffMtxPrecom
    for iGPs = 1:numGPsEl
        %% 2i. Compute the basis functions matrices pagewise
        %
        % for iCPs = 1:noCPsEl
        %     % dR/dxi
        %     dRdxiMatrix(1, 2*iCPs - 1) = dR(iCPs, 2);
        %     dRdxiMatrix(2, 2*iCPs) = dR(iCPs, 2);
        %     
        %     % dR/deta
        %     dRdetaMatrix(1, 2*iCPs - 1) = dR(iCPs, 3);
        %     dRdetaMatrix(2, 2*iCPs) = dR(iCPs, 3);
        % end
        %
        % dR/dxi
        dRdxiMatrixPage(:, 1, numDOFsPerCP*(1:numCPsEl) - numDOFsPerCP + 1) = ...
            BSplinePatch.dRdXi(:, iGPs, :);
        dRdxiMatrixPage(:, 2, numDOFsPerCP*(1:numCPsEl) - numDOFsPerCP + 2) = ...
            BSplinePatch.dRdXi(:, iGPs, :);

        % dR/deta
        dRdetaMatrixPage(:, 1, numDOFsPerCP*(1:numCPsEl) - numDOFsPerCP + 1) = ...
            BSplinePatch.dRdEta(:, iGPs, :);
        dRdetaMatrixPage(:, 2, numDOFsPerCP*(1:numCPsEl) - numDOFsPerCP + 2) = ...
            BSplinePatch.dRdEta(:, iGPs, :);

        %% 2ii. Get the covariant basis of the reference configuration pagewise
        G1GPPage = BSplinePatch.GXi(:, iGPs, :);
        G2GPPage = BSplinePatch.GEta(:, iGPs, :);
        GGPCovariantPage = cat(2, G1GPPage, G2GPPage);

        %% 2iii. Compute the surface normal base vector pagewise
        elementAreaOnGP = BSplinePatch.elementAreaOnGP(:, iGPs);

        %% 2iv. Compute the covariant basis of the reference configuration pagewise
        G1Page = pmtimes(dRdxiMatrixPage, CPVctPage);
        G2Page = pmtimes(dRdetaMatrixPage, CPVctPage);

        %% 2v. Compute element tangent stiffness matrix
        stiffMtxElOnGPPage = ...
            computeElStiffMtxAndLoadVctPlateInMembraneActionLinearPagewise ...
            (dRdxiMatrixPage, dRdetaMatrixPage, GGPCovariantPage, G1Page, ...
            G2Page, DmPage);

        %% 2vi. Add the contributions from the Gauss Point
        stiffMtxElPage = stiffMtxElPage + ...
            pstimes(stiffMtxElOnGPPage, elementAreaOnGP);
    end
end

%% 3. Assemble to the global tangent matrix
if ~isStiffMtxPrecom
    if isPrecomputed
        stiffMtx = stiffMtx + assembleSparseMatricies ...
            (BSplinePatch.EFT, numDOFs, numDOFsEl, stiffMtxElPage);
    else
        stiffMtx = assembleSparseMatricies ...
            (BSplinePatch.EFT, numDOFs, numDOFsEl, stiffMtxElPage);
    end
end

%% 4. Compute the exernally applied load vector
if isfield(BSplinePatch, 'FGamma')
    F = BSplinePatch.FGamma;
else
    F = 'undefined';
end

%% 5. Check output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end
