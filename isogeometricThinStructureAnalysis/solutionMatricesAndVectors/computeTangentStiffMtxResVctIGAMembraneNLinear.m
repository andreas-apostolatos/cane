function [tanStiffMtx,resVct,BSplinePatch,propCoupling,minElAreaSize] = ...
    computeTangentStiffMtxResVctIGAMembraneNLinear...
    (constMtx,tanMtxLoad,dHat,dHatSaved,dHatDot,dHatDotSaved,BSplinePatch,...
    connections,propCoupling,loadFactor,noPatch,noTimeStep,...
    noNonlinearIteration,noWeakDBCCnd,t,propTransientAnalysis,...
    isReferenceUpdated,tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangent stiffness matrix, mass matrix and internal force 
% vector for the nonlinear analysis case to the isogeometrc membrane and if 
% transient for the appropriate time integration scheme.
%
%                 Input :
%              constMtx : The constant part of the tangent stiffness matrix
%                         from the application of weak boundary conditions
%                         using the penalty method
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
%           tanStiffMtx : The master tangent stiffness matrix
%                resVct : The global residual vector
%          BSplinePatch : The updated with the stabilization parameters
%                         BSpline patch array when weak enforcement of the
%                         Dirichlet boundary conditions is chosen
%          propCoupling : The updated with the stabilization terms coupling 
%                         properties array
%         minElAreaSize : Dummy output for this function
%             
% Function layout :
%
% 0. Read input
%
% 1. Loop over the elements in the B-Spline patch and reshape necessary 3D arrays
% ->
%    1i. Assign the material matrices into a pagewise array
%
%   1ii. Assign the 3D array for the reference Control Point array
%
%  1iii. Assign the 3D array for the displaced Control Point array
% <-
%
% 2. Get the function handle for the computation of the tangent stiffness and the residual load vector corresponding to the application of weak Dirichlet boundary conditions
%
% 3. Loop over all Gauss points
% ->
%    3i. Compute the basis functions matrices pagewise
%
%   3ii. Compute the prestress Voigt vectors pagewise
%
%  3iii. Get the covariant basis of the reference configuration pagewise
%
%   3iv. Compute the surface normal base vector pagewise
%
%    3v. Compute the covariant basis of the current configuration pagewise
%
%   3vi. Compute element tangent stiffness matrix, the residual internal load vector and the stiffness matrix needed for estimation of the stabilization parameter corresponding to the application of weak boundary conditions with the Nitsche method
%
%  3vii. Add the contributions from the Gauss Point
% <-
%
% 4. Assemble to the global tangent matrix and the stiffness matrix needed for estimation of the stabilization parameter corresponding to the application of weak boundary conditions with the Nitsche method
%
% 5. Add existing tangent matrix from the application of follower loads
%
% 6. Assemble to the global residual vector
%
% 7. Update the residual load vector to account for the external loading
%
% 8. Updated the tangent stiffness matrix and the residual vector with the contributions resulting from the constant matrix
%
% 9. Compute the tangent stiffness matrix and residual load vector resulting from embedded into the patch cables
%
% 10. Compute the tangent matrix and residual vector contributions due to the application of weak boundary conditions
%
% 11. Check output
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

% Check the tangent matrix from the application of a follower load
if iscell(tanMtxLoad)
    if length(tanMtxLoad) > 1
        error('Multipatch tangent matrix due to follower loading is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        tanMtxLoad = tanMtxLoad{1};
    end
end

% Get the patch parameters
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
parameters = BSplinePatch.parameters;
thickness = parameters.t;

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Number of DOFs per Control Point
noDOFsPerCP = 3;

% Number of degrees of freedom for the whole structure
noDOFs = length(dHat);

% Number of Gauss Points in the element level
noGPsEl = BSplinePatch.noGPsEl;

% Number of Control Points and DOFs in the element level
noCPsEl = (p+1)*(q+1);
noDOFsEl = noDOFsPerCP*noCPsEl;

% Initialize the pagewise material matrix
DmPage = zeros(BSplinePatch.noElmnts,3,3);

% Initialize the pagewise reference Control Point array
if isReferenceUpdated
    CPVctPage = zeros(BSplinePatch.noElmnts,noDOFsEl);
end

% Initialize the pagewise displaced Control Point array
CPdVctPage = zeros(BSplinePatch.noElmnts,noDOFsEl);

% On the application of weak Dirichlet boundary conditions
isWeakDBC = false;
if ~isempty(BSplinePatch.weakDBC)
    if BSplinePatch.weakDBC.noCnd > 0
        isWeakDBC = true;
    end
end

% Minimum element edge size is returned as a dummy array
minElAreaSize = 'undefined';

% Initialize pagewise basis function matrices
dRdxiMatrixPage = zeros(BSplinePatch.noElmnts,3,noDOFsEl);
dRdetaMatrixPage = zeros(BSplinePatch.noElmnts,3,noDOFsEl);

% Initialize pagewise element stiffness matrices and residual vectors
tangStiffMtxElPage = zeros(BSplinePatch.noElmnts,noDOFsEl,noDOFsEl);
if ~ischar(loadFactor)
    resVctElPage = zeros(BSplinePatch.noElmnts,noDOFsEl);
end

%% 1. Loop over the elements in the B-Spline patch and reshape necessary 3D arrays
for iElmnts = 1:BSplinePatch.noElmnts
    %% 1i. Assign the material matrices into a pagewise array
    DmPage(iElmnts,:,:) = parameters.E*parameters.t/(1-parameters.nue^2)*...
        [1              parameters.nue 0
         parameters.nue 1              0
         0              0              (1-parameters.nue)/2];
     
    %% 1ii. Assign the 3D array for the reference Control Point array
    if isReferenceUpdated
        CPPage = BSplinePatch.CP(BSplinePatch.xiIndexCP(iElmnts,:),BSplinePatch.etaIndexCP(iElmnts,:),1:3);
        CPVctPage(iElmnts,:) = reshape(reshape(CPPage,[(p+1)*(q+1),3])',[3*(p+1)*(q+1),1]);
    end
    
    %% 1iii. Assign the 3D array for the displaced Control Point array
    CPdPage = BSplinePatch.CPd(BSplinePatch.xiIndexCP(iElmnts,:),BSplinePatch.etaIndexCP(iElmnts,:),1:3);
    CPdVctPage(iElmnts,:) = reshape(reshape(CPdPage,[(p+1)*(q+1),3])',[3*(p+1)*(q+1),1]);
end

%% 2. Get the function handle for the computation of the tangent stiffness and the residual load vector corresponding to the application of weak Dirichlet boundary conditions
if isWeakDBC
    if isfield(BSplinePatch.weakDBC,'method')
        if isfield(BSplinePatch.weakDBC,'imposedMotion')    
            if strcmp(BSplinePatch.weakDBC.method,'penalty')
                computeTangMtxResVct = @computeWeakDBCTangMtxResVctPenaltyIGAMembrane;
            elseif strcmp(BSplinePatch.weakDBC.method,'lagrangeMultipliers')
                computeTangMtxResVct = @computeWeakDBCTangMtxResVctLagrangeMultipliersIGAMembrane;
            end
        else
            if strcmp(BSplinePatch.weakDBC.method,'penalty') || ...
                    strcmp(BSplinePatch.weakDBC.method,'lagrangeMultipliers')
                computeTangMtxResVct = 'undefined';
            end
        end
        if strcmp(BSplinePatch.weakDBC.method,'nitsche')
            computeTangMtxResVct = @computeWeakDBCTangMtxResVctNitscheIGAMembrane;
        elseif ~strcmp(BSplinePatch.weakDBC.method,'penalty') && ...
               ~strcmp(BSplinePatch.weakDBC.method,'lagrangeMultipliers') && ...
               ~strcmp(BSplinePatch.weakDBC.method,'nitsche')
            error('Define a valid method for the application of weak Dirichlet boundary conditions in BSplinePatch.weakDBC.method')
        end
    else
        error('Define method for the application of weak Dirichlet boundary conditions in BSplinePatch.weakDBC');
    end
end
isWeakDBCTangent = false;
if isWeakDBC && isa(computeTangMtxResVct,'function_handle')
    isWeakDBCTangent = true;
end

%% 3. Loop over all Gauss points
for iGPs = 1:noGPsEl
    %% 3i. Compute the basis functions matrices pagewise
    %
    % for iCPs = 1:noCPsEl
    %     % dR/dxi
    %     dRdxiMatrix(1,3*iCPs-2) = dR(iCPs,2);
    %     dRdxiMatrix(2,3*iCPs-1) = dR(iCPs,2);
    %     dRdxiMatrix(3,3*iCPs) = dR(iCPs,2);
    %     
    %     % dR/deta
    %     dRdetaMatrix(1,3*iCPs-2) = dR(iCPs,3);
    %     dRdetaMatrix(2,3*iCPs-1) = dR(iCPs,3);
    %     dRdetaMatrix(3,3*iCPs) = dR(iCPs,3);
    % end
    %
    % dR/dxi
    dRdxiMatrixPage(:,1,noDOFsPerCP*(1:noCPsEl) - noDOFsPerCP + 1) = BSplinePatch.dRdXi(:,iGPs,:);
    dRdxiMatrixPage(:,2,noDOFsPerCP*(1:noCPsEl) - noDOFsPerCP + 2) = BSplinePatch.dRdXi(:,iGPs,:);
    dRdxiMatrixPage(:,3,noDOFsPerCP*(1:noCPsEl) - noDOFsPerCP + 3) = BSplinePatch.dRdXi(:,iGPs,:);
    
    % dR/deta
    dRdetaMatrixPage(:,1,noDOFsPerCP*(1:noCPsEl) - noDOFsPerCP + 1) = BSplinePatch.dRdEta(:,iGPs,:);
    dRdetaMatrixPage(:,2,noDOFsPerCP*(1:noCPsEl) - noDOFsPerCP + 2) = BSplinePatch.dRdEta(:,iGPs,:);
    dRdetaMatrixPage(:,3,noDOFsPerCP*(1:noCPsEl) - noDOFsPerCP + 3) = BSplinePatch.dRdEta(:,iGPs,:);
    
    %% 3ii. Compute the prestress Voigt vectors pagewise
    prestressActPage = squeeze(BSplinePatch.prestressVoigtVector(:,iGPs,:));
    
    %% 3iii. Get the covariant basis of the reference configuration pagewise
    if isReferenceUpdated
        G1GP = pmtimes(dRdxiMatrixPage,CPVctPage);
        G2GP = pmtimes(dRdetaMatrixPage,CPVctPage);
        G1GPPage = reshape(G1GP,[BSplinePatch.noElmnts,1,3]);
        G2GPPage = reshape(G2GP,[BSplinePatch.noElmnts,1,3]);
    else
        G1GPPage = BSplinePatch.GXi(:,iGPs,:);
        G2GPPage = BSplinePatch.GEta(:,iGPs,:);
    end
    GGPCovariantPage = cat(2,G1GPPage,G2GPPage);
    
    %% 3iv. Compute the surface normal base vector pagewise
    if isReferenceUpdated
        G3GP = cross(G1GP,G2GP);
        G3GP_old = cross(squeeze(BSplinePatch.GXi(:,iGPs,:)),squeeze(BSplinePatch.GEta(:,iGPs,:)));
        dA = sqrt(G3GP(:,1).^2 + G3GP(:,2).^2 + G3GP(:,3).^2);
        dA_old = sqrt(G3GP_old(:,1).^2 + G3GP_old(:,2).^2 + G3GP_old(:,3).^2);
        correctionFactorPage = dA_old.^(-1).*dA;
        elementAreaOnGP = BSplinePatch.elementAreaOnGP(:,iGPs).*correctionFactorPage;
    else
        elementAreaOnGP = BSplinePatch.elementAreaOnGP(:,iGPs);
    end
    
    %% 3v. Compute the covariant basis of the current configuration pagewise
    g1Page = pmtimes(dRdxiMatrixPage,CPdVctPage);
    g2Page = pmtimes(dRdetaMatrixPage,CPdVctPage);

    %% 3vi. Compute element tangent stiffness matrix, the residual internal load vector and the stiffness matrix needed for estimation of the stabilization parameter corresponding to the application of weak boundary conditions with the Nitsche method
    [tangStiffMtxElOnGPPage,resVctElOnGPPage] = computeIGAElTangentStiffMtxResVctMembraneNLinear...
        (dRdxiMatrixPage,dRdetaMatrixPage,GGPCovariantPage,g1Page,g2Page,thickness,DmPage,prestressActPage);
    
    %% 3vii. Add the contributions from the Gauss Point
    tangStiffMtxElPage = tangStiffMtxElPage + pstimes(tangStiffMtxElOnGPPage,elementAreaOnGP);
    if ~ischar(loadFactor)
        resVctElPage = resVctElPage + pstimes(resVctElOnGPPage,elementAreaOnGP);
    end
end

%% 4. Assemble to the global tangent matrix and the stiffness matrix needed for estimation of the stabilization parameter corresponding to the application of weak boundary conditions with the Nitsche method
tanStiffMtx = assembleSparseMatricies(BSplinePatch.EFT,noDOFs,noDOFsEl,tangStiffMtxElPage);

%% 5. Add existing tangent matrix from the application of follower loads
if ~ischar(tanMtxLoad)
    tanStiffMtx = tanStiffMtx + tanMtxLoad;
end

%% 6. Assemble to the global residual vector
if ~ischar(loadFactor)
    resVct = assembleVectors(BSplinePatch.EFT,noDOFs,noDOFsEl,resVctElPage);
else
    resVct = 'undefined';
end

%% 7. Update the residual load vector to account for the external loading
if ~ischar(loadFactor)
    resVct = resVct - loadFactor*BSplinePatch.FGamma;
end

%% 8. Updated the tangent stiffness matrix and the residual vector with the contributions resulting from the constant matrix
if ~ischar(constMtx)
    tanStiffMtx = tanStiffMtx + constMtx;
end
if ~ischar(constMtx) && ~ischar(loadFactor)
    resVct = resVct + constMtx*dHat;
end

%% 9. Compute the tangent stiffness matrix and residual load vector resulting from embedded into the patch cables
if isfield(BSplinePatch,'cables')
    if BSplinePatch.cables.No > 0
        [tanMtxCables,resVctCables] = ...
            computeTangentStiffMtxResVctCablesInThinStructureAnalysis...
            (BSplinePatch,noDOFs);
    end
else
    tanMtxCables = 'undefined';
    resVctCables = 'undefined';
end
if isfield(BSplinePatch,'cables')
    if BSplinePatch.cables.No > 0
        tanStiffMtx = tanStiffMtx + tanMtxCables;
    end
end
if isfield(BSplinePatch,'cables') && ~ischar(loadFactor)
    if BSplinePatch.cables.No > 0
        resVct = resVct + resVctCables;
    end
end

%% 10. Compute the tangent matrix and residual vector contributions due to the application of weak boundary conditions
if isWeakDBCTangent
    [tanMtxWeakDBC,resVctWeakDBC,BSplinePatch] = computeTangMtxResVct...
        (BSplinePatch,dHat,connections,noDOFs,propCoupling,tanStiffMtx,...
        noPatch,noTimeStep,noNonlinearIteration,noWeakDBCCnd,...
        thickness,t,propTransientAnalysis,tab,outMsg);
else
    tanMtxWeakDBC = 'undefined';
    resVctWeakDBC = 'undefined';
end
if ~ischar(tanMtxWeakDBC)
    tanStiffMtx = tanStiffMtx + tanMtxWeakDBC;
end
if ~ischar(resVctWeakDBC) && ~ischar(loadFactor)
    resVct = resVct + resVctWeakDBC;
end

%% 11. Check output
if isBSplinePatchCell
    BSplinePatch = {BSplinePatch};
end

end
