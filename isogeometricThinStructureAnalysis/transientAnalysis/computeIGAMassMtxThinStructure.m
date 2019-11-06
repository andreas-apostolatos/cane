function massMtx = computeIGAMassMtxThinStructure...
    (BSplinePatches,noDOFs)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the mass matrix corresponding to multipatch isogeometric thin
% structure which is assumed to be constant throughout the simulation.
%
%                 input :
%              constMtx : The constant matrix of the system
%        BSplinePatches : Its an array of structures {patch1,patch2,...} each 
%                         of the patch structures containing the
%                            following information
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
%                         .DOFNumbering: Numbering of the DOFs sorted into
%                                        a 3D array
%                           .parameters: material parameters of the shell
%                                  .int: On the numerical integration
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
%               .cables : Array related to cables that are embedded onto
%                         the B-Spline patch containing the 
%                              following information :
%                                   .noCables : Number of embedded cables
%                                .xiExtension : The xi-extension of the 
%                                               embedded cables in the 
%                                               B-Spline patch
%                               .etaExtension : The eta-extension of the 
%                                               embedded cables in the 
%                                               B-Spline patch
%                            .parametersCable : Array containing the 
%                                               technical and geometrical 
%                                               parameters that a cable
%                                               must have :
%                                                       .E : The Young's 
%                                                            modulus of
%                                                            the cable
%                                                .radiusCS : The cross 
%                                                            sectional
%                                                            radius of the 
%                                                            cable
%                                                     .rho : The density of 
%                                                            the cable
%                                               .prestress : The prestress 
%                                                            of the cable
%                                                     .int : On the
%                                                            integration of
%                                                            the cables
%           connections : Define the connection between the patches:
%                            .No : Number of connections
%                     .xiEtaCoup : [patchID1 patchID2 xi12 eta12 xi21 eta21
%                                   ...      ...    ...   ...  ...   ...]
%                noDOFs : Number of DOFs for the mulitpatch system
%                         including also Lagrange Multipliers DOFs if they
%                         are employed
%          propCoupling : Coupling properties for the multipatch geometries
%                           .method : Coupling method
%
%                output :
%               massMtx : The mass matrix corresponding to the multipatch
%                         thin isogeometric structure
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the patches
% ->
%    1i. Read input for the current patch
%
%   1ii. Choose an integration rule
%
% 1iii. loops over elements
%       ->
%           1iii.1. Create an element freedom table for the current patch
%
%           1iii.2. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1]
%
%           1iii.3. Loop over all Gauss points
%                   ->
%                       1iii.3i. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
%
%                      1iii.3ii. Compute the NURBS basis functions and up to their first derivatives at the Gauss Point
%
%                     1iii.3iii. Compute the covariant base vectors of the reference configuration
%
%                      1iii.3iv. Compute the surface normal of the reference configuration (third covariant base vector not normalized)
%
%                       1iii.3v. Compute the legth of G3Tilde (= area dA of the undeformed configuration)
%
%                      1iii.3vi. Form the basis function matrix
%
%                     1iii.3vii. Compute the mass matrix at the Gauss Point and add the contribution
%                   <-
%       <-
%
%  1iv. Loop over all the cables that are ambedded in the patch
%  ->
%       1iv.1. Get the extensions of the boundary where the cable is embedded and the corresponding parameters
%
%       1iv.2. Get the running and the fixed parameters on the patch boundary where the cables is embedded
%
%       1iv.3. Get the parameter space of the cable
%
%       1iv.4. Issue Gauss Point coordinates and weights
%
%       1iv.5. Loop over the elements on the parameter space where the cable is embedded
%       ->
%              1iv.5i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%              1iv.5ii. Loop over the Gauss points
%              ->
%                       1iv.5ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%                       1iv.5ii.2. Compute the NURBS basis functions and their derivatives
%
%                       1iv.5ii.3. Create the element freedom table
%
%                       1iv.5ii.4. Compute the matrices containing the basis functions
%
%                       1iv.5ii.5. Compute the base vectors of the reference configuration along the cable
%
%                       1iv.5ii.6. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%                       1iv.5ii.7. Compute the element length at the GP
%
%                       1iv.5ii.8. Compute the mass matrix at the Gauss point and add it to the global matrix
%              <-
%       <-
%  <-
% <-
%
% 2. Loop over all the connections in case the mortar method is used for the multipatch coupling
% ->
%    2i. Get the ID of the slave patch
%
%   2ii. Get the numbering of the necessary Lagrange Multilpiers, interface and domain DOFs
%
%  2iii. Get the mortar transformation matrix for the given patch connection
%
%   2iv. Re-arrange the stiffness entries by elimination of the Lagrange Multipliers DOFs  
% <-
%
%% Function main body

%% 0. Read input

% Number of patches
noPatches = length(BSplinePatches);

% Initialize output array
massMtx = zeros(noDOFs,noDOFs);

%% 1. Loop over all the patches
for iPatches = 1:noPatches
    %% 1i. Read input for the current patch
    p = BSplinePatches{iPatches}.p;
    q = BSplinePatches{iPatches}.q;
    Xi = BSplinePatches{iPatches}.Xi;
    Eta = BSplinePatches{iPatches}.Eta;
    CP = BSplinePatches{iPatches}.CP;
    isNURBS = BSplinePatches{iPatches}.isNURBS;
    DOFNumbering = BSplinePatches{iPatches}.DOFNumbering;
    parameters = BSplinePatches{iPatches}.parameters;
    rho = parameters.rho;
    thickness = parameters.t;
    int = BSplinePatches{iPatches}.int;
    mxi = length(Xi);
    meta = length(Eta);
    nxi = length(CP(:,1,1));
    neta = length(CP(1,:,1));
    checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);
    noDOFsPatch = BSplinePatches{iPatches}.noDOFs;
    noCPsEl = (p + 1)*(q + 1);
    noDOFsEl = 3*noCPsEl;
    EFT = zeros(1,noDOFsEl);
    RMtx = zeros(3,noDOFsEl);
    massMtxPatch = zeros(noDOFsPatch,noDOFsPatch);
    if isfield(BSplinePatches{iPatches},'cables')
        if BSplinePatches{iPatches}.cables.No > 0
            MassMtxPatchCables = zeros(noDOFsPatch,noDOFsPatch);
        end
    end
    
    %% 1ii. Choose an integration rule

    % Select the integration scheme
    if strcmp(int.type,'default')
        noGPXi = p + 1;
        noGPEta = q + 1;
    elseif strcmp(int.type,'user')
        noGPXi = int.xiNGP;
        noGPEta = int.etaNGP;
    end

    % Issue the Gauss Point coordinates and weights
    [xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(noGPXi);
    [etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(noGPEta);
    
    %% 1iii. loops over elements
    for iSpanEta = q + 1:meta - q - 1
        for iSpanXi = p + 1:mxi - p - 1
            % check if element is greater than zero
            if (Xi(iSpanXi + 1) ~= Xi(iSpanXi) && Eta(iSpanEta + 1) ~= Eta(iSpanEta))
                %% 1iii.1. Create an element freedom table for the current patch
                k = 1;
                for iCPEta = iSpanEta - q:iSpanEta
                    for iCPXi = iSpanXi - p:iSpanXi
                        EFT(k) = DOFNumbering(iCPXi,iCPEta,1);
                        EFT(k + 1) = DOFNumbering(iCPXi,iCPEta,2);
                        EFT(k + 2) = DOFNumbering(iCPXi,iCPEta,3);
                        k = k + 3;
                    end
                end
                
                %% 1iii.2. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1]
                %
                %         | xi_i+1 - xi_i                    |
                %         | -------------            0       |
                %         |        2                         |
                %  xi,u = |                                  |
                %         |                  eta_j+1 - eta_j |
                %         |        0         --------------- |
                %         |                          2       |
                detJxiu = (Xi(iSpanXi + 1) - Xi(iSpanXi))*(Eta(iSpanEta + 1) - Eta(iSpanEta))/4;
                
                %% 1iii.3. Loop over all Gauss points
                for iGPEta = 1:noGPEta
                    for iGPXi = 1:noGPXi
                        %% 1iii.3i. Compute the NURBS coordinates xi,eta of the Gauss Point coordinates in the bi-unit interval [-1, 1]
                        xi = (Xi(iSpanXi + 1) + Xi(iSpanXi) + xiGP(iGPXi)*(Xi(iSpanXi + 1) - Xi(iSpanXi)))/2;
                        eta = (Eta(iSpanEta + 1) + Eta(iSpanEta) + etaGP(iGPEta)*(Eta(iSpanEta + 1) - Eta(iSpanEta)))/2;
                        
                        %% 1iii.3ii. Compute the NURBS basis functions and up to their first derivatives at the Gauss Point
                        nDrvBasis = 1;
                        dR = computeIGABasisFunctionsAndDerivativesForSurface...
                            (iSpanXi,p,xi,Xi,iSpanEta,q,eta,Eta,CP,isNURBS,nDrvBasis);
                        
                        %% 1iii.3iii. Compute the covariant base vectors of the reference configuration
                        nDrvBaseVct = 0;
                        [G1,G2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                            (iSpanXi,p,iSpanEta,q,CP,nDrvBaseVct,dR);
                        
                        %% 1iii.3iv. Compute the surface normal of the reference configuration (third covariant base vector not normalized)
                        G3Tilde = cross(G1(:,1),G2(:,1));
                        
                        %% 1iii.3v. Compute the legth of G3Tilde (= area dA of the undeformed configuration)
                        dA = norm(G3Tilde);
                        
                        %% 1iii.3vi. Form the basis function matrix
                        RMatrixEl = zeros(3,noDOFsEl);
                        for iCP = 1:noCPsEl
                            RMatrixEl(1,3*iCP - 2) = dR(iCP,1);
                            RMatrixEl(2,3*iCP - 1) = dR(iCP,1);
                            RMatrixEl(3,3*iCP) = dR(iCP,1);
                        end
                        
                        %% 1iii.3vii. Compute the mass matrix at the Gauss Point and add the contribution
                        massMtxPatch(EFT,EFT) = massMtxPatch(EFT,EFT) + ...
                            rho*(RMatrixEl'*RMatrixEl)*...
                            thickness*dA*detJxiu*xiGW(iGPXi)*etaGW(iGPEta);
                    end
                end
            end
        end
    end
    
    %% 1iv. Loop over all the cables that are ambedded in the patch
    for counterCables = 1:BSplinePatches{iPatches}.cables.No
        %% 1iv.1. Get the extensions of the boundary where the cable is embedded and the corresponding parameters
        xiExtension = BSplinePatches{iPatches}.cables.xiExtension{counterCables};
        etaExtension = BSplinePatches{iPatches}.cables.etaExtension{counterCables};
        parametersCable = BSplinePatches{iPatches}.cables.parameters{counterCables};
        areaCS = parametersCable.areaCS;
        rhoCable = parametersCable.rho;
        
        %% 1iv.2. Get the running and the fixed parameters on the patch boundary where the cables is embedded
        if etaExtension(1) == etaExtension(2)
            % Coupled region in xi-direction
            couplingRegion = xiExtension;

            % Find the correct spans for the coupled region
            spanStart = findKnotSpan(couplingRegion(1),Xi,nxi);
            spanEnd = findKnotSpan(couplingRegion(2),Xi,nxi) + 1;

            % Corresponding to the coupled region knot span
            cableRegionOnKnotVector = Xi(spanStart:spanEnd);

            % Fixed parameter on the parametric net
            eta = etaExtension(1);

            % Find the span where xiEta it lies in
            etaSpan = findKnotSpan(eta,Eta,neta);

            % Flag on whether the coupling line is over xi
            isOnXi = true;
        else
            % Coupled region in eta-direction
            couplingRegion = etaExtension;

            % Find the correct spans for the coupled region
            spanStart = findKnotSpan(couplingRegion(1),Eta,neta);   
            spanEnd = findKnotSpan(couplingRegion(2),Eta,neta) + 1;

            % Corresponding to the coupled region knot span
            cableRegionOnKnotVector = Eta(spanStart:spanEnd);

            % Fixed parameter on the parametric net
            xi = xiExtension(1);

            % Find the span where uv it lies in
            xiSpan = findKnotSpan(xi,Xi,nxi);

            % Flag on whether the coupling line is over eta
            isOnXi = false;
        end
        
        %% 1iv.3. Get the parameter space of the cable
        cableRegionOnKnotVector = unique(cableRegionOnKnotVector);
        
        %% 1iv.4. Issue Gauss Point coordinates and weights
        if strcmp(BSplinePatches{iPatches}.cables.int.type,'default')
            if isOnXi
                pDegree = p + 1;
            else
                pDegree = q + 1;
            end
            noGPs = ceil((pDegree + 1)/2);
        elseif strcmp(BSplinePatches{iPatches}.cables.int.type,'user')
            noGPs = BSplinePatches{iPatches}.cables.int.noGPs;
        end
        [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);
        
        %% 1iv.5. Loop over the elements on the parameter space where the cable is embedded
        for iSpan = 1:length(cableRegionOnKnotVector)-1
            if cableRegionOnKnotVector(iSpan) ~= cableRegionOnKnotVector(iSpan + 1)
                %% 1iv.5i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
                detJxizeta = (cableRegionOnKnotVector(iSpan + 1) - cableRegionOnKnotVector(iSpan))/2;
                
                %% 1iv.5ii. Loop over the Gauss points
                for iGPs = 1:noGPs
                    %% 1iv.5ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                    xiEta = ((1 - GP(iGPs))*cableRegionOnKnotVector(iSpan) + (1 + GP(iGPs))*cableRegionOnKnotVector(iSpan + 1))/2;
                    
                    %% 1iv.5ii.2. Compute the NURBS basis functions and their derivatives
                    if isOnXi
                        xi = xiEta;
                        xiSpan = findKnotSpan(xi,Xi,nxi);
                    else
                        eta = xiEta;
                        etaSpan = findKnotSpan(eta,Eta,neta);
                    end
                    dR = computeIGABasisFunctionsAndDerivativesForSurface...
                        (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);
                    
                    %% 1iv.5ii.3. Create the element freedom table
                    r = 1;
                    for iCPEta = etaSpan - q:etaSpan
                        for iCPXi = xiSpan - p:xiSpan
                            EFT(r)   = DOFNumbering(iCPXi,iCPEta,1);
                            EFT(r + 1) = DOFNumbering(iCPXi,iCPEta,2);
                            EFT(r + 2) = DOFNumbering(iCPXi,iCPEta,3);
                            r = r + 3;
                        end
                    end
                    
                    %% 1iv.5ii.4. Compute the matrices containing the basis functions
                    for iCPs = 1:noCPsEl
                        RMtx(1,3*iCPs - 2) = dR(iCPs,1);
                        RMtx(2,3*iCPs - 1) = dR(iCPs,1);
                        RMtx(3,3*iCPs) = dR(iCPs,1);
                    end
                    
                    %% 1iv.5ii.5. Compute the base vectors of the reference configuration along the cable
                    [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                        (xiSpan,p,etaSpan,q,CP,0,dR);
                    
                    %% 1iv.5ii.6. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                    if isOnXi
                        detJxxi = norm(A1(:,1));
                    else
                        detJxxi = norm(A2(:,1));
                    end
                    
                    %% 1iv.5ii.7. Compute the element length at the GP
                    elementLengthOnGP = detJxxi*detJxizeta*GW(iGPs);
                    
                    %% 1iv.5ii.8. Compute the mass matrix at the Gauss point and add it to the global matrix
                    MassMtxPatchCables(EFT,EFT) = MassMtxPatchCables(EFT,EFT) + ...
                        rhoCable*areaCS*(RMtx'*RMtx)*elementLengthOnGP;
                end
            end
        end
    end
    
    %% 1v. Assemble to the mass matrix of the mulipatch system
    if isfield(BSplinePatches{iPatches},'cables')
        if BSplinePatches{iPatches}.cables.No > 0
            massMtxPatch = massMtxPatch + MassMtxPatchCables;
        end
    end
    massMtx(BSplinePatches{iPatches}.EFTPatches,BSplinePatches{iPatches}.EFTPatches) = massMtxPatch;
end

end
