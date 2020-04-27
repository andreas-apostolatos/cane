function [tanMtxCables, resVctCables] = ...
    computeTangentStiffMtxResVctCablesInThinStructureAnalysis ...
    (BSplinePatch, noDOFs)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangent stiffness matrix and residual load vector
% corresponding to cables which are embedded into a B-Spline patch
% representing a thin-walled structure in isogeometric analysis.
%
%                Input :
%         BSplinePatch : The B-Spline patch array containing:
%                       .p,q : The polynomial orders of the B-Spline 
%                              surface in both parametric directions
%                    .Xi,Eta : The knot vectors in both parametric
%                              directions
%                        .CP : The set of control points and weights
%                   .isNURBS : Flag on whether the basis is a NURBS or a
%                              B-Spline
%                .parameters : Technical parameters for the structure
%                   .homDOFs : The global numbering of the DOFs where 
%                              homogeneous Dirichlet boundary conditions 
%                              are applied
%                 .inhomDOFs : The global numbering of the DOFs where
%                              inhomogeneous Dirichlet boundary conditions
%                              are applied
%           .valuesInhomDOFs : The values on the DOFs corresponding to the
%                              application of inhomogeneous Dirichlet
%                              boundary conditions
%                       .NBC : Structure containing information on the 
%                              application of the Neumann boundary 
%                              conditions
%                                        .noCnd : Number of Neumann 
%                                                 boundary conditions
%                              .xiLoadExtension : Cell array {.noCnd} 
%                                                 containing the load 
%                                                 extensions in the xi-
%                                                 direction
%                             .etaLoadExtension : Cell array {.noCnd} 
%                                                 containing the load 
%                                                 extensions in the eta-
%                                                 direction
%                                .loadAmplitude : Array (1,.noCnd) 
%                                                 containing the load 
%                                                 amplitudes
%                                .loadDirection : Array (1,.noCnd) 
%                                                 containing the load 
%                                                 directions
%                               .computeLoadVct : Cell array {.noCnd} 
%                                                 containing the function 
%                                                 name for the computation 
%                                                 of the load vector
%                               .isConservative : Array (1,.noCnd) of flags 
%                                                 indicating whether the 
%                                                 load is conservative or 
%                                                 not
%                   .weakDBC : Structure containing information on the
%                              application of weak boundary conditions:
%                                        .alpha : The penalty parameter
%                                  .xiExtension : xi-extension of the
%                                                 boundary where the weak
%                                                 boundary conditions are
%                                                 applied
%                                 .etaExtension : eta-extension of the
%                                                 boundary where the weak
%                                                 boundary conditions are
%                                                 applied
%                    .cables : Array related to cables that are embedded 
%                              onto the B-Spline patch containing the 
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
%               noDOFs : The number of DOFs for the B-Spline patch
%
%               Output :
%         tanMtxCables : The contribution of the tangent stiffness matrix
%                        due to the embedded cables
%         resVctCables : The residual load vector due to the embedded
%                        cables
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the cables that are ambedded in the patch
% ->
%    1i. Get the extensions of the boundary where the cable is embedded and the corresponding parameters
%
%   1ii. Get the running and the fixed parameters on the patch boundary where the cables is embedded
%
%  1iii. Get the parameter space of the cable
%
%   1iv. Issue Gauss Point coordinates and weights
%
%    1v. Loop over the elements on the parameter space where the cable is embedded
%    ->
%        1v.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        1v.2. Loop over the Gauss points
%        ->
%              1v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%             1v.2ii. Compute the NURBS basis functions and their derivatives
%
%            1v.2iii. Create the element freedom table
%
%             1v.2iv. Compute the matrices containing the basis functions and their derivatives
%
%              1v.2v. Compute the base vectors of the reference and the current configuration along the cable
%
%             1v.2vi. Compute the Green-Lagrange strain in the contravariant basis
%
%            1v.2vii. Transform the Green-Lagrange strain in the local Cartesian basis
%
%           1v.2viii. Compute the 2nd Piola-Kirchhoff stress in the local Cartesian basis
%
%             1v.2ix. Compute the first variation of the strain with respect to the DOFs in the contravariant basis pagewise
%
%              1v.2x. Transform the first variation of the strain with respect to the DOFs in the local Cartesian system
%
%             1v.2xi. Compute the first variation of the 2nd Piola-Kirchhoff stress wrt the DOFs in the local Cartesian basis
%
%            1v.2xii. Compute the second variation of the strain with respect to the DOFs in the contravariant basis
%
%           1v.2xiii. Transform the second variation of the strain with respect to the DOFs in the local Cartesian basis
%
%            1v.2xiv. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%             1v.2xv. Compute the element length at the GP
%
%            1v.2xvi. Compute the tangent stiffness matrix at the Gauss Point and add it to the global matrix
%
%           1v.2xvii. Compute the residual vector at the Gauss Point and add it to the global matrix
%        <-
%    <-
% <-
%
%% Function main body

%% 0. Read input

% Check the NURBS geometry input
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
end

% Reassign the analysis arrays
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
CPd = BSplinePatch.CPd;
isNURBS = BSplinePatch.isNURBS;

% Get the DOF numbering
DOFNumbering = BSplinePatch.DOFNumbering;

% Number of Control Points in xi-,eta- directions
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of local DOFs
noCPsEl = (p+1)*(q+1);
noDOFsEl = 3*noCPsEl;

% Initialize the element freedom table
EFT = zeros(1,noDOFsEl);

% Initialize auxiliary arrays
dRdxiTildeMtx = zeros(3,noDOFsEl);

% Initialize the output arrays
tanMtxCables = zeros(noDOFs,noDOFs);
resVctCables = zeros(noDOFs,1);

%% 1. Loop over all the cables that are ambedded in the patch
for counterCables = 1:BSplinePatch.cables.No
    %% 1i. Get the extensions of the boundary where the cable is embedded and the corresponding parameters
    xiExtension = BSplinePatch.cables.xiExtension{counterCables};
    etaExtension = BSplinePatch.cables.etaExtension{counterCables};
    parametersCable = BSplinePatch.cables.parameters{counterCables};
    areaCS = parametersCable.areaCS;
    prestress = parametersCable.prestress;
    
    %% 1ii. Get the running and the fixed parameters on the patch boundary where the cables is embedded
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
    
    %% 1iii. Get the parameter space of the cable
    cableRegionOnKnotVector = unique(cableRegionOnKnotVector);
    
    %% 1iv. Issue Gauss Point coordinates and weights
    if strcmp(BSplinePatch.cables.int.type,'default')
        if isOnXi
            pDegree = p + 1;
        else
            pDegree = q + 1;
        end
        noGPs = ceil((pDegree + 1)/2);
    elseif strcmp(BSplinePatch.cables.int.type,'user')
        noGPs = BSplinePatch.cables.int.noGPs;
    end
    [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);
    
    %% 1v. Loop over the elements on the parameter space where the cable is embedded
    for i = 1:length(cableRegionOnKnotVector)-1
        if cableRegionOnKnotVector(i) ~= cableRegionOnKnotVector(i+1)
            %% 1v.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
            detJxizeta = (cableRegionOnKnotVector(i+1) - cableRegionOnKnotVector(i))/2;
            
            %% 1v.2. Loop over the Gauss points
            for j = 1:noGPs
                %% 1v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                xiEta = ((1-GP(j))*cableRegionOnKnotVector(i) + (1+GP(j))*cableRegionOnKnotVector(i+1))/2;

                %% 1v.2ii. Compute the NURBS basis functions and their derivatives
                if isOnXi
                    xi = xiEta;
                    xiSpan = findKnotSpan(xi,Xi,nxi);
                else
                    eta = xiEta;
                    etaSpan = findKnotSpan(eta,Eta,neta);
                end
                dR = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);
                
                %% 1v.2iii. Create the element freedom table
                r = 1;
                for cpj = etaSpan-q:etaSpan
                    for cpi = xiSpan-p:xiSpan
                        EFT(r)   = DOFNumbering(cpi,cpj,1);
                        EFT(r+1) = DOFNumbering(cpi,cpj,2);
                        EFT(r+2) = DOFNumbering(cpi,cpj,3);
                        r = r + 3;
                    end
                end

                %% 1v.2iv. Compute the matrices containing the basis functions and their derivatives
                for iCPs = 1:noCPsEl
                    if isOnXi
                        deriv = 2;
                    else
                        deriv = 3;
                    end
                    dRdxiTildeMtx(1,3*iCPs - 2) = dR(iCPs,deriv);
                    dRdxiTildeMtx(2,3*iCPs - 1) = dR(iCPs,deriv);
                    dRdxiTildeMtx(3,3*iCPs) = dR(iCPs,deriv);
                end
                
                %% 1v.2v. Compute the base vectors of the reference and the current configuration along the cable
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CP,0,dR);
                [a1,a2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CPd,0,dR);
                if isOnXi
                    A = A1;
                    a = a1;
                else
                    A = A2;
                    a = a2;
                end
                
                %% 1v.2vi. Compute the Green-Lagrange strain in the contravariant basis
                EpsilonContravariant = .5*(a'*a - A'*A);
                
                %% 1v.2vii. Transform the Green-Lagrange strain in the local Cartesian basis
                EpsilonLC = 1/(A'*A)*EpsilonContravariant;
                
                %% 1v.2viii. Compute the 2nd Piola-Kirchhoff stress in the local Cartesian basis
                NLC = prestress*areaCS + parametersCable.E*areaCS*EpsilonLC;
                
                %% 1v.2ix. Compute the first variation of the strain with respect to the DOFs in the contravariant basis pagewise
                dEpsilonContravariant = a'*dRdxiTildeMtx;
                
                %% 1v.2x. Transform the first variation of the strain with respect to the DOFs in the local Cartesian system
                dEpsilonLC = 1/(A'*A)*dEpsilonContravariant;
                
                %% 1v.2xi. Compute the first variation of the 2nd Piola-Kirchhoff stress wrt the DOFs in the local Cartesian basis
                dNLC = parametersCable.E*areaCS*dEpsilonLC;
                
                %% 1v.2xii. Compute the second variation of the strain with respect to the DOFs in the contravariant basis
                ddEpsilonContravariant = dRdxiTildeMtx'*dRdxiTildeMtx;
                
                %% 1v.2xiii. Transform the second variation of the strain with respect to the DOFs in the local Cartesian basis
                ddEpsilonLC = 1/(A'*A)*ddEpsilonContravariant;
                
                %% 1v.2xiv. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                detJxxi = norm(A);
                
                %% 1v.2xv. Compute the element length at the GP
                elementLengthOnGP = detJxxi*detJxizeta*GW(j);
                
                %% 1v.2xvi. Compute the tangent stiffness matrix at the Gauss Point and add it to the global matrix
                tanMtxCables(EFT,EFT) = tanMtxCables(EFT,EFT) + ...
                    (dEpsilonLC'*dNLC + ddEpsilonLC*NLC)*elementLengthOnGP;
                
                %% 1v.2xvii. Compute the residual vector at the Gauss Point and add it to the global matrix
                resVctCables(EFT,1) = resVctCables(EFT,1) + ...
                    dEpsilonLC'*NLC*elementLengthOnGP;
            end
        end
    end
end
