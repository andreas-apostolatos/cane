function QMtx = computeWeakDBCMtxStabilEigenvalueProblemNitscheIGAMembrane...
    (BSplinePatch,xiExtension,etaExtension)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function main body
%
% Returns the Q-matrix resulting from the discretization of the L2-norm of
% the traction forces corresponding to the linear membrane problem needed
% for setting up the eigenvalue problem which's solution can be used for
% the estimation of the stabilization parameter for the weak imposition of
% Dirichlet boundary conditions in the isogeometric membrane, see,
%
% "Apostolatos et al. A nitsche-type formulation and comparison of the 
% most common domain decomposition methods in isogeometric analysis, 
% International Journal for Numerical Methods in Engineering 97 (7) (2014) 
% 473--504."
%
%               Input :
%        BSplinePatch : The B-Spline patch array containing:
%                      .p,q : The polynomial orders of the B-Spline surface
%                             in both parametric directions
%                   .Xi,Eta : The knot vectors in both parametric
%                             directions
%                       .CP : The set of control points and weights
%                  .isNURBS : Flag on whether the basis is a NURBS or a
%                             B-Spline
%               .parameters : Technical parameters for the structure
%                  .homDOFs : The global numbering of the DOFs where 
%                             homogeneous Dirichlet boundary conditions are
%                             applied
%                .inhomDOFs : The global numbering of the DOFs where
%                             inhomogeneous Dirichlet boundary conditions
%                             are applied
%          .valuesInhomDOFs : The values on the DOFs corresponding to the
%                             application of inhomogeneous Dirichlet
%                             boundary conditions
%                      .NBC : Structure containing information on the 
%                             application of the Neumann boundary 
%                             conditions
%                                       .noCnd : Number of Neumann boundary 
%                                                conditions
%                             .xiLoadExtension : Cell array {.noCnd} 
%                                                containing the load 
%                                                extensions in the xi-
%                                                direction
%                            .etaLoadExtension : Cell array {.noCnd} 
%                                                containing the load 
%                                                extensions in the eta-
%                                                direction
%                               .loadAmplitude : Array (1,.noCnd) 
%                                                containing the load 
%                                                amplitudes
%                               .loadDirection : Array (1,.noCnd) 
%                                                containing the load 
%                                                directions
%                              .computeLoadVct : Cell array {.noCnd} 
%                                                containing the function 
%                                                name for the computation 
%                                                of the load vector
%                              .isConservative : Array (1,.noCnd) of flags 
%                                                indicating whether the 
%                                                load is conservative or 
%                                                not
%                  .weakDBC : Structure containing information on the
%                             application of weak boundary conditions:
%                                       .alpha : The penalty parameter
%                                 .xiExtension : xi-extension of the
%                                                boundary where the weak
%                                                boundary conditions are
%                                                applied
%                                .etaExtension : eta-extension of the
%                                                boundary where the weak
%                                                boundary conditions are
%                                                applied
%         xiExtension : The extension of the Dirichlet boundary in
%                       xi-direction
%        etaExtension : The extension of the Dirichlet boundary in
%                       eta-direction
%
% Function layout :
%
% 0. Read input
%
% 1. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
%
% 2. Get the parameter space of the application of the weak boundary conditions using the knot vector information
%
% 3. Issue Gauss Point coordinates and weights
%
% 4. Loop over the elements on the parameter space where the weak boundary conditions are applied
% ->
%    4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%   4ii. Loop over the Gauss points
%   ->
%        4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%        4ii.2. Compute the NURBS basis functions
%
%        4ii.3. Create the element freedom table
%
%        4ii.4. Compute the matrices containing the basis functions and their derivatives
%
%        4ii.5. Compute the covariant base vectors of the reference configuration
%
%        4ii.6. Compute the normal to the boundary vector and transform it into the contavariant basis
%
%        4ii.7. Compute the contravariant basis
%
%        4ii.8. Compute the local Cartesian basis
%
%        4ii.9. Compute the transformation matrix from the contravariant basis to the local Cartesian one
%
%        4ii.10. Compute the transformation matrix from the local Cartesian basis to the covariant one
%
%        4ii.11. Compute the 2nd Piola-kirchhoff stress in the local Cartesian system
%
%        4ii.12.  Transform the 2nd Piola-kirchhoff stress in the covariant system
%
%        4ii.13. Compute the first variation of the Green-Lagrange strain in the contravariant basis with respect to the DOFs
%
%        4ii.14. Transform the first variation of the Green-Lagrange strain at the local Cartesian basis
%
%        4ii.15. Compute the first variation of the 2nd Piola-Kichhoff stress in the local Cartesian basis
%
%        4ii.16. Transform the first variation of the 2nd Piola-Kichhoff stress at the covariant basis
%
%        4ii.17. Compute the first variation of the traction vector
%
%        4ii.18. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%        4ii.19. Compute the element length at the GP
%
%        4ii.20. Compute the element Q-matrix contribution corresponding to the eigenvalue problem for the estimation of the stabilization parameter for the weak application of the Dirichlet boundary conditions with the Nitsche method and add it to the global matrix
%   <-
%
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
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;
prestress = parameters.prestress;

% Get the DOF numbering
DOFNumbering = BSplinePatch.DOFNumbering;

% Compute the material matrix
Dm = parameters.E*parameters.t/(1-parameters.nue^2)*...
     [1              parameters.nue 0
      parameters.nue 1              0
      0              0              (1-parameters.nue)/2];

% Number of Control Points in xi-,eta- directions
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of local DOFs
noCPsEl = (p+1)*(q+1);
noDOFsEl = 3*noCPsEl;

% Number of DOFs
noDOFs = 3*nxi*neta;

% Initialize the element freedom table
EFT = zeros(1,noDOFsEl);

% Initialize auxiliary arrays
RMtx = zeros(3,noDOFsEl);
dRdxiMtx = zeros(3,noDOFsEl);
dRdetaMtx = zeros(3,noDOFsEl);

% Initialize the output arrays
QMtx = zeros(noDOFs,noDOFs);

%% 1. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
if etaExtension(1) == etaExtension(2)
    % Coupled region in xi-direction
    couplingRegion = xiExtension;

    % Find the correct spans for the coupled region
    spanStart = findKnotSpan(couplingRegion(1),Xi,nxi);
    spanEnd = findKnotSpan(couplingRegion(2),Xi,nxi) + 1;

    % Corresponding to the coupled region knot span
    weakDBCRgionOnKnotVector = Xi(spanStart:spanEnd);

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
    weakDBCRgionOnKnotVector = Eta(spanStart:spanEnd);

    % Fixed parameter on the parametric net
    xi = xiExtension(1);

    % Find the span where uv it lies in
    xiSpan = findKnotSpan(xi,Xi,nxi);

    % Flag on whether the coupling line is over eta
    isOnXi = false;
end

%% 2. Get the parameter space of the application of the weak boundary conditions using the knot vector information
weakDBCRgionOnKnotVector = unique(weakDBCRgionOnKnotVector);

%% 3. Issue Gauss Point coordinates and weights
if strcmp(BSplinePatch.weakDBC.int.type,'default')
    if isOnXi
        pDegree = p + 1;
    else
        pDegree = q + 1;
    end
    noGPs = ceil((pDegree + 1)/2);
elseif strcmp(BSplinePatch.weakDBC.int.type,'user')
    noGPs = BSplinePatch.weakDBC.int.noGPs;
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);

%% 4. Loop over the elements on the parameter space where the weak boundary conditions are applied
for i = 1:length(weakDBCRgionOnKnotVector)-1
    if weakDBCRgionOnKnotVector(i) ~= weakDBCRgionOnKnotVector(i+1)
        %% 4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
        detJxizeta = (weakDBCRgionOnKnotVector(i+1)-weakDBCRgionOnKnotVector(i))/2;
        
        %% 4ii. Loop over the Gauss points
        for j = 1:noGPs
            %% 4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
            xiEta = ((1-GP(j))*weakDBCRgionOnKnotVector(i) + (1+GP(j))*weakDBCRgionOnKnotVector(i+1))/2;
            
            %% 4ii.2. Compute the NURBS basis functions
            if isOnXi
                xi = xiEta;
                xiSpan = findKnotSpan(xi,Xi,nxi);
            else
                eta = xiEta;
                etaSpan = findKnotSpan(eta,Eta,neta);
            end
            dR = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);
            
            %% 4ii.3. Create the element freedom table
            r = 1;
            for cpj = etaSpan-q:etaSpan
                for cpi = xiSpan-p:xiSpan
                    EFT(r)   = DOFNumbering(cpi,cpj,1);
                    EFT(r+1) = DOFNumbering(cpi,cpj,2);
                    EFT(r+2) = DOFNumbering(cpi,cpj,3);
                    r = r + 3;
                end
            end
            
            %% 4ii.4. Compute the matrices containing the basis functions and their derivatives
            for iCPs = 1:noCPsEl
                % R
                RMtx(1,3*iCPs - 2) = dR(iCPs,1);
                RMtx(2,3*iCPs - 1) = dR(iCPs,1);
                RMtx(3,3*iCPs) = dR(iCPs,1);

                % dR/dxi
                dRdxiMtx(1,3*iCPs - 2) = dR(iCPs,2);
                dRdxiMtx(2,3*iCPs - 1) = dR(iCPs,2);
                dRdxiMtx(3,3*iCPs) = dR(iCPs,2);

                % dR/deta
                dRdetaMtx(1,3*iCPs - 2) = dR(iCPs,3);
                dRdetaMtx(2,3*iCPs - 1) = dR(iCPs,3);
                dRdetaMtx(3,3*iCPs) = dR(iCPs,3);
            end
            
            %% 4ii.5. Compute the covariant base vectors of the reference configuration
            [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                (xiSpan,p,etaSpan,q,CP,0,dR);
            A3 = cross(A1,A2)/norm(cross(A1,A2));
            
            %% 4ii.6. Compute the normal to the boundary vector and transform it into the contavariant basis
            [uGC,~] = computeNormalAndTangentVectorsToBSplineBoundary...
                (xi,Xi,eta,Eta,A1,A2,A3,isOnXi);
            uContravariant = [A1'
                              A2']*uGC;
                          
            %% 4ii.7. Compute the contravariant basis
            AabCovariant = [A1 A2]'*[A1 A2];
            AContravariant = AabCovariant\[A1 A2]';
            AContravariant = AContravariant';
            
            %% 4ii.8. Compute the local Cartesian basis
            eLC = computeLocalCartesianBasis4BSplineSurface...
                ([A1 A2],AContravariant);
            
            %% 4ii.9. Compute the transformation matrix from the contravariant basis to the local Cartesian one
            TFromContraToLC4VoigtStrain = ...
                computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell...
                (eLC,AContravariant);
            
            %% 4ii.10. Compute the transformation matrix from the local Cartesian basis to the covariant one
            TFromLCToCov = computeTFromLocalCartesian2CovariantBasis4BSplineSurface...
                (eLC,AContravariant);
            
            %% 4ii.11. Compute the 2nd Piola-kirchhoff stress in the local Cartesian system
            NLC = prestress;
            
            %% 4ii.12.  Transform the 2nd Piola-kirchhoff stress in the covariant system
            NCovariant = TFromLCToCov*NLC;
            
            %% 4ii.13. Compute the first variation of the Green-Lagrange strain in the contravariant basis with respect to the DOFs
            dEpsilonContra = [A1(:,1)'*dRdxiMtx
                              A2(:,1)'*dRdetaMtx
                              .5*(A2(:,1)'*dRdxiMtx + A1(:,1)'*dRdetaMtx)];
                              
            %% 4ii.14. Transform the first variation of the Green-Lagrange strain at the local Cartesian basis
            dEpsilonCartesian = TFromContraToLC4VoigtStrain*dEpsilonContra;
            
            %% 4ii.15. Compute the first variation of the 2nd Piola-Kichhoff stress in the local Cartesian basis
            dNCartesian = Dm*dEpsilonCartesian;
            
            %% 4ii.16. Transform the first variation of the 2nd Piola-Kichhoff stress at the covariant basis
            dNCovariant = TFromLCToCov*dNCartesian;
            
            %% 4ii.17. Compute the first variation of the traction vector
            dtractionVct = uContravariant(1,1)*A1*dNCovariant(1,:) + ...
                uContravariant(2,1)*A2*dNCovariant(2,:) + ...
                (uContravariant(2,1)*A1 + uContravariant(1,1)*A2)*dNCovariant(3,:) + ...
                (uContravariant(1,1)*NCovariant(1,1) + uContravariant(2,1)*NCovariant(3,1))*dRdxiMtx + ...
                (uContravariant(2,1)*NCovariant(2,1) + uContravariant(1,1)*NCovariant(3,1))*dRdetaMtx;
            
            %% 4ii.18. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
            if isOnXi
                detJxxi = norm(A1(:,1));
            else
                detJxxi = norm(A2(:,1));
            end
            
            %% 4ii.19. Compute the element length at the GP
            elementLengthOnGP = detJxxi*detJxizeta*GW(j);
            
            %% 4ii.20. Compute the element Q-matrix contribution corresponding to the eigenvalue problem for the estimation of the stabilization parameter for the weak application of the Dirichlet boundary conditions with the Nitsche method and add it to the global matrix
            QMtx(EFT,EFT) = QMtx(EFT,EFT) + (dtractionVct'*dtractionVct)*elementLengthOnGP;
        end
    end
end

end
