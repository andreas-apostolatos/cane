function [tangMtxWeakDBCNitsche, resVctWeakDBCNitsche] = ...
    computeWeakDBCTangMtxResVctNitscheIGAMembraneOutdatedKiendl ...
    (BSplinePatch, dHat, connections, numDOFs, propCoupling)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangent stiffness matrix and the residual vector for a single
% B-Spline patch correspondingo to the application of the Nitsche method
% for enfocing the Dirichlet boundary conditions weakly.
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
%                dHat : The displacement solution vector from the previous
%                       nonlinear iteration step
%         connections : Dummy variable for this function
%             numDOFs : Dummy variable for this function
%        propCoupling : Dummy variable for this function
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the conditions
% ->
%    1i. Get the extensions of the Dirichlet boundary where weak boundary conditions using the Nitsche method are to be applied
%
%   1ii. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
%
%  1iii. Get the parameter space of the application of the weak boundary conditions using the knot vector information
%
%   1iv. Issue Gauss Point coordinates and weights
%
%    1v. Loop over the elements on the parameter space where the weak boundary conditions are applied
%    ->
%        1v.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        1v.2. Loop over the Gauss points
%        ->
%              1v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%             1v.2ii. Compute the NURBS basis functions
%
%            1v.2iii. Create the element freedom table
%
%             1v.2iv. Get the displacement vector of the previous iteration step at the Gauss point
%
%              1v.2v. Compute the matrices of the basis functions and their derivatives
%
%             1v.2vi. Compute the covariant base vectors of the reference configuration
%
%            1v.2vii. Compute the normal to the boundary vector and transform it into the contavariant basis
%
%           1v.2viii. Compute the covariant base vectors of the current configuration
%
%             1v.2ix. Compute the covariant metric coefficients of the current configuration
%
%              1v.2x. Compute the contravariant basis
%
%             1v.2xi. Compute the local Cartesian basis
%
%            1v.2xii. Compute the transformation matrix from the contravariant basis to the local Cartesian one
%
%           1v.2xiii. Compute the transformation matrix from the local Cartesian basis to the covariant one
%
%            1v.2xiv. Compute the Green-Lagrange strain in the contravariant basis
%
%             1v.2xv. Transform the Green-Lagrange strain at the local Cartesian basis
%
%            1v.2xvi. Compute the 2nd Piola-Kirchhoff stress in the local Cartesian system
%
%           1v.2xvii. Transform the 2nd Piola-Kirchhoff stress at the covariant basis
%
%          1v.2xviii. Arrange the covariant 2nd Piola-Kirchhoff stress tensor in a matrix
%
%            1v.2xix. Compute the traction vector at the Dirichlet boundary
%
%             1v.2xx. Loop over the r-DOFs and compute the first variations
%             ->
%                     1v.2xx.1. Get local node number kr and dof direction dirr
%
%                     1v.2xx.2. Compute the variation of the covariant base vectors wrt the DOFs
%
%                     1v.2xx.3. Compute first variation of the Green-Lagrange strain in the contravariant system wrt the DOFS
%
%                     1v.2xx.4. Transform the variation of the Green-Lagrange strain at the local Cartesian system
%
%                     1v.2xx.5. Compute the variation of the 2nd Piola-Kirchhoff stress in the local Cartesian basis
%
%                     1v.2xx.6. Transform the variation of the 2nd Piola-Kirchhoff stress at the covariant basis
%
%                     1v.2xx.7. Arrange the variation of the 2nd Piola-Kirchhoff stress tensor at the covariant basis in a matrix
%
%                     1v.2xx.8. Compute the first variation of the traction vector wrt the DOFs
%             <-
%
%            1v.2xxi. Loop over the r-DOFs and compute the first variations
%            ->
%                     1v.2xxi.1. Get local node number kr and dof direction dirr
%
%                     1v.2xxi.2. Compute the variation of the covariant base vectors wrt the DOFs
%
%                     1v.2xxi.3. Loop over the s-DOFs and compute the second variations
%                     ->
%                                1v.2xxi.3i. Get local node number ks and dof direction dirs for DoF us
%
%                               1v.2xxi.3ii. Compute the variation of the covariant base vectors wrt the DOFs
%
%                              1v.2xxi.3iii. Compute second variation of the Green-Lagrange strain in the contravariant system wrt DOFs
%
%                               1v.2xxi.3iv. Transform the variation of the Green-Lagrange strain at the local Cartesian system
%
%                                1v.2xxi.3v. Compute the second variation of the 2nd Piola-Kirchhoff stress in the local Cartesian basis
%
%                               1v.2xxi.3vi. Transform the second variation of the 2nd Piola-Kirchhoff stress at the covariant basis
%
%                              1v.2xxi.3vii. Arrange the second variation of the 2nd Piola-Kirchhoff stress at the covariant basis in a matrix
%
%                             1v.2xxi.3viii. Compute the second variation of the traction vector
%                     <-
%            <-
%
%           1v.2xxii. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%          1v.2xxiii. Compute the element length at the GP
%
%           1v.2xxiv. Compute the element tangent matrix and add it to the master tangent stiffness matrix corresponding to the Nitsche method for enforcing the Dirichlet boundary conditions weakly
%
%            1v.2xxv. Compute the element residual vector and ass it to the global residual vector corresponding to the Nitsche method for enforcing the Dirichlet boundary conditions weakly
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
numDOFs = 3*nxi*neta;

% Initialize the element freedom table
EFT = zeros(1,noDOFsEl);

% Initialize auxiliary arrays
RMtx = zeros(3,noDOFsEl);
dNalphaBeta = zeros(noDOFsEl,2,2);
dtractionVct = zeros(3,noDOFsEl);
ddtractionVct = zeros(3,noDOFsEl,noDOFsEl);

% Initialize the output arrays
tangMtxWeakDBCNitsche = zeros(numDOFs,numDOFs);
resVctWeakDBCNitsche = zeros(numDOFs,1);

%% 1. Loop over all the conditions
for counterCnd = 1:BSplinePatch.weakDBC.noCnd
    %% 1i. Get the extensions of the Dirichlet boundary where weak boundary conditions using the Nitsche method are to be applied
    xiExtension = BSplinePatch.weakDBC.xiExtension{counterCnd};
    etaExtension = BSplinePatch.weakDBC.etaExtension{counterCnd};
    
    %% 1ii. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
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

    %% 1iii. Get the parameter space of the application of the weak boundary conditions using the knot vector information
    weakDBCRgionOnKnotVector = unique(weakDBCRgionOnKnotVector);

    %% 1iv. Issue Gauss Point coordinates and weights
    if strcmp(BSplinePatch.weakDBC.int.type,'default')
        if isOnXi
            pDegree = p + 1;
        else
            pDegree = q + 1;
        end
        noGPs = ceil((pDegree + 1)/2);
    elseif strcmp(BSplinePatch.weakDBC.int,'user')
        noGPs = BSplinePatch.weakDBC.noGPs;
    end
    [GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);

    %% 1v. Loop over the elements on the parameter space where the weak boundary conditions are applied
    for i = 1:length(weakDBCRgionOnKnotVector)-1
        if weakDBCRgionOnKnotVector(i) ~= weakDBCRgionOnKnotVector(i+1)
            %% 1v.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
            detJxizeta = (weakDBCRgionOnKnotVector(i+1)-weakDBCRgionOnKnotVector(i))/2;

            %% 1v.2. Loop over the Gauss points
            for j = 1:noGPs
                %% 1v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                xiEta = ((1-GP(j))*weakDBCRgionOnKnotVector(i) + (1+GP(j))*weakDBCRgionOnKnotVector(i+1))/2;

                %% 1v.2ii. Compute the NURBS basis functions
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
                
                %% 1v.2iv. Get the displacement vector of the previous iteration step at the Gauss point
                dHatEl = dHat(EFT);
                dispVct = computePostprocDisplacementIGAKirchhoffLoveShell...
                    (p,q,dR(:,1),dHatEl);
                
                %% 1v.2v. Compute the matrices of the basis functions and their derivatives
                for iCPs = 1:noCPsEl
                    RMtx(1,3*iCPs - 2) = dR(iCPs,1);
                    RMtx(2,3*iCPs - 1) = dR(iCPs,1);
                    RMtx(3,3*iCPs) = dR(iCPs,1);
                end
                
                %% 1v.2vi. Compute the covariant base vectors of the reference configuration
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CP,0,dR);
                A3 = cross(A1,A2)/norm(cross(A1,A2));
                
                %% 1v.2vii. Compute the normal to the boundary vector and transform it into the contavariant basis
                [uGC,~] = computeNormalAndTangentVectorsToBSplineBoundary...
                    (xi,Xi,eta,Eta,A1,A2,A3,isOnXi);
                uContravariant = [A1'
                                  A2']*uGC;
                              
                %% 1v.2viii. Compute the covariant base vectors of the current configuration
                [a1,a2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CPd,0,dR);
                
                %% 1v.2ix. Compute the covariant metric coefficients of the current configuration
                aabCovariant = [a1 a2]'*[a1 a2];
                
                %% 1v.2x. Compute the contravariant basis
                AabCovariant = [A1 A2]'*[A1 A2];
                AContravariant = AabCovariant\[A1 A2]';
                AContravariant = AContravariant';
                
                %% 1v.2xi. Compute the local Cartesian basis
                eLC = computeLocalCartesianBasis4BSplineSurface...
                    ([A1 A2],AContravariant);
                
                %% 1v.2xii. Compute the transformation matrix from the contravariant basis to the local Cartesian one
                TFromContraToLC4VoigtStrain = ...
                    computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell...
                    (eLC,AContravariant);
                
                %% 1v.2xiii. Compute the transformation matrix from the local Cartesian basis to the covariant one
                TFromLCToCov = computeTFromLocalCartesian2CovariantBasis4BSplineSurface...
                    (eLC,AContravariant);
                
                %% 1v.2xiv. Compute the Green-Lagrange strain in the contravariant basis
                EpsilonContra = .5*[aabCovariant(1,1) - AabCovariant(1,1)
                                    aabCovariant(2,2) - AabCovariant(2,2)
                                    aabCovariant(1,2) - AabCovariant(1,2)];
                                
                %% 1v.2xv. Transform the Green-Lagrange strain at the local Cartesian basis
                EpsilonLC = TFromContraToLC4VoigtStrain*EpsilonContra;
                
                %% 1v.2xvi. Compute the 2nd Piola-Kirchhoff stress in the local Cartesian system
                NLC = prestress +  Dm*EpsilonLC;
                
                %% 1v.2xvii. Transform the 2nd Piola-Kirchhoff stress at the covariant basis
                NCovariant = TFromLCToCov*NLC;
                
                %% 1v.2xviii. Arrange the covariant 2nd Piola-Kirchhoff stress tensor in a matrix
                PalphaBeta = [NCovariant(1,1) NCovariant(3,1)
                              NCovariant(3,1) NCovariant(2,1)];
                
                %% 1v.2xix. Compute the traction vector at the Dirichlet boundary
                tractionVct = [a1 a2]*PalphaBeta*uContravariant;
                
                %% 1v.2xx. Loop over the r-DOFs and compute the first variations
                for rDOFs = 1:noDOFsEl
                    %% 1v.2xx.1. Get local node number kr and dof direction dirr
                    kr = ceil(rDOFs/3);
                    dirr = rDOFs-3*(kr-1);
                    
                    %% 1v.2xx.2. Compute the variation of the covariant base vectors wrt the DOFs
                    daCovariant_r = zeros(2,3);
                    daCovariant_r(1,dirr) = dR(kr,2);
                    daCovariant_r(2,dirr) = dR(kr,3);
                    
                    %% 1v.2xx.3. Compute first variation of the Green-Lagrange strain in the contravariant system wrt the DOFS
                    dEpsilonContra = [dR(kr,2)*a1(dirr,1)
                                      dR(kr,3)*a2(dirr,1)
                                      .5*(dR(kr,2)*a2(dirr,1) + a1(dirr,1)*dR(kr,3))];

                    %% 1v.2xx.4. Transform the variation of the Green-Lagrange strain at the local Cartesian system
                    dEpsilonLC = TFromContraToLC4VoigtStrain*dEpsilonContra;
                    
                    %% 1v.2xx.5. Compute the variation of the 2nd Piola-Kirchhoff stress in the local Cartesian basis
                    dNLC = Dm*dEpsilonLC;
                    
                    %% 1v.2xx.6. Transform the variation of the 2nd Piola-Kirchhoff stress at the covariant basis
                    dNCovariant = TFromLCToCov*dNLC;
                    
                     %% 1v.2xx.7. Arrange the variation of the 2nd Piola-Kirchhoff stress tensor at the covariant basis in a matrix
                    dNalphaBeta(rDOFs,:,:) = [dNCovariant(1,1) dNCovariant(3,1)
                                              dNCovariant(3,1) dNCovariant(2,1)];

                    %% 1v.2xx.8. Compute the first variation of the traction vector wrt the DOFs
                    dtractionVct(:,rDOFs) = [a1 a2]*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant + ...
                        daCovariant_r'*PalphaBeta*uContravariant;
                end

                %% 1v.2xxi. Loop over the r-DOFs and compute the first variations
                for rDOFs = 1:noDOFsEl
                     %% 1v.2xxi.1. Get local node number kr and dof direction dirr
                    kr = ceil(rDOFs/3);
                    dirr = rDOFs - 3*(kr-1);
                    
                    %% 1v.2xxi.2. Compute the variation of the covariant base vectors wrt the DOFs
                    daCovariant_r = zeros(2,3);
                    daCovariant_r(1,dirr) = dR(kr,2);
                    daCovariant_r(2,dirr) = dR(kr,3);
                    
                    %% 1v.2xxi.3. Loop over the s-DOFs and compute the second variations
                    for sDOFs = 1:noDOFsEl
                        %% 1v.2xxi.3i. Get local node number ks and dof direction dirs for DoF us
                        ks = ceil(sDOFs/3);
                        dirs = sDOFs-3*(ks-1);
                        
                        %% 1v.2xxi.3ii. Compute the variation of the covariant base vectors wrt the DOFs
                        daCovariant_s = zeros(2,3);
                        daCovariant_s(1,dirs) = dR(ks,2);
                        daCovariant_s(2,dirs) = dR(ks,3);
                        
                        %% 1v.2xxi.3iii. Compute second variation of the Green-Lagrange strain in the contravariant system wrt DOFs
                        if dirr == dirs
                            ddEpsilonContra = [dR(kr,2)*dR(ks,2)
                                               dR(kr,3)*dR(ks,3)
                                               .5*(dR(kr,2)*dR(ks,3) + dR(kr,3)*dR(ks,2))];
                        else
                            ddEpsilonContra = zeros(3,1);
                        end
                        
                        %% 1v.2xxi.3iv. Transform the variation of the Green-Lagrange strain at the local Cartesian system
                        ddEpsilonLC = TFromContraToLC4VoigtStrain*ddEpsilonContra;
                        
                       %% 1v.2xxi.3v. Compute the second variation of the 2nd Piola-Kirchhoff stress in the local Cartesian basis
                       ddNLC = Dm*ddEpsilonLC;
                       
                       %% 1v.2xxi.3vi. Transform the second variation of the 2nd Piola-Kirchhoff stress at the covariant basis
                       ddNCovariant = TFromLCToCov*ddNLC;
                       
                       %% 1v.2xxi.3vii. Arrange the second variation of the 2nd Piola-Kirchhoff stress at the covariant basis in a matrix
                       ddNalphaBeta = [ddNCovariant(1,1) ddNCovariant(3,1)
                                       ddNCovariant(3,1) ddNCovariant(2,1)];

                       %% 1v.2xxi.3viii. Compute the second variation of the traction vector
                       ddtractionVct(:,rDOFs,sDOFs) = [a1 a2]*ddNalphaBeta*uContravariant + ...
                            daCovariant_s'*squeeze(dNalphaBeta(rDOFs,:,:))*uContravariant + ...
                            daCovariant_r'*squeeze(dNalphaBeta(sDOFs,:,:))*uContravariant;
                    end
                end
                
                %% 1v.2xxii. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                if isOnXi
                    detJxxi = norm(A1(:,1));
                else
                    detJxxi = norm(A2(:,1));
                end

                %% 1v.2xxiii. Compute the element length at the GP
                elementLengthOnGP = detJxxi*detJxizeta*GW(j);
                
                %% 1v.2xxiv. Compute the element tangent matrix and add it to the master tangent stiffness matrix corresponding to the Nitsche method for enforcing the Dirichlet boundary conditions weakly
                tangMtxWeakDBCNitsche(EFT,EFT) = tangMtxWeakDBCNitsche(EFT,EFT) - ...
                    (dtractionVct'*RMtx + RMtx'*dtractionVct + ...
                    squeeze(dispVct(1,1)*ddtractionVct(1,:,:) + ...
                    dispVct(2,1)*ddtractionVct(2,:,:) + ...
                    dispVct(3,1)*ddtractionVct(3,:,:)))*elementLengthOnGP;
                
                %% 1v.2xxv. Compute the element residual vector and ass it to the global residual vector corresponding to the Nitsche method for enforcing the Dirichlet boundary conditions weakly
                resVctWeakDBCNitsche(EFT,1) = resVctWeakDBCNitsche(EFT,1) - ...
                    (dtractionVct'*dispVct + RMtx'*tractionVct)*elementLengthOnGP;
            end
        end
    end
end

end
