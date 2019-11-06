function KWeakDBCLMMembrane = computeWeakDBCMtxLagrangeMultipliersIGAMembrane...
    (BSplinePatch,connections,noDOFs,propCoupling)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the complement of the stiffness matrix corresponding to an
% isogeometric membrane which accounts for the application of the Dirichlet
% boundary conditions weakly using the Penalty method.
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
%         connections : Dummy variable for this function
%              noDOFs : The total number of DOFs for the multipatch
%                       structure with the number of DOFs for the Lagrange 
%                       Multipliers 
%        propCoupling : Dummy variable for this function
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the conditions
% ->
%    1i. Get the extensions of the Dirichlet boundary where weak boundary conditions using the penalty method are to be applied
%
%   1ii. Get the discretization of the Lagrange Multipliers discretization for the given condition
%
%  1iii. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
%
%   1iv. Compute the merged knot vector from the patch and the Lagrange Multipliers fields on the interface
%
%    1v. Issue Gauss Point coordinates and weights
%
%   1vi. Loop over the elements on the parameter space where the weak boundary conditions are applied
%   ->
%        1vi.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%        1vi.2. Loop over the Gauss points
%        ->
%               1v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%             1vi.2ii. Compute the NURBS basis functions
%
%            1vi.2iii. Create the element freedom table for the patch and the Lagrange Multipliers discretizations
%
%             1vi.2iv. Compute the covariant base vectors
%
%              1vi.2v. Compute the B-operator matrix for the displacement and the Lagrange Multipliers field
%
%             1vi.2vi. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%            1vi.2vii. Compute the element length at the GP
%
%           1vi.2viii. Compute the element stiffness contribution matrix for the weak application of the Dirichlet boundary conditions and add it to the global matrix
%        <-
%   <-
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

% Define a tolerance for the length of a curve
tolCurveLength = 1e-4;

% Reassign the analysis arrays
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Get the DOF numbering
DOFNumbering = BSplinePatch.DOFNumbering;

% Number of Control Points in xi-,eta- directions
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of local DOFs
noNodesLoc = (p+1)*(q+1);
noDOFsLoc = 3*noNodesLoc;

% Initialize element freedom table
EFT = zeros(1,noDOFsLoc);

% Initialize auxiliary arrays
RMtx = zeros(3,noDOFsLoc);

% % Initialize the output arrays
KWeakDBCLMMembrane = zeros(noDOFs,noDOFs);

%% 1. Loop over all the conditions
for iCnd = 1:BSplinePatch.weakDBC.noCnd
    %% 1i. Get the extensions of the Dirichlet boundary where weak boundary conditions using the penalty method are to be applied
    xiExtension = BSplinePatch.weakDBC.xiExtension{iCnd};
    etaExtension = BSplinePatch.weakDBC.etaExtension{iCnd};
    isWeakDBCOverPoint = false;
    
    %% 1ii. Get the discretization of the Lagrange Multipliers discretization for the given condition
    pLambda = BSplinePatch.weakDBC.lambda{iCnd}.p;
    XiLambda = BSplinePatch.weakDBC.lambda{iCnd}.Xi;
    CPLambda = BSplinePatch.weakDBC.lambda{iCnd}.CP;
    isNURBSLambda = BSplinePatch.weakDBC.lambda{iCnd}.isNURBS;
    mLambda = length(XiLambda);
    nxiLambda = length(CPLambda(:,1));
    DOFNumberingCnd = BSplinePatch.weakDBC.lambda{iCnd}.DOFNumbering;
    EFTOfCnd = BSplinePatch.weakDBC.lambda{iCnd}.EFT;
    noCPsLambdaLoc = pLambda + 1;
    noDOFsLambdaLoc = 3*noCPsLambdaLoc;
    EFTInCnd = zeros(1,noDOFsLambdaLoc);
    BLambdaMtx = zeros(3,noDOFsLambdaLoc);
    
    %% 1iii. Get the running and the fixed parameters on the patch boundary where the weak boundary conditions are applied
    if etaExtension(1) == etaExtension(2) && xiExtension(1) ~= xiExtension(2)
        % Coupled region in xi-direction
        weakDBCRegionPatch = xiExtension;

        % Find the correct spans for the coupled region
        spanStart = findKnotSpan(weakDBCRegionPatch(1),Xi,nxi);
        spanEnd = findKnotSpan(weakDBCRegionPatch(2),Xi,nxi) + 1;

        % Find the corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = Xi(spanStart:spanEnd);

        % Fixed parameter on the parametric net
        eta = etaExtension(1);

        % Find the span where xiEta it lies in
        etaSpan = findKnotSpan(eta,Eta,neta);

        % Flag on whether the coupling line is over xi
        isOnXi = true;
    elseif xiExtension(1) == xiExtension(2) && etaExtension(1) ~= etaExtension(2)
        % Coupled region in eta-direction
        weakDBCRegionPatch = etaExtension;

        % Find the correct spans for the coupled region
        spanStart = findKnotSpan(weakDBCRegionPatch(1),Eta,neta);   
        spanEnd = findKnotSpan(weakDBCRegionPatch(2),Eta,neta) + 1;

        % Find the corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = Eta(spanStart:spanEnd);

        % Fixed parameter on the parametric net
        xi = xiExtension(1);

        % Find the span where uv it lies in
        xiSpan = findKnotSpan(xi,Xi,nxi);

        % Flag on whether the coupling line is over eta
        isOnXi = false;
    elseif etaExtension(1) == etaExtension(2) && xiExtension(1) == xiExtension(2)
        % Initialize flag on whether the application of weak Dirichlet
        % boundary conditions is over a point or over a line
        isWeakDBCOverPoint = true;
        
        % Fixed parameters on the parametric net
        xi = xiExtension(1);
        eta = etaExtension(1);
        
        % Find the correct spans for the point
        xiSpan = findKnotSpan(xi,Xi,nxi);
        etaSpan = findKnotSpan(eta,Eta,neta);
        
        % Find the corresponding to the coupled region knot span
        weakDBCRgionOnKnotVector = zeros(2,1);
    else
        error('Boundary over which weak Dirichlet boundary conditions are to be imposed wrongly defined');
    end

    %% 1iv. Compute the merged knot vector from the patch and the Lagrange Multipliers fields on the interface
    if ~isWeakDBCOverPoint
        % Projected knot vector for the Lagrange multipliers field on the
        % structural coupled field:
        % If UL = [v1 v2] and US = [u1 u2] the two knot vectors, then the
        % transformation rule is t(v) = u1 + (u1-u2)*(v-v1)/(v1-v2)

        % Initialization of the projected knot vector
        XiLambdaProjected = XiLambda;

        % Assign the new entries of the projected Lagrange multipliers knot vector
        for i = 1:length(XiLambda)
           XiLambdaProjected(i) = weakDBCRegionPatch(1) + (weakDBCRegionPatch(1)-weakDBCRegionPatch(2))*(XiLambda(i)-XiLambda(1))/(XiLambda(1)-XiLambda(length(XiLambda)));
        end

        % Weak Dirichlet boundary condition region on the Lagrange Multipliers
        % knot span
        weakDBCRegionLambda = XiLambdaProjected(pLambda+2:mLambda-pLambda-1);

        % Merged Dirichlet boundary condition region
        weakDBCRegion = mergesorted(weakDBCRegionPatch,weakDBCRegionLambda);
        weakDBCRgionOnKnotVector = unique(weakDBCRegion);
    end
    
    %% 1v. Issue Gauss Point coordinates and weights
    if ~isWeakDBCOverPoint
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
    else
        noGPs = 1;
    end

    %% 1vi. Loop over the elements on the parameter space where the weak boundary conditions are applied
    for i = 1:length(weakDBCRgionOnKnotVector)-1
        if (weakDBCRgionOnKnotVector(i) ~= weakDBCRgionOnKnotVector(i+1)) || isWeakDBCOverPoint
            %% 1vi.1. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
            if ~isWeakDBCOverPoint
                detJxizeta = (weakDBCRgionOnKnotVector(i+1)-weakDBCRgionOnKnotVector(i))/2;
            end

            %% 1vi.2. Loop over the Gauss points
            for j = 1:noGPs
                %% 1v.2i. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
                
                % For the patch :
                % _______________
                
                if ~isWeakDBCOverPoint
                    xiEta = ((1-GP(j))*weakDBCRgionOnKnotVector(i) + (1+GP(j))*weakDBCRgionOnKnotVector(i+1))/2;
                end
                
                % For the Lagrange Multipliers field :
                % ____________________________________
                
                xiLambda = (XiLambda(1)-XiLambda(length(XiLambda)))*(xiEta-weakDBCRgionOnKnotVector(1))/...
                    (weakDBCRgionOnKnotVector(1)-weakDBCRgionOnKnotVector(length(weakDBCRgionOnKnotVector)))+XiLambda(1);

                %% 1vi.2ii. Compute the NURBS basis functions
                
                % For patch :
                % ___________
                
                if ~isWeakDBCOverPoint
                    if isOnXi
                        xi = xiEta;
                        xiSpan = findKnotSpan(xi,Xi,nxi);
                    else
                        eta = xiEta;
                        etaSpan = findKnotSpan(eta,Eta,neta);
                    end
                end
                dR = computeIGABasisFunctionsAndDerivativesForSurface...
                    (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,1);
                
                % For the Lagrange Multipliers field :
                % ____________________________________
                
                % Find the span on the Lagrange multipliers knot span
                xiSpanLambda = findKnotSpan(xiLambda,XiLambda,nxiLambda);

                % Compute the NURBS basis functions
                RLambda = computeIGABasisFunctionsAndDerivativesForCurve...
                    (xiSpanLambda,pLambda,xiLambda,XiLambda,CPLambda,0,isNURBSLambda);

                %% 1vi.2iii. Create the element freedom table for the patch and the Lagrange Multipliers discretizations
                
                % For patch :
                % ___________

                % Initialize of the counter
                r = 1;

                % Relation global-local DoFs
                for cpj = etaSpan-q:etaSpan
                    for cpi = xiSpan-p:xiSpan
                        EFT(r)   = DOFNumbering(cpi,cpj,1);
                        EFT(r+1) = DOFNumbering(cpi,cpj,2);
                        EFT(r+2) = DOFNumbering(cpi,cpj,3);

                        % update counter
                        r = r + 3;
                    end
                end

                % For the Lagrange Multipliers field :
                % ____________________________________
                
                % Initialize of the counter
                r = 1;

                % Relation global-local DoFs
                for cpi = xiSpanLambda-pLambda:xiSpanLambda
                    EFTInCnd(r) = DOFNumberingCnd(cpi,1);
                    EFTInCnd(r+1) = DOFNumberingCnd(cpi,2);
                    EFTInCnd(r+2) = DOFNumberingCnd(cpi,3);

                    % Update counter
                    r = r + 3;
                end

                %% 1vi.2iv. Compute the covariant base vectors
                [A1,A2] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CP,0,dR);

                %% 1vi.2v. Compute the B-operator matrix for the displacement and the Lagrange Multipliers field

                % For patch :
                % ___________
                    
                % initialize counter
                k = 0;

                % Loop over all the non-zero contributions at the span
                % under study
                for c = 0:q
                    for b = 0:p
                        % Update counter
                        k = k + 1;

                        % Matrix containing the basis functions
                        RMtx(1,3*k-2) = dR(k,1);
                        RMtx(2,3*k-1) = dR(k,1);
                        RMtx(3,3*k) = dR(k,1);
                    end
                end
                
                % For the Lagrange Multipliers field :
                % ____________________________________
                
                for k = 1:noCPsLambdaLoc
                    BLambdaMtx(1,3*k-2) = RLambda(k);
                    BLambdaMtx(2,3*k-1) = RLambda(k);
                    BLambdaMtx(3,3*k) = RLambda(k);
                end

                %% 1vi.2vi. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
                if ~isWeakDBCOverPoint
                    if isOnXi
                        detJxxi = norm(A1(:,1));
                    else
                        detJxxi = norm(A2(:,1));
                    end
                end

                %% 1vi.2vii. Compute the element length at the GP
                if ~isWeakDBCOverPoint
                    elementLengthOnGP = detJxxi*detJxizeta*GW(j);
                else
                    elementLengthOnGP = 1;
                end

                %% 1vi.2viii. Compute the element stiffness contribution matrix for the weak application of the Dirichlet boundary conditions and add it to the global matrix
                if BSplinePatch.weakDBC.alpha ~= 0
                    KWeakDBCLMMembrane(EFT,EFT) = KWeakDBCLMMembrane(EFT,EFT) + ...
                        BSplinePatch.weakDBC.alpha*(RMtx'*RMtx)*elementLengthOnGP;
                end
                KWeakDBCLMMembrane(EFT,EFTOfCnd(EFTInCnd)) = KWeakDBCLMMembrane(EFT,EFTOfCnd(EFTInCnd)) + ...
                    RMtx'*BLambdaMtx*elementLengthOnGP;
                KWeakDBCLMMembrane(EFTOfCnd(EFTInCnd),EFT) = KWeakDBCLMMembrane(EFTOfCnd(EFTInCnd),EFT) + ...
                    BLambdaMtx'*RMtx*elementLengthOnGP;
            end
        end
    end
end



end
