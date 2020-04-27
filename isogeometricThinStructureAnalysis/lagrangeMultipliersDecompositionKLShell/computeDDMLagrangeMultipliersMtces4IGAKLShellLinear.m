function [Lambda, Mu] = ...
    computeDDMLagrangeMultipliersMtces4IGAKLShellLinear ...
    (BSplinePatch, lambda, mu, isMaster, isSameOrientation, propCoupling)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the Lagrange mutlipliers matrices for the multipatch coupling of
% the Kirchhoff-Love shell problem using the Lagrange Multipliers method.
%
%               Input :
%        BSplinePatch : The B-Spline patch containing polynomial orders, 
%                       knot vectors, Control Point coordinates etc.
%              lambda : B-Spline discretization of the Lagrange Multipliers 
%                       field for the traction forces
%                  mu : B-Spline discretization of the Lagrange Multipliers 
%                       field for the traction moments
%            isMaster : Flag on whether the patch is the master or the 
%                       slave one
%   isSameOrientation : Flag on whether the interfaces from both patches
%                       are oriented in the same direction
%        propCoupling : Properties of the multipatch coupling
%                           .alphaD : penalty factor for the displacement
%                                     coupling
%                           .alphaR : penalty factor for the rotation
%                                     coupling
%                             .intC : On the integration of the coupling
%                                     interface 
%
%              Output :
%              Lambda : The Lagrange Multipliers matrix for the traction
%                       forces
%                  Mu : The Lagrange Multipliers matrix for the traction
%                       moments
%
% Function layout :
%
% 0. Read input
%
% 1. Get the running and the fixed parameters on the patch interface and the coupling region
%
% 2. Compute the merged knot vector from the patch and the Lagrange Multipliers fields on the interface
%
% 3. Issue Gauss Point coordinates and weights
%
% 4. Loop over all the elements on the coupling interface
% ->
%    4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%  4ii. Loop over all Gauss points
%  ->
%       4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%       4ii.2. Compute the NURBS basis functions and their derivatives
%
%       4ii.3. Create the element freedom tables for the patch and the Lagrange Multipliers fields
%
%       4ii.4. Compute the covariant base vectors
%
%       4ii.5. Compute the surface normal vectors
%
%       4ii.6. Compute the derivatives of the surface normal vector
%
%       4ii.7. Compute the normal to the boundary vector
%
%       4ii.8. Compute the covariant metric coefficients
%
%       4ii.9. Compute the contravariant base vectors
%
%       4ii.10. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%       4ii.11. Transform the normal and the tangent vectors to the covariant bases
%
%       4ii.12. Compute the curvature coefficients
%
%       4ii.13. Compute the B-operator matrices for the rotations
%
%       4ii.14. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%       4ii.15. Compute the element length at the GP
%
%       4ii.16. Compute the Lagrange Multipliers matrices and add them to the global matrices
%  <-
% <-
%
%% Function main body

%% 0. Read input

% Re-assign the patch properties
p = BSplinePatch.p;
q = BSplinePatch.q;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
xicoup = BSplinePatch.xicoup;
etacoup = BSplinePatch.etacoup;
noDOFs = BSplinePatch.noDOFs;

% Number of Control Points in xi-,eta- directions
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Local number of control points that influence the element for the
% membrane surface
noCPsLoc = (p+1)*(q+1);

% Local number of degrees of freedom for the membrane surface
noDOFsLoc = 3*noCPsLoc;

% Check if Lagrange Multipliers field for the rotational coupling is
% enforced
isLagrangeMultipliersRotationsEnabled = true;
if ischar(mu)
    isLagrangeMultipliersRotationsEnabled = false;
end

% Control Points for the Lagrange multipliers field
CPLambda = lambda.CP;
if isLagrangeMultipliersRotationsEnabled
    CPMu = mu.CP;
end

% Flag on whether the basis of the Lagrange Multipliers fields is a NURBS
% or a B-Spline
isNURBSLambda = lambda.isNURBS;
if isLagrangeMultipliersRotationsEnabled
    isNURBSMu = mu.isNURBS;
end

% Number of degrees of freedom for the Lagrange multipliers field (note we
% are given only the weights since not actual NURBS curve is needed)
nxiLambda = length(CPLambda(:,1));
if isLagrangeMultipliersRotationsEnabled
    nxiMu = length(CPMu(:,1));
end

% Polynomial degree for the Lagrange multipliers field
pLambda = lambda.p;
if isLagrangeMultipliersRotationsEnabled
    pMu = mu.p;
end

% Knot vector and its length, for the Lagrange multipliers field
XiLambda = lambda.Xi;
mLambda = length(XiLambda);
if isLagrangeMultipliersRotationsEnabled
    XiMu = mu.Xi;
    mMu = length(XiMu);
end

% Check if the given paramters for the Lagrange multipliers field are
% compatible
checkInputForBSplineCurve(pLambda,mLambda,nxiLambda);
if isLagrangeMultipliersRotationsEnabled
    checkInputForBSplineCurve(pMu,mMu,nxiMu);
end

% Number of nodes (Control Points) and degrees of freedom for the Lagrange
% multipliers field
noDOFsLambda = 3*nxiLambda;
if isLagrangeMultipliersRotationsEnabled
    noDOFsMu = 2*nxiMu;
end

% Local number of control points that affect the coupling surface
noCPsLambdaLoc = pLambda + 1;
noDOFsLambdaLoc = 3*noCPsLambdaLoc;
if isLagrangeMultipliersRotationsEnabled
    noCPsMuLoc = pMu + 1;
    noDOFsMuLoc = 2*noCPsMuLoc;
end

% Get the DOF numbering for the patch as well as the Lagrange Multipliers
% fields
DOFNumbering = BSplinePatch.DOFNumbering;
DOFNumberingLambda = lambda.DOFNumbering;
if isLagrangeMultipliersRotationsEnabled
    DOFNumberingMu = mu.DOFNumbering;
end

% Initialize the element freedom tables
EFT = zeros(1,noDOFsLoc);
EFTLambda = zeros(1,noDOFsLambdaLoc);
if isLagrangeMultipliersRotationsEnabled
    EFTMu = zeros(1,noDOFsMuLoc);
end

% Initialize auxiliary arrays
RMatrix = zeros(3,noDOFsLoc);
if isLagrangeMultipliersRotationsEnabled
    DRsDxiMatrix = zeros(3,noDOFsLoc);
    DRsDetaMatrix = zeros(3,noDOFsLoc);
    RMuMatrix = zeros(2,noCPsMuLoc);
end
RLambdaMatrix = zeros(3,noCPsLambdaLoc);

% Initialize the output arrays
Lambda = zeros(noDOFs,noDOFsLambda);
if isLagrangeMultipliersRotationsEnabled
    Mu = zeros(noDOFs,noDOFsMu);
else
    Mu = 'undefined';
end

% Define a singularity tolerance
tolArea = 1e-8;

%% 1. Get the running and the fixed parameters on the patch interface and the coupling region
if etacoup(1) == etacoup(2)
    % Coupled region in xi-direction
    couplingRegion = xicoup;
    
    % Find the correct spans for the coupled region
    spanStart = findKnotSpan(couplingRegion(1),Xi,nxi);
    spanEnd = findKnotSpan(couplingRegion(2),Xi,nxi)+1;
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorPatch = Xi(spanStart:spanEnd);
    
    % Fixed parameter on the parametric net
    eta = etacoup(1);
    
    % Find the span where xiEta it lies in
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    % Flag on whether the coupling line is over xi
    isOnXi = true;
else
    % Coupled region in eta-direction
    couplingRegion = etacoup;
    
    % Find the correct spans for the coupled region
    spanStart = findKnotSpan(couplingRegion(1),Eta,neta);   
    spanEnd = findKnotSpan(couplingRegion(2),Eta,neta)+1;   
    
    % Corresponding to the coupled region knot span
    couplingRegionOnKnotVectorPatch = Eta(spanStart:spanEnd);
    
    % Fixed parameter on the parametric net
    xi = xicoup(1);
    
    % Find the span where xiEta it lies in
    xiSpan = findKnotSpan(xi,Xi,nxi);
    
    % Flag on whether the coupling line is over eta
    isOnXi = false;
end

%% 2. Compute the merged knot vector from the patch and the Lagrange Multipliers fields on the interface

% Projected knot vector for the Lagrange multipliers field on the
% structural coupled field:
% If UL = [v1 v2] and US = [u1 u2] the two knot vectors, then the
% transformation rule is t(v) = u1 + (u1-u2)*(v-v1)/(v1-v2)

% Lagrange multipliers field for the force tractions :
% ____________________________________________________

% Initialization of the projected knot vector
XiLambdaProjected = XiLambda;

% Assign the new entries of the projected Lagrange multipliers knot vector
for i = 1:length(XiLambda)
   XiLambdaProjected(i) = couplingRegion(1) + (couplingRegion(1)-couplingRegion(2))*(XiLambda(i)-XiLambda(1))/(XiLambda(1)-XiLambda(length(XiLambda)));
end

% Coupled parameter knot vector for the Lagrange multipliers field
coupledRegionLambda = XiLambdaProjected(pLambda+2:mLambda-pLambda-1);

% Lagrange multipliers field for the moments tractions :
% ______________________________________________________

if isLagrangeMultipliersRotationsEnabled
    % Initialization of the projected knot vector
    XiMuProjected = XiMu;

    % Assign the new entries of the projected Lagrange multipliers knot vector
    for i = 1:length(XiMu)
       XiMuProjected(i) = couplingRegion(1) + (couplingRegion(1)-couplingRegion(2))*(XiMu(i)-XiMu(1))/(XiMu(1)-XiMu(length(XiMu)));
    end

    % Coupled parameter knot vector for the Lagrange multipliers field
    coupledRegionMu = XiMuProjected(pMu+2:mMu-pMu-1);
end

% Coupled parameter knot vector for the integration
couplingRegionOnKnotVectorTemp = mergesorted(couplingRegionOnKnotVectorPatch,coupledRegionLambda);
if isLagrangeMultipliersRotationsEnabled
    couplingRegionOnKnotVectorTemp = mergesorted(couplingRegionOnKnotVectorTemp,coupledRegionMu);
end
couplingRegionOnKnotVector = unique(couplingRegionOnKnotVectorTemp);

%% 3. Issue Gauss Point coordinates and weights
if strcmp(propCoupling.intC.type,'default')
    if isOnXi
        pDegree = p + 1;
    else
        pDegree = q + 1;
    end
    pLM = pLambda;
    if isLagrangeMultipliersRotationsEnabled
        pLM = max(pLM,pMu);
    end
    noGPs = ceil((pDegree + pLM + 1)/2);
elseif strcmp(propCoupling.intC.type,'user')
    noGPs = propCoupling.intC.noGPs;
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGPs);

%% 4. Loop over all the elements on the coupling interface
for i = 1:length(couplingRegionOnKnotVector)-1
    %% 4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
    detJxizeta = (couplingRegionOnKnotVector(i+1) - couplingRegionOnKnotVector(i))/2;
    
    %% 4ii. Loop over all Gauss points
    for j = 1:length(GP)
        %% 4ii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
        
        % For patch :
        % ___________
        
        xiEta = ((1-GP(j))*couplingRegionOnKnotVector(i)+(1+GP(j))*couplingRegionOnKnotVector(i+1))/2;
        
        % For the Lagrange Multipliers field for the traction forces :
        % ____________________________________________________________
        
        xiLambda = (XiLambda(1)-XiLambda(length(XiLambda)))*(xiEta-couplingRegionOnKnotVector(1))/...
            (couplingRegionOnKnotVector(1)-couplingRegionOnKnotVector(length(couplingRegionOnKnotVector)))+XiLambda(1);
        
        % For the Lagrange Multipliers field for the traction moments :
        % _____________________________________________________________
        
        if isLagrangeMultipliersRotationsEnabled
            xiMu = (XiMu(1)-XiMu(length(XiMu)))*(xiEta-couplingRegionOnKnotVector(1))/...
                (couplingRegionOnKnotVector(1)-couplingRegionOnKnotVector(length(couplingRegionOnKnotVector)))+XiMu(1);
        end
        
        %% 4ii.2. Compute the NURBS basis functions and their derivatives
        
        % For patch :
        % ___________
        
        if ~isMaster && ~isSameOrientation
            if isOnXi
                xi = xiEta;
                xiSpan = findKnotSpan(xi,Xi,nxi);
            else
                eta = xiEta;
                etaSpan = findKnotSpan(eta,Eta,neta);
            end
        else
            if isOnXi
                xi = Xi(length(Xi)) - xiEta;
                xiSpan = findKnotSpan(xi,Xi,nxi);
            else
                eta = Eta(length(Eta)) - xiEta;
                etaSpan = findKnotSpan(eta,Eta,neta);
            end
        end
        if isLagrangeMultipliersRotationsEnabled
            noDrvsBasis = 2;
        else
            noDrvsBasis = 1;
        end
        dR = computeIGABasisFunctionsAndDerivativesForSurface...
            (xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,noDrvsBasis);
       
        % For the Lagrange Multipliers field for the traction forces :
        % ____________________________________________________________
        
        % Find the span on the Lagrange multipliers knot span
        xiSpanLambda = findKnotSpan(xiLambda,XiLambda,nxiLambda);
        
        % Compute the NURBS basis functions
        RLambda = computeIGABasisFunctionsAndDerivativesForCurve...
            (xiSpanLambda,pLambda,xiLambda,XiLambda,CPLambda,0,isNURBSLambda);
        
        % For the Lagrange Multipliers field for the traction moments :
        % _____________________________________________________________
        
        if isLagrangeMultipliersRotationsEnabled
            % Find the span on the Lagrange multipliers knot span
            xiSpanMu = findKnotSpan(xiMu,XiMu,nxiMu);

            % Compute the NURBS basis functions
            RMu = computeIGABasisFunctionsAndDerivativesForCurve...
                (xiSpanMu,pMu,xiMu,XiMu,CPMu,0,isNURBSMu);
        end
        
        %% 4ii.3. Create the element freedom tables for the patch and the Lagrange Multipliers fields
        
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
        
        % For the Lagrange Multipliers field for the traction forces :
        % ____________________________________________________________

        % Initialize of the counter
        r = 1;

        % Relation global-local DoFs
        for cpi = xiSpanLambda-pLambda:xiSpanLambda
            EFTLambda(r) = DOFNumberingLambda(cpi,1);
            EFTLambda(r+1) = DOFNumberingLambda(cpi,2);
            EFTLambda(r+2) = DOFNumberingLambda(cpi,3);

            % Update counter
            r = r + 3;
        end
        
        % For the Lagrange Multipliers field for the traction moments :
        % ____________________________________________________________
        
        if isLagrangeMultipliersRotationsEnabled
            % Initialize of the counter
            r = 1;

            % Relation global-local DoFs
            for cpi = xiSpanMu-pMu:xiSpanMu
                EFTMu(r) = DOFNumberingMu(cpi,1);
                EFTMu(r+1) = DOFNumberingMu(cpi,2);

                % Update counter
                r = r + 2;
            end
        end
        
        %% 4ii.4. Compute the covariant base vectors
        if isLagrangeMultipliersRotationsEnabled
            noDrvsBaseVct = 1;
        else
            noDrvsBaseVct = 0;
        end
        [dA1,dA2] = computeBaseVectorsAndDerivativesForBSplineSurface...
            (xiSpan,p,etaSpan,q,CP,noDrvsBaseVct,dR);
        if norm(cross(dA1(:,1),dA2(:,1))) < tolArea
            continue;
        end

        %% 4ii.5. Compute the surface normal vectors
        if isLagrangeMultipliersRotationsEnabled
            A3Tilde = cross(dA1(:,1),dA2(:,1));
            A3 = A3Tilde/norm(A3Tilde);
        end
        
        %% 4ii.6. Compute the derivatives of the surface normal vector
        if isLagrangeMultipliersRotationsEnabled
            dA3 = computeParametricDrvsSurfaceNormalOnBSplineSurface...
                ([dA1(:,1) dA2(:,1)],[dA1(:,2) dA2(:,2) dA1(:,3)],...
                A3,norm(A3Tilde));
        end

        %% 4ii.7. Compute the normal to the boundary vector
        if isLagrangeMultipliersRotationsEnabled
            [n,t] = computeNormalAndTangentVectorsToBSplineBoundary...
                (xi,Xi,eta,Eta,dA1(:,1),dA2(:,1),A3,isOnXi);
        end
        
        %% 4ii.8. Compute the covariant metric coefficients
        if isLagrangeMultipliersRotationsEnabled
            AabCov = [dA1(:,1) dA2(:,1)]'*[dA1(:,1) dA2(:,1)];
        end

        %% 4ii.9. Compute the contravariant base vectors
        if isLagrangeMultipliersRotationsEnabled
            AContravariant = (AabCov\[dA1(:,1) dA2(:,1)]')';
        end
        
        %% 4ii.10. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
        
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
                RMatrix(1,3*k-2) = dR(k,1);
                RMatrix(2,3*k-1) = dR(k,1);
                RMatrix(3,3*k) = dR(k,1);

                if isLagrangeMultipliersRotationsEnabled
                    % Matrix containing the derivatives of the basis functions
                    % With respect to xi:
                    DRsDxiMatrix(1,3*k-2) = dR(k,2);
                    DRsDxiMatrix(2,3*k-1) = dR(k,2);
                    DRsDxiMatrix(3,3*k) = dR(k,2);

                    % With respect to eta:
                    DRsDetaMatrix(1,3*k-2) = dR(k,4);
                    DRsDetaMatrix(2,3*k-1) = dR(k,4);
                    DRsDetaMatrix(3,3*k) = dR(k,4);
                end
            end
        end
        
        % For the Lagrange Multipliers field for the traction forces :
        % ____________________________________________________________
        
        for k = 1:noCPsLambdaLoc
            RLambdaMatrix(1,3*k-2) = RLambda(k);
            RLambdaMatrix(2,3*k-1) = RLambda(k);
            RLambdaMatrix(3,3*k) = RLambda(k);
        end
        
        % For the Lagrange Multipliers field for the traction forces :
        % ____________________________________________________________
        if isLagrangeMultipliersRotationsEnabled
            for k = 1:noCPsMuLoc
                RMuMatrix(1,2*k-1) = RMu(k);
                RMuMatrix(2,2*k) = RMu(k);
            end
        end
        
        %% 4ii.11. Transform the normal and the tangent vectors to the covariant bases
        if isLagrangeMultipliersRotationsEnabled
            nCovariant = AContravariant'*n;
            tCovariant = AContravariant'*t;
        end
        
        %% 4ii.12. Compute the curvature coefficients
        if isLagrangeMultipliersRotationsEnabled
            BV = [dA1(:,2) dA2(:,2) dA1(:,3)]'*A3;
        end
        
        %% 4ii.13. Compute the B-operator matrices for the rotations
        if isLagrangeMultipliersRotationsEnabled
            [Bt,Bn,~,~] = computeBOperatorMatrix4RotationsIGAKirchhoffLoveShell(RMatrix,DRsDxiMatrix,DRsDetaMatrix,...
                A3,dA3,AContravariant,BV,nCovariant,tCovariant);
            BRotations = [Bt
                          Bn];
        end

        %% 4ii.14. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
        if isOnXi
            detJxxi = norm(dA1(:,1));
        else
            detJxxi = norm(dA2(:,1));
        end
        
        %% 4ii.15. Compute the element length at the GP
        elementLengthOnGP = detJxxi*detJxizeta*GW(j);
        
        %% 4ii.16. Compute the Lagrange Multipliers matrices and add them to the global matrices
        
        % Lagrange Multipliers matrix for the traction forces :
        % _____________________________________________________
        
        if isMaster
            Lambda(EFT,EFTLambda) = Lambda(EFT,EFTLambda) + RMatrix'*RLambdaMatrix*elementLengthOnGP;
        else
            Lambda(EFT,EFTLambda) = Lambda(EFT,EFTLambda) - RMatrix'*RLambdaMatrix*elementLengthOnGP;
        end
        
        % Lagrange Multipliers matrix for the traction moments :
        % ______________________________________________________
        
        if isLagrangeMultipliersRotationsEnabled
            Mu(EFT,EFTMu) = Mu(EFT,EFTMu) + BRotations'*RMuMatrix*elementLengthOnGP;
        end
    end
end

end
